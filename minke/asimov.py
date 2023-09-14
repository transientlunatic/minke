import asimov.pipeline
import importlib
import os
import glob
import htcondor
from asimov.utils import set_directory
from asimov import config
import configparser

class Asimov(asimov.pipeline.Pipeline):

    name = 'minke'
    with importlib.resources.path(f"minke", f"{name}_template.yml") as template_file:
        config_template = template_file
    _pipeline_command = "minke"

    def build_dag(self, dryrun=False):
        """
        Create a condor submission description.
        """
        name = self.production.name
        ini = self.production.event.repository.find_prods(name, self.category)[0]
        executable = os.path.join(config.get('pipelines', 'environment'), 'bin', self._pipeline_command)
        command = ["--settings", ini]
        full_command = executable + " " + " ".join(command)
        self.logger.info(full_command)

        description = {
            "executable": f"{executable}",
            "arguments": f"{' '.join(command)}",
            "output": f"{name}.out",
            "error": f"{name}.err",
            "log": f"{name}.log",
            "request_disk": "1024",
            "request_memory": "1024",
            "batch_name": f"{self.name}/{self.production.event.name}/{name}",
            "+flock_local": "True",
            "+DESIRED_Sites": htcondor.classad.quote("nogrid"),
        }

        accounting_group = self.production.meta.get("scheduler", {}).get("accounting group", None)

        if accounting_group:
            description["accounting_group_user"] = config.get('condor', 'user')
            description["accounting_group"] = accounting_group
        else:
            self.logger.warning(
                "This job does not supply any accounting information, which may prevent it running on some clusters."
            )

        
        job = htcondor.Submit(description)
        os.makedirs(self.production.rundir, exist_ok=True)
        with set_directory(self.production.rundir):
            with open(f"{name}.sub", "w") as subfile:
                subfile.write(job.__str__())

            full_command = f"""#! /bin/bash
{ full_command }
"""

            with open(f"{name}.sh", "w") as bashfile:
                bashfile.write(str(full_command))

        with set_directory(self.production.rundir):
            try:
                schedulers = htcondor.Collector().locate(
                    htcondor.DaemonTypes.Schedd, config.get("condor", "scheduler")
                )
            except configparser.NoOptionError:
                schedulers = htcondor.Collector().locate(htcondor.DaemonTypes.Schedd)
            schedd = htcondor.Schedd(schedulers)
            with schedd.transaction() as txn:
                cluster_id = job.queue(txn)
                self.logger.info("Submitted to htcondor job queue.")

        self.production.job_id = int(cluster_id)
        self.clusterid = cluster_id

    def submit_dag(self, dryrun=False):
        self.production.status = "running"
        self.production.job_id = int(self.clusterid)
        return self.clusterid

    def detect_completion(self):
        self.logger.info("Checking for completion.")
        assets = self.collect_assets()
        if len(list(assets.keys())) > 0:
            self.logger.info("Outputs detected, job complete.")
            return True
        else:
            self.logger.info(f"{self.name} job completion was not detected.")
            return False

    def after_completion(self):
        self.production.status = "uploaded"
        self.production.event.update_data()

    def collect_assets(self):
        """
        Collect the assets for this job.
        """
        outputs = {}
        if os.path.exists(os.path.join(self.production.rundir, "frames")):
            results_dir = glob.glob(os.path.join(self.production.rundir, "frames", "*"))
            frames = {}

            for frame in results_dir:
                ifo = frame.split("/")[-1].split("-")[0]
                frames[ifo] = frame

            outputs["frames"] = frames

            self.production.event.meta['data']['data files'] = frames

        if os.path.exists(os.path.join(self.production.rundir, "cache")):
            results_dir = glob.glob(os.path.join(self.production.rundir, "cache", "*"))
            frames = {}

            for frame in results_dir:
                ifo = frame.split("/")[-1].split(".")[0]
                frames[ifo] = frame

            outputs["cache"] = frames

            self.production.event.meta['data']['cache files'] = frames

        if os.path.exists(os.path.join(self.production.rundir, "psds")):
            results_dir = glob.glob(os.path.join(self.production.rundir, "psds", "*"))
            frames = {}

            for frame in results_dir:
                ifo = frame.split("/")[-1].split("-")[0]
                frames[ifo] = frame

            outputs["psds"] = frames

            
        return outputs

    def html(self):
        """Return the HTML representation of this pipeline."""
        out = ""
        if self.production.status in {"finished", "uploaded"}:
            out += """<div class="asimov-pipeline">"""
            pp = pprint.PrettyPrinter(indent=4)
            out += f"<pre>{ pp.pprint(self.collect_assets()) }</pre>"
            out += """</div>"""

        return out


