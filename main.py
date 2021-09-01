from logger_factory import LoggerFactory
import RE
import user_IO


logger = LoggerFactory.create_logger("main")


def run_modules():
    modules = [user_IO.UserInputModule, RE.REModule, user_IO.UserOutputModule]

    initial_input = None
    module_input = initial_input
    for module in modules:
        out = module.run_module(module_input)
        # logger.info("Output of module %s is: %s", module.get_name(), out)
        module_input = out

    final_output = out
    logger.info("Final output: %s", final_output)


if __name__ == "__main__":
    run_modules()
