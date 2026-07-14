/// <reference types="vite/client" />

declare module "molstar" {
  export * from "molstar/lib/mol-plugin/spec";
}

declare module "molstar/lib/mol-plugin/ui" {
  export function createPluginAsync(container: HTMLElement, options?: any): Promise<any>;
}

declare module "molstar/lib/mol-plugin/spec" {
  export const PluginSpec: any;
}

declare module "molstar/lib/mol-plugin/config" {
  export const PluginConfig: any;
}
