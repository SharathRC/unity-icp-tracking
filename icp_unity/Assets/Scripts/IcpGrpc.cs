// <auto-generated>
//     Generated by the protocol buffer compiler.  DO NOT EDIT!
//     source: icp.proto
// </auto-generated>
#pragma warning disable 0414, 1591
#region Designer generated code

using grpc = global::Grpc.Core;

namespace Protoicp {
  /// <summary>
  ///Services define the different communication scenarios available
  ///A simple service for generating color
  /// </summary>
  public static partial class CalcTransform
  {
    static readonly string __ServiceName = "protoicp.CalcTransform";

    static readonly grpc::Marshaller<global::Protoicp.ObjectClouds> __Marshaller_protoicp_ObjectClouds = grpc::Marshallers.Create((arg) => global::Google.Protobuf.MessageExtensions.ToByteArray(arg), global::Protoicp.ObjectClouds.Parser.ParseFrom);
    static readonly grpc::Marshaller<global::Protoicp.Transform> __Marshaller_protoicp_Transform = grpc::Marshallers.Create((arg) => global::Google.Protobuf.MessageExtensions.ToByteArray(arg), global::Protoicp.Transform.Parser.ParseFrom);

    static readonly grpc::Method<global::Protoicp.ObjectClouds, global::Protoicp.Transform> __Method_getTransform = new grpc::Method<global::Protoicp.ObjectClouds, global::Protoicp.Transform>(
        grpc::MethodType.Unary,
        __ServiceName,
        "getTransform",
        __Marshaller_protoicp_ObjectClouds,
        __Marshaller_protoicp_Transform);

    /// <summary>Service descriptor</summary>
    public static global::Google.Protobuf.Reflection.ServiceDescriptor Descriptor
    {
      get { return global::Protoicp.IcpReflection.Descriptor.Services[0]; }
    }

    /// <summary>Base class for server-side implementations of CalcTransform</summary>
    [grpc::BindServiceMethod(typeof(CalcTransform), "BindService")]
    public abstract partial class CalcTransformBase
    {
      /// <summary>
      ///Generates a random color on each request
      /// </summary>
      /// <param name="request">The request received from the client.</param>
      /// <param name="context">The context of the server-side call handler being invoked.</param>
      /// <returns>The response to send back to the client (wrapped by a task).</returns>
      public virtual global::System.Threading.Tasks.Task<global::Protoicp.Transform> getTransform(global::Protoicp.ObjectClouds request, grpc::ServerCallContext context)
      {
        throw new grpc::RpcException(new grpc::Status(grpc::StatusCode.Unimplemented, ""));
      }

    }

    /// <summary>Client for CalcTransform</summary>
    public partial class CalcTransformClient : grpc::ClientBase<CalcTransformClient>
    {
      /// <summary>Creates a new client for CalcTransform</summary>
      /// <param name="channel">The channel to use to make remote calls.</param>
      public CalcTransformClient(grpc::ChannelBase channel) : base(channel)
      {
      }
      /// <summary>Creates a new client for CalcTransform that uses a custom <c>CallInvoker</c>.</summary>
      /// <param name="callInvoker">The callInvoker to use to make remote calls.</param>
      public CalcTransformClient(grpc::CallInvoker callInvoker) : base(callInvoker)
      {
      }
      /// <summary>Protected parameterless constructor to allow creation of test doubles.</summary>
      protected CalcTransformClient() : base()
      {
      }
      /// <summary>Protected constructor to allow creation of configured clients.</summary>
      /// <param name="configuration">The client configuration.</param>
      protected CalcTransformClient(ClientBaseConfiguration configuration) : base(configuration)
      {
      }

      /// <summary>
      ///Generates a random color on each request
      /// </summary>
      /// <param name="request">The request to send to the server.</param>
      /// <param name="headers">The initial metadata to send with the call. This parameter is optional.</param>
      /// <param name="deadline">An optional deadline for the call. The call will be cancelled if deadline is hit.</param>
      /// <param name="cancellationToken">An optional token for canceling the call.</param>
      /// <returns>The response received from the server.</returns>
      public virtual global::Protoicp.Transform getTransform(global::Protoicp.ObjectClouds request, grpc::Metadata headers = null, global::System.DateTime? deadline = null, global::System.Threading.CancellationToken cancellationToken = default(global::System.Threading.CancellationToken))
      {
        return getTransform(request, new grpc::CallOptions(headers, deadline, cancellationToken));
      }
      /// <summary>
      ///Generates a random color on each request
      /// </summary>
      /// <param name="request">The request to send to the server.</param>
      /// <param name="options">The options for the call.</param>
      /// <returns>The response received from the server.</returns>
      public virtual global::Protoicp.Transform getTransform(global::Protoicp.ObjectClouds request, grpc::CallOptions options)
      {
        return CallInvoker.BlockingUnaryCall(__Method_getTransform, null, options, request);
      }
      /// <summary>
      ///Generates a random color on each request
      /// </summary>
      /// <param name="request">The request to send to the server.</param>
      /// <param name="headers">The initial metadata to send with the call. This parameter is optional.</param>
      /// <param name="deadline">An optional deadline for the call. The call will be cancelled if deadline is hit.</param>
      /// <param name="cancellationToken">An optional token for canceling the call.</param>
      /// <returns>The call object.</returns>
      public virtual grpc::AsyncUnaryCall<global::Protoicp.Transform> getTransformAsync(global::Protoicp.ObjectClouds request, grpc::Metadata headers = null, global::System.DateTime? deadline = null, global::System.Threading.CancellationToken cancellationToken = default(global::System.Threading.CancellationToken))
      {
        return getTransformAsync(request, new grpc::CallOptions(headers, deadline, cancellationToken));
      }
      /// <summary>
      ///Generates a random color on each request
      /// </summary>
      /// <param name="request">The request to send to the server.</param>
      /// <param name="options">The options for the call.</param>
      /// <returns>The call object.</returns>
      public virtual grpc::AsyncUnaryCall<global::Protoicp.Transform> getTransformAsync(global::Protoicp.ObjectClouds request, grpc::CallOptions options)
      {
        return CallInvoker.AsyncUnaryCall(__Method_getTransform, null, options, request);
      }
      /// <summary>Creates a new instance of client from given <c>ClientBaseConfiguration</c>.</summary>
      protected override CalcTransformClient NewInstance(ClientBaseConfiguration configuration)
      {
        return new CalcTransformClient(configuration);
      }
    }

    /// <summary>Creates service definition that can be registered with a server</summary>
    /// <param name="serviceImpl">An object implementing the server-side handling logic.</param>
    public static grpc::ServerServiceDefinition BindService(CalcTransformBase serviceImpl)
    {
      return grpc::ServerServiceDefinition.CreateBuilder()
          .AddMethod(__Method_getTransform, serviceImpl.getTransform).Build();
    }

    /// <summary>Register service method with a service binder with or without implementation. Useful when customizing the  service binding logic.
    /// Note: this method is part of an experimental API that can change or be removed without any prior notice.</summary>
    /// <param name="serviceBinder">Service methods will be bound by calling <c>AddMethod</c> on this object.</param>
    /// <param name="serviceImpl">An object implementing the server-side handling logic.</param>
    public static void BindService(grpc::ServiceBinderBase serviceBinder, CalcTransformBase serviceImpl)
    {
      serviceBinder.AddMethod(__Method_getTransform, serviceImpl == null ? null : new grpc::UnaryServerMethod<global::Protoicp.ObjectClouds, global::Protoicp.Transform>(serviceImpl.getTransform));
    }

  }
}
#endregion
