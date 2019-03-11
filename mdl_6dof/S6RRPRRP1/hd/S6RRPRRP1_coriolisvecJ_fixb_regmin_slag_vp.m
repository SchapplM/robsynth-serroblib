% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:15
% EndTime: 2019-03-09 11:41:24
% DurationCPUTime: 3.32s
% Computational Cost: add. (7578->348), mult. (19527->463), div. (0->0), fcn. (14954->8), ass. (0->207)
t188 = sin(qJ(5));
t191 = cos(qJ(5));
t189 = sin(qJ(4));
t192 = cos(qJ(4));
t237 = qJD(4) * t189;
t186 = sin(pkin(10));
t187 = cos(pkin(10));
t193 = cos(qJ(2));
t231 = qJD(1) * qJD(2);
t223 = t193 * t231;
t190 = sin(qJ(2));
t224 = t190 * t231;
t142 = -t186 * t224 + t187 * t223;
t275 = -qJ(3) - pkin(7);
t219 = qJD(2) * t275;
t147 = t193 * qJD(3) + t190 * t219;
t129 = t147 * qJD(1);
t148 = -t190 * qJD(3) + t193 * t219;
t130 = t148 * qJD(1);
t89 = -t186 * t129 + t187 * t130;
t72 = -t142 * pkin(8) + t89;
t163 = t186 * t193 + t187 * t190;
t153 = t163 * qJD(2);
t141 = qJD(1) * t153;
t92 = t187 * t129 + t186 * t130;
t73 = -t141 * pkin(8) + t92;
t173 = t275 * t193;
t168 = qJD(1) * t173;
t157 = t186 * t168;
t171 = t275 * t190;
t167 = qJD(1) * t171;
t264 = qJD(2) * pkin(2);
t161 = t167 + t264;
t109 = t187 * t161 + t157;
t154 = t163 * qJD(1);
t278 = t154 * pkin(8);
t84 = qJD(2) * pkin(3) + t109 - t278;
t247 = t187 * t168;
t110 = t186 * t161 - t247;
t162 = -t186 * t190 + t187 * t193;
t152 = t162 * qJD(1);
t279 = t152 * pkin(8);
t87 = t110 + t279;
t14 = (qJD(4) * t84 + t73) * t192 + t189 * t72 - t87 * t237;
t230 = qJD(2) + qJD(4);
t43 = t189 * t84 + t192 * t87;
t40 = pkin(9) * t230 + t43;
t107 = t192 * t152 - t189 * t154;
t227 = -t193 * pkin(2) - pkin(1);
t211 = t227 * qJD(1);
t169 = qJD(3) + t211;
t115 = -t152 * pkin(3) + t169;
t207 = t189 * t152 + t192 * t154;
t47 = -pkin(4) * t107 - pkin(9) * t207 + t115;
t20 = t188 * t47 + t191 * t40;
t177 = pkin(2) * t224;
t116 = t141 * pkin(3) + t177;
t236 = qJD(4) * t192;
t199 = -t189 * t141 + t192 * t142 + t152 * t236 - t154 * t237;
t215 = t192 * t141 + t189 * t142;
t288 = qJD(4) * t207;
t63 = t215 + t288;
t27 = t63 * pkin(4) - pkin(9) * t199 + t116;
t26 = t191 * t27;
t196 = -qJD(5) * t20 - t188 * t14 + t26;
t212 = t191 * t230;
t235 = qJD(5) * t188;
t35 = -qJD(5) * t212 - t191 * t199 + t207 * t235;
t93 = t188 * t230 + t191 * t207;
t1 = t63 * pkin(5) + t35 * qJ(6) - t93 * qJD(6) + t196;
t90 = t188 * t207 - t212;
t11 = -t90 * qJ(6) + t20;
t304 = qJD(5) - t107;
t287 = t11 * t304 + t1;
t234 = qJD(5) * t191;
t201 = t191 * t14 + t188 * t27 + t47 * t234 - t235 * t40;
t36 = qJD(5) * t93 + t188 * t199;
t3 = -t36 * qJ(6) - t90 * qJD(6) + t201;
t315 = -t287 * t188 + t3 * t191;
t253 = t304 * t188;
t314 = pkin(5) * t253;
t179 = t187 * pkin(2) + pkin(3);
t280 = pkin(2) * t186;
t291 = t192 * t179 - t189 * t280;
t113 = -t186 * t167 + t247;
t95 = t113 - t279;
t114 = t187 * t167 + t157;
t96 = t114 - t278;
t295 = -t291 * qJD(4) + t189 * t95 + t192 * t96;
t292 = t107 * t230;
t312 = t199 - t292;
t307 = t107 * t191;
t33 = t35 * t188;
t311 = -t33 + (t234 - t307) * t93;
t308 = t107 * t188;
t310 = qJ(6) * t308 + t191 * qJD(6);
t257 = t93 * t207;
t60 = t188 * t63;
t267 = t234 * t304 + t60;
t309 = -t304 * t307 - t257 + t267;
t306 = t207 * t107;
t303 = -t107 ^ 2 + t207 ^ 2;
t67 = pkin(4) * t207 - t107 * pkin(9);
t302 = -t115 * t107 - t14;
t183 = t191 * qJ(6);
t301 = -t207 * pkin(5) + t107 * t183;
t203 = -t191 * t35 - t93 * t235;
t259 = t191 * t90;
t261 = t188 * t93;
t209 = t259 + t261;
t271 = -t188 * t36 - t90 * t234;
t300 = t107 * t209 + t203 + t271;
t298 = -0.2e1 * t231;
t263 = t207 * t90;
t204 = t189 * t179 + t192 * t280;
t256 = t204 * qJD(4) - t189 * t96 + t192 * t95;
t62 = t191 * t63;
t294 = t235 * t304 - t62;
t15 = t189 * t73 - t192 * t72 + t87 * t236 + t84 * t237;
t42 = -t189 * t87 + t192 * t84;
t39 = -pkin(4) * t230 - t42;
t222 = -t15 * t191 + t39 * t235;
t117 = t187 * t171 + t186 * t173;
t101 = -t163 * pkin(8) + t117;
t118 = t186 * t171 - t187 * t173;
t102 = t162 * pkin(8) + t118;
t58 = -t192 * t101 + t189 * t102;
t293 = t304 * t207;
t233 = t207 * qJD(2);
t290 = t233 - t215;
t238 = qJD(1) * t190;
t124 = pkin(2) * t238 + t154 * pkin(3);
t52 = t124 + t67;
t289 = t188 * t52 + t295 * t191;
t286 = -t115 * t207 - t15;
t19 = -t188 * t40 + t191 * t47;
t285 = -t19 * t207 + t222;
t284 = t15 * t188 + t20 * t207 + t39 * t234;
t283 = t93 ^ 2;
t10 = -t93 * qJ(6) + t19;
t7 = pkin(5) * t304 + t10;
t282 = t10 - t7;
t146 = pkin(9) + t204;
t240 = -qJ(6) - t146;
t214 = qJD(5) * t240;
t49 = t191 * t52;
t281 = t191 * t214 - t49 + (-qJD(6) + t295) * t188 + t301;
t277 = t191 * pkin(5);
t276 = t191 * t7;
t274 = -qJ(6) - pkin(9);
t218 = qJD(5) * t274;
t221 = -t188 * t42 + t191 * t67;
t272 = -t188 * qJD(6) + t191 * t218 - t221 + t301;
t270 = t188 * t67 + t191 * t42;
t59 = t189 * t101 + t192 * t102;
t56 = t191 * t59;
t112 = t189 * t162 + t192 * t163;
t132 = -t162 * pkin(3) + t227;
t206 = t192 * t162 - t189 * t163;
t57 = -pkin(4) * t206 - t112 * pkin(9) + t132;
t268 = t188 * t57 + t56;
t266 = t188 * t214 - t289 + t310;
t156 = t162 * qJD(2);
t68 = qJD(4) * t206 - t189 * t153 + t192 * t156;
t260 = t191 * t68;
t258 = t191 * t93;
t255 = t188 * t218 - t270 + t310;
t248 = t112 * t188;
t195 = qJD(1) ^ 2;
t243 = t193 * t195;
t194 = qJD(2) ^ 2;
t242 = t194 * t190;
t241 = t194 * t193;
t100 = t187 * t147 + t186 * t148;
t239 = t190 ^ 2 - t193 ^ 2;
t99 = -t186 * t147 + t187 * t148;
t75 = -t156 * pkin(8) + t99;
t76 = -t153 * pkin(8) + t100;
t23 = -t58 * qJD(4) + t189 * t75 + t192 * t76;
t181 = t190 * t264;
t125 = t153 * pkin(3) + t181;
t69 = qJD(4) * t112 + t192 * t153 + t189 * t156;
t30 = t69 * pkin(4) - t68 * pkin(9) + t125;
t229 = t188 * t30 + t191 * t23 + t57 * t234;
t226 = t112 * t234;
t216 = pkin(1) * t298;
t213 = t304 * t191;
t145 = -pkin(4) - t291;
t210 = -t39 * t107 - t146 * t63;
t208 = -qJ(6) * t68 - qJD(6) * t112;
t205 = t304 * t308 - t294;
t6 = t36 * pkin(5) + t15;
t202 = t188 * t68 + t226;
t24 = qJD(4) * t59 + t189 * t76 - t192 * t75;
t172 = t191 * pkin(9) + t183;
t170 = t274 * t188;
t120 = t191 * t146 + t183;
t119 = t240 * t188;
t88 = t90 ^ 2;
t55 = t191 * t57;
t31 = t90 * pkin(5) + qJD(6) + t39;
t29 = t191 * t30;
t21 = -qJ(6) * t248 + t268;
t17 = -pkin(5) * t206 - t112 * t183 - t188 * t59 + t55;
t5 = -qJ(6) * t226 + (-qJD(5) * t59 + t208) * t188 + t229;
t4 = t69 * pkin(5) - t188 * t23 + t29 + t208 * t191 + (-t56 + (qJ(6) * t112 - t57) * t188) * qJD(5);
t2 = [0, 0, 0, 0.2e1 * t190 * t223, t239 * t298, t241, -t242, 0, -pkin(7) * t241 + t190 * t216, pkin(7) * t242 + t193 * t216, t100 * t152 - t109 * t156 - t110 * t153 - t117 * t142 - t118 * t141 - t99 * t154 + t92 * t162 - t89 * t163, t110 * t100 + t109 * t99 + t89 * t117 + t92 * t118 + (t169 + t211) * t181, t112 * t199 + t207 * t68, t107 * t68 - t112 * t63 + t199 * t206 - t207 * t69, t68 * t230, -t69 * t230, 0, -t107 * t125 + t115 * t69 - t116 * t206 + t132 * t63 - t230 * t24, t116 * t112 + t115 * t68 + t125 * t207 + t132 * t199 - t23 * t230, t112 * t203 + t258 * t68, -t209 * t68 + (t33 - t191 * t36 + (t188 * t90 - t258) * qJD(5)) * t112, t112 * t62 + t35 * t206 + t93 * t69 + (-t112 * t235 + t260) * t304, -t202 * t304 + t206 * t36 - t248 * t63 - t90 * t69, -t206 * t63 + t304 * t69 (-t234 * t59 + t29) * t304 + t55 * t63 - (-t234 * t40 + t26) * t206 + t19 * t69 + t24 * t90 + t58 * t36 + t39 * t226 + ((-qJD(5) * t57 - t23) * t304 - t59 * t63 - (-qJD(5) * t47 - t14) * t206 + t15 * t112 + t39 * t68) * t188 -(-t235 * t59 + t229) * t304 - t268 * t63 + t201 * t206 - t20 * t69 + t24 * t93 - t58 * t35 + t39 * t260 - t222 * t112, t17 * t35 - t21 * t36 - t4 * t93 - t5 * t90 + (-t11 * t188 - t276) * t68 + (-t1 * t191 - t3 * t188 + (-t11 * t191 + t188 * t7) * qJD(5)) * t112, t3 * t21 + t11 * t5 + t1 * t17 + t7 * t4 + t6 * (pkin(5) * t248 + t58) + t31 * (pkin(5) * t202 + t24); 0, 0, 0, -t190 * t243, t239 * t195, 0, 0, 0, t195 * pkin(1) * t190, pkin(1) * t243 (t110 + t113) * t154 + (t109 - t114) * t152 + (-t141 * t186 - t142 * t187) * pkin(2), -t109 * t113 - t110 * t114 + (-t169 * t238 + t186 * t92 + t187 * t89) * pkin(2), -t306, t303, t312, t290, 0, t107 * t124 - t256 * t230 + t286, -t124 * t207 + t295 * t230 + t302, t311, t300, t309, t205 + t263, -t293, t145 * t36 + t256 * t90 + t210 * t188 + (-t146 * t234 + t295 * t188 - t49) * t304 + t285, -t145 * t35 + t256 * t93 + t210 * t191 + (t146 * t235 + t289) * t304 + t284, t119 * t35 - t120 * t36 - t7 * t213 - t266 * t90 - t281 * t93 + t315, t3 * t120 + t1 * t119 + t6 * (t145 - t277) + t281 * t7 + (t256 + t314) * t31 + t266 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 ^ 2 - t154 ^ 2, t109 * t154 - t110 * t152 + t177, 0, 0, 0, 0, 0, t215 + t233 + 0.2e1 * t288, t199 + t292, 0, 0, 0, 0, 0, t205 - t263, -t213 * t304 - t257 - t60 (t259 - t261) * t107 - t203 + t271, -t31 * t207 + t287 * t191 + (-t304 * t7 + t3) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, t303, t312, t290, 0, t230 * t43 + t286, t230 * t42 + t302, t311, t300, t309, -t253 * t304 + t263 + t62, -t293, -pkin(4) * t36 - t267 * pkin(9) - t221 * t304 - t308 * t39 - t43 * t90 + t285, pkin(4) * t35 + t294 * pkin(9) + t270 * t304 - t307 * t39 - t43 * t93 + t284, t170 * t35 - t172 * t36 - t255 * t90 - t272 * t93 - t304 * t276 + t315, t3 * t172 + t1 * t170 + t6 * (-pkin(4) - t277) + t272 * t7 + (-t43 + t314) * t31 + t255 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t90, -t88 + t283, t304 * t90 - t35, t304 * t93 - t36, t63, t20 * t304 - t39 * t93 + t196, t19 * t304 + t39 * t90 - t201, pkin(5) * t35 + t282 * t90, -t282 * t11 + (-t31 * t93 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 - t283, t11 * t90 + t7 * t93 + t6;];
tauc_reg  = t2;
