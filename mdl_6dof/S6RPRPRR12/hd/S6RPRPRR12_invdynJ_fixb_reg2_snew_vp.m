% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 20:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:46:37
% EndTime: 2019-05-05 20:46:46
% DurationCPUTime: 3.08s
% Computational Cost: add. (12747->340), mult. (25817->421), div. (0->0), fcn. (15166->8), ass. (0->219)
t262 = pkin(1) + pkin(7);
t181 = sin(qJ(6));
t183 = sin(qJ(3));
t234 = qJD(1) * qJD(3);
t170 = t183 * t234;
t187 = cos(qJ(3));
t172 = t187 * qJDD(1);
t149 = t172 - t170;
t140 = qJDD(5) + t149;
t136 = qJDD(6) + t140;
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t236 = qJD(1) * t183;
t141 = qJD(3) * t182 - t186 * t236;
t143 = qJD(3) * t186 + t182 * t236;
t185 = cos(qJ(6));
t112 = t185 * t141 + t143 * t181;
t114 = -t141 * t181 + t143 * t185;
t85 = t114 * t112;
t266 = -t85 + t136;
t274 = t181 * t266;
t120 = t143 * t141;
t265 = -t120 + t140;
t273 = t182 * t265;
t272 = t185 * t266;
t271 = t186 * t265;
t179 = t183 ^ 2;
t190 = qJD(1) ^ 2;
t173 = t179 * t190;
t189 = qJD(3) ^ 2;
t160 = -t173 - t189;
t237 = t187 * t190;
t229 = t183 * t237;
t155 = -qJDD(3) + t229;
t244 = t155 * t187;
t123 = -t160 * t183 + t244;
t270 = t262 * t123;
t180 = t187 ^ 2;
t174 = t180 * t190;
t162 = -t174 - t189;
t156 = qJDD(3) + t229;
t243 = t156 * t183;
t124 = -t162 * t187 + t243;
t269 = t262 * t124;
t268 = t149 - t170;
t224 = t187 * t234;
t231 = t183 * qJDD(1);
t148 = t224 + t231;
t178 = qJDD(1) * qJ(2);
t184 = sin(qJ(1));
t188 = cos(qJ(1));
t218 = t188 * g(1) + t184 * g(2);
t211 = -t178 + t218;
t199 = -pkin(3) * t224 + qJ(4) * t268 + t211;
t230 = t262 * t190;
t267 = t148 * pkin(3) - t199 - t230;
t107 = -t141 * qJD(5) + t186 * qJDD(3) + t182 * t148;
t235 = t187 * qJD(1);
t167 = qJD(5) + t235;
t245 = t141 * t167;
t91 = -t107 - t245;
t233 = qJD(2) * qJD(1);
t176 = 0.2e1 * t233;
t264 = -0.2e1 * qJD(4) * t235 + t176;
t238 = t187 * qJ(4);
t216 = t183 * pkin(3) - t238;
t144 = t216 * qJD(1);
t263 = qJDD(3) * qJ(4) - t144 * t236;
t110 = t112 ^ 2;
t111 = t114 ^ 2;
t138 = t141 ^ 2;
t139 = t143 ^ 2;
t158 = qJD(6) + t167;
t157 = t158 ^ 2;
t163 = t167 ^ 2;
t261 = pkin(3) + pkin(8);
t260 = pkin(3) * t189;
t207 = pkin(4) * t235 - qJD(3) * pkin(8);
t68 = -t207 * t235 + t261 * t148 + (-t179 * pkin(4) - t262) * t190 - t199 + t264;
t217 = t184 * g(1) - t188 * g(2);
t205 = qJDD(2) - t217;
t200 = -t190 * qJ(2) + t205;
t197 = -t262 * qJDD(1) + t200;
t195 = t187 * t197;
t194 = t144 * t235 + qJDD(4) - t195;
t192 = t189 * qJ(4) - t194;
t74 = t149 * pkin(4) - t261 * qJDD(3) + (pkin(4) * t234 + pkin(8) * t237 - g(3)) * t183 - t192;
t37 = t182 * t68 - t186 * t74;
t33 = t265 * pkin(5) + t91 * pkin(9) - t37;
t222 = qJDD(3) * t182 - t186 * t148;
t106 = -qJD(5) * t143 - t222;
t219 = pkin(5) * t167 - pkin(9) * t143;
t38 = t182 * t74 + t186 * t68;
t36 = -t138 * pkin(5) + t106 * pkin(9) - t167 * t219 + t38;
t14 = t181 * t36 - t185 * t33;
t15 = t181 * t33 + t185 * t36;
t8 = -t14 * t185 + t15 * t181;
t259 = t182 * t8;
t258 = t183 * g(3);
t257 = t186 * t8;
t127 = t183 * t197;
t117 = g(3) * t187 - t127;
t202 = -t117 + t263;
t198 = -t202 + t260;
t232 = qJD(4) * qJD(3);
t73 = pkin(4) * t148 + pkin(8) * t173 - qJD(3) * t207 + t198 - 0.2e1 * t232;
t39 = pkin(5) * t106 + pkin(9) * t138 - t143 * t219 + t73;
t256 = t181 * t39;
t76 = t85 + t136;
t255 = t181 * t76;
t254 = t182 * t73;
t253 = t185 * t39;
t252 = t185 * t76;
t251 = t186 * t73;
t201 = (-qJD(5) + t167) * t143 - t222;
t63 = t182 * t201 + t186 * t91;
t250 = t187 * t63;
t249 = qJDD(1) * pkin(1);
t100 = t120 + t140;
t248 = t100 * t182;
t247 = t100 * t186;
t246 = t112 * t158;
t242 = t158 * t181;
t241 = t158 * t185;
t240 = t167 * t182;
t239 = t167 * t186;
t228 = t187 * t85;
t227 = t187 * t120;
t226 = t112 * qJD(6) - t181 * t106 - t185 * t107;
t9 = t14 * t181 + t185 * t15;
t223 = -t185 * t106 + t107 * t181;
t80 = -t157 - t110;
t42 = t181 * t80 + t272;
t221 = pkin(5) * t42 - t14;
t152 = (t179 + t180) * qJDD(1);
t154 = -t174 - t173;
t220 = qJ(2) * t154 + t152 * t262;
t21 = t182 * t38 - t186 * t37;
t215 = t182 * t37 + t186 * t38;
t214 = t183 * t73 + t187 * t21;
t116 = t195 + t258;
t86 = t116 * t187 - t117 * t183;
t213 = (-t174 + t189) * t183 + t244;
t212 = -(t173 - t189) * t187 + t243;
t210 = -t226 - t246;
t209 = t107 - t245;
t95 = -t111 - t157;
t56 = t185 * t95 - t255;
t208 = pkin(5) * t56 - t15;
t206 = qJ(2) + t216;
t204 = t183 * t261 + qJ(2) - t238;
t203 = (-qJD(6) + t158) * t114 - t223;
t193 = t264 + t267;
t191 = -qJDD(3) * pkin(3) - t192;
t175 = 0.2e1 * t232;
t153 = t174 - t173;
t150 = t172 - 0.2e1 * t170;
t147 = -0.2e1 * t224 - t231;
t131 = -t200 + t249;
t130 = -t139 + t163;
t129 = t138 - t163;
t128 = t211 + t230 - 0.2e1 * t233;
t126 = t268 * t187;
t121 = (t148 + t224) * t183;
t119 = t139 - t138;
t118 = -t139 - t163;
t115 = t187 * t147 - t183 * t150;
t109 = -t163 - t138;
t98 = -t138 - t139;
t97 = -t111 + t157;
t96 = t110 - t157;
t94 = -t191 + t258;
t93 = t175 - t198;
t87 = (qJD(5) + t167) * t143 + t222;
t84 = t111 - t110;
t81 = t118 * t186 - t248;
t78 = t109 * t182 + t271;
t71 = (-t112 * t185 + t114 * t181) * t158;
t70 = (-t112 * t181 - t114 * t185) * t158;
t69 = -t110 - t111;
t66 = -qJD(6) * t114 - t223;
t65 = t183 * t93 + t187 * t94;
t62 = t185 * t96 - t255;
t61 = -t181 * t97 + t272;
t60 = t181 * t96 + t252;
t59 = t185 * t97 + t274;
t58 = t183 * t209 - t187 * t81;
t57 = -t181 * t95 - t252;
t54 = t183 * t87 - t187 * t78;
t52 = t226 - t246;
t48 = (qJD(6) + t158) * t114 + t223;
t47 = -t114 * t242 - t185 * t226;
t46 = t114 * t241 - t181 * t226;
t45 = t112 * t241 - t181 * t66;
t44 = t112 * t242 + t185 * t66;
t43 = t185 * t80 - t274;
t40 = t183 * t98 - t250;
t34 = t182 * t57 + t186 * t56;
t31 = -t181 * t52 + t185 * t203;
t30 = -t181 * t210 - t185 * t48;
t29 = t181 * t203 + t185 * t52;
t28 = -t181 * t48 + t185 * t210;
t27 = pkin(5) * t29;
t26 = -pkin(9) * t56 - t253;
t24 = t182 * t43 + t186 * t42;
t23 = -pkin(9) * t42 - t256;
t20 = t183 * t210 - t187 * t34;
t19 = -pkin(5) * t210 + pkin(9) * t57 - t256;
t17 = t183 * t48 - t187 * t24;
t16 = -pkin(5) * t48 + pkin(9) * t43 + t253;
t11 = t182 * t31 + t186 * t29;
t10 = -t187 * t11 + t183 * t69;
t7 = pkin(5) * t8;
t6 = pkin(5) * t39 + pkin(9) * t9;
t5 = -pkin(9) * t29 - t8;
t4 = -pkin(5) * t69 + pkin(9) * t31 + t9;
t2 = t182 * t9 + t257;
t1 = -t183 * t39 - t187 * t2;
t3 = [0, 0, 0, 0, 0, qJDD(1), t217, t218, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t205 - 0.2e1 * t249, t176 + 0.2e1 * t178 - t218, pkin(1) * t131 + qJ(2) * (-t190 * pkin(1) + t176 - t211), t126, t115, -t213, t121, -t212, 0, -qJ(2) * t147 - t128 * t183 + t270, qJ(2) * t150 - t128 * t187 + t269, t220 - t86, -qJ(2) * t128 - t262 * t86, 0, t213, t212, t126, t115, t121, t187 * (-qJ(4) * t154 + t191) - t183 * (-pkin(3) * t154 + t127 + t175 - t260 + t263) + t220, t206 * t147 - t183 * t193 - t270, t187 * (0.2e1 * (qJD(4) * t187 - qJD(2)) * qJD(1) - t267) - t206 * t150 - t269, t206 * t193 - t262 * t65, t227 - t183 * (-t107 * t182 - t143 * t239), t187 * t119 - t183 * (t182 * t87 - t186 * t209), -t187 * t91 - t183 * (-t130 * t186 - t273), -t227 - t183 * (-t106 * t186 - t141 * t240), t187 * t201 - t183 * (-t129 * t182 - t247), t187 * t140 - t183 * (t141 * t182 + t143 * t186) * t167, t187 * (pkin(4) * t78 - t37) - t183 * (pkin(4) * t87 - t251) + t204 * (t109 * t186 - t273) - t262 * t54, t187 * (pkin(4) * t81 - t38) - t183 * (pkin(4) * t209 + t254) + t204 * (-t118 * t182 - t247) - t262 * t58, pkin(4) * t250 - t183 * (pkin(4) * t98 - t215) + t204 * (-t182 * t91 + t186 * t201) - t262 * t40, t204 * t215 + (pkin(4) + t262) * t214, t228 - t183 * (-t182 * t47 - t186 * t46), t187 * t84 - t183 * (-t182 * t30 - t186 * t28), -t187 * t52 - t183 * (-t182 * t61 - t186 * t59), -t228 - t183 * (-t182 * t45 - t186 * t44), t187 * t203 - t183 * (-t182 * t62 - t186 * t60), t187 * t136 - t183 * (-t182 * t71 - t186 * t70), t187 * (pkin(4) * t24 + t221) - t183 * (pkin(4) * t48 - t16 * t186 - t182 * t23) + t204 * (-t182 * t42 + t186 * t43) - t262 * t17, t187 * (pkin(4) * t34 + t208) - t183 * (pkin(4) * t210 - t182 * t26 - t186 * t19) + t204 * (-t182 * t56 + t186 * t57) - t262 * t20, t187 * (pkin(4) * t11 + t27) - t183 * (pkin(4) * t69 - t182 * t5 - t186 * t4) + t204 * (-t182 * t29 + t186 * t31) - t262 * t10, t187 * (pkin(4) * t2 + t7) - t183 * (-pkin(4) * t39 + pkin(9) * t259 - t186 * t6) + t204 * (t186 * t9 - t259) - t262 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t190, -t131, 0, 0, 0, 0, 0, 0, -t123, -t124, -t152, t86, 0, 0, 0, 0, 0, 0, -t152, t123, t124, t65, 0, 0, 0, 0, 0, 0, t54, t58, t40, -t214, 0, 0, 0, 0, 0, 0, t17, t20, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t153, t172, -t229, -t231, qJDD(3), t116, t117, 0, 0, qJDD(3), -t172, t231, t229, t153, -t229, (-pkin(3) * t187 - qJ(4) * t183) * qJDD(1), -t258 + (-t189 - t160) * qJ(4) + (-qJDD(3) + t155) * pkin(3) + t194, qJ(4) * t156 + t175 + (-t162 - t189) * pkin(3) + t202, pkin(3) * t94 + qJ(4) * t93, t107 * t186 - t143 * t240, -t182 * t209 - t186 * t87, -t130 * t182 + t271, -t106 * t182 + t141 * t239, t129 * t186 - t248, (-t141 * t186 + t143 * t182) * t167, qJ(4) * t87 - t261 * t78 - t254, qJ(4) * t209 - t261 * t81 - t251, qJ(4) * t98 - t261 * t63 - t21, -qJ(4) * t73 - t261 * t21, -t182 * t46 + t186 * t47, -t182 * t28 + t186 * t30, -t182 * t59 + t186 * t61, -t182 * t44 + t186 * t45, -t182 * t60 + t186 * t62, -t182 * t70 + t186 * t71, qJ(4) * t48 - t16 * t182 + t186 * t23 - t261 * t24, qJ(4) * t210 - t182 * t19 + t186 * t26 - t261 * t34, qJ(4) * t69 - t261 * t11 - t182 * t4 + t186 * t5, -pkin(9) * t257 - qJ(4) * t39 - t182 * t6 - t261 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t155, t162, -t94, 0, 0, 0, 0, 0, 0, t78, t81, t63, t21, 0, 0, 0, 0, 0, 0, t24, t34, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, -t91, -t120, t201, t140, -t37, -t38, 0, 0, t85, t84, -t52, -t85, t203, t136, t221, t208, t27, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t84, -t52, -t85, t203, t136, -t14, -t15, 0, 0;];
tauJ_reg  = t3;
