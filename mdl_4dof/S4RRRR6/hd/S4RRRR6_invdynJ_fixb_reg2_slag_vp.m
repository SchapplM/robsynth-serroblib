% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:54
% EndTime: 2019-12-31 17:31:04
% DurationCPUTime: 4.89s
% Computational Cost: add. (5079->475), mult. (13333->672), div. (0->0), fcn. (10184->10), ass. (0->217)
t164 = cos(qJ(3));
t301 = pkin(7) * t164;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t158 = cos(pkin(4));
t249 = qJD(1) * t158;
t208 = qJD(2) + t249;
t165 = cos(qJ(2));
t161 = sin(qJ(2));
t157 = sin(pkin(4));
t250 = qJD(1) * t157;
t225 = t161 * t250;
t233 = pkin(1) * t249;
t97 = -pkin(6) * t225 + t165 * t233;
t74 = -pkin(2) * t208 - t97;
t160 = sin(qJ(3));
t203 = t160 * t225;
t85 = -t164 * t208 + t203;
t87 = t160 * t208 + t164 * t225;
t28 = t85 * pkin(3) - t87 * pkin(8) + t74;
t143 = t165 * t250;
t192 = t143 - qJD(3);
t259 = t157 * t165;
t149 = pkin(6) * t259;
t284 = pkin(1) * t161;
t75 = qJD(2) * pkin(7) + (t149 + (pkin(7) + t284) * t158) * qJD(1);
t187 = -pkin(2) * t165 - pkin(7) * t161 - pkin(1);
t79 = t187 * t250;
t34 = t160 * t79 + t164 * t75;
t30 = -pkin(8) * t192 + t34;
t190 = t159 * t30 - t163 * t28;
t239 = qJD(1) * qJD(2);
t220 = t161 * t239;
t201 = t157 * t220;
t237 = qJDD(1) * t165;
t142 = t157 * t237;
t235 = qJDD(3) - t142;
t174 = t201 + t235;
t244 = qJD(3) * t164;
t246 = qJD(3) * t160;
t238 = qJDD(1) * t158;
t144 = qJDD(2) + t238;
t205 = qJD(2) * t233;
t232 = pkin(1) * t238;
t229 = -pkin(6) * t142 - t161 * t232 - t165 * t205;
t61 = -pkin(6) * t201 - t229;
t47 = pkin(7) * t144 + t61;
t198 = pkin(2) * t161 - pkin(7) * t165;
t181 = t198 * qJD(2);
t54 = (qJD(1) * t181 + qJDD(1) * t187) * t157;
t182 = -t160 * t54 - t164 * t47 - t79 * t244 + t246 * t75;
t5 = pkin(8) * t174 - t182;
t185 = qJD(3) * t208;
t236 = t161 * qJDD(1);
t218 = t157 * t236;
t219 = t165 * t239;
t298 = -t157 * t219 - t218;
t35 = qJD(3) * t203 - t160 * t144 + (-t185 + t298) * t164;
t240 = t160 * qJD(2);
t223 = t165 * t240;
t36 = -t164 * t144 + (qJD(1) * (t161 * t244 + t223) + t160 * t236) * t157 + t160 * t185;
t62 = t298 * pkin(6) - t161 * t205 + t165 * t232;
t48 = -pkin(2) * t144 - t62;
t8 = pkin(3) * t36 + pkin(8) * t35 + t48;
t1 = -t190 * qJD(4) + t159 * t8 + t163 * t5;
t80 = qJD(4) + t85;
t300 = t190 * t80 + t1;
t285 = cos(qJ(1));
t226 = t285 * t165;
t162 = sin(qJ(1));
t257 = t161 * t162;
t109 = -t158 * t226 + t257;
t227 = t285 * t161;
t256 = t162 * t165;
t111 = t158 * t256 + t227;
t171 = -g(1) * t111 - g(2) * t109 + g(3) * t259;
t209 = qJD(3) * t192;
t299 = pkin(7) * t209 - t171 - t48;
t33 = -t160 * t75 + t164 * t79;
t297 = t192 * t33 - t182;
t262 = t157 * t161;
t145 = pkin(6) * t262;
t283 = pkin(1) * t165;
t113 = t158 * t283 - t145;
t101 = qJD(2) * t113;
t110 = t158 * t227 + t256;
t228 = t157 * t285;
t68 = t110 * t164 - t160 * t228;
t296 = -t109 * t163 + t159 * t68;
t295 = t109 * t159 + t163 * t68;
t151 = t158 * t284;
t196 = pkin(3) * t160 - pkin(8) * t164;
t294 = -qJD(4) * t301 + t196 * qJD(3) - (t151 + (pkin(6) + t196) * t259) * qJD(1);
t154 = t157 ^ 2;
t293 = 0.2e1 * t154;
t12 = t159 * t28 + t163 * t30;
t2 = -qJD(4) * t12 - t159 * t5 + t163 * t8;
t291 = -t12 * t80 - t2;
t289 = t192 * t85;
t288 = t192 * t87;
t60 = -t159 * t192 + t163 * t87;
t266 = qJD(4) * t60;
t14 = -t159 * t35 - t163 * t174 + t266;
t197 = pkin(3) * t164 + pkin(8) * t160;
t253 = t149 + t151;
t102 = t253 * qJD(2);
t287 = -pkin(7) * t174 - t192 * t74;
t214 = t160 * t47 - t164 * t54;
t10 = -qJD(3) * t34 - t214;
t95 = pkin(7) * t158 + t253;
t254 = pkin(2) * t259 + pkin(7) * t262;
t96 = -pkin(1) * t157 - t254;
t272 = t160 * t96 + t164 * t95;
t99 = t157 * t181;
t24 = -qJD(3) * t272 - t101 * t160 + t164 * t99;
t166 = qJD(1) ^ 2;
t211 = t163 * t192;
t58 = t159 * t87 + t211;
t278 = t58 * t80;
t277 = t60 * t58;
t276 = t60 * t80;
t275 = t87 * t85;
t127 = -pkin(2) - t197;
t242 = qJD(4) * t163;
t245 = qJD(3) * t163;
t98 = t198 * t250;
t56 = t160 * t98 + t164 * t97;
t41 = pkin(8) * t225 + t56;
t274 = -t160 * t245 * pkin(7) + t127 * t242 + t294 * t159 - t163 * t41;
t243 = qJD(4) * t159;
t273 = -t127 * t243 + (t246 * pkin(7) + t41) * t159 + t294 * t163;
t13 = qJD(4) * t211 - t159 * t174 + t163 * t35 + t243 * t87;
t271 = t13 * t159;
t270 = t14 * t163;
t32 = qJDD(4) + t36;
t269 = t159 * t32;
t268 = t159 * t80;
t267 = t163 * t32;
t263 = t154 * t166;
t261 = t157 * t162;
t260 = t157 * t164;
t258 = t159 * t164;
t255 = t163 * t165;
t252 = t285 * pkin(1) + pkin(6) * t261;
t155 = t161 ^ 2;
t156 = t165 ^ 2;
t251 = t155 - t156;
t248 = qJD(2) * t161;
t247 = qJD(2) * t165;
t231 = t165 * t263;
t230 = t159 * t259;
t224 = t157 * t248;
t222 = t157 * t158 * t166;
t221 = pkin(1) * t293;
t216 = -pkin(1) * t162 + pkin(6) * t228;
t213 = t163 * t80;
t212 = -t110 * t160 - t164 * t228;
t210 = t165 * t192;
t207 = qJD(2) + 0.2e1 * t249;
t206 = t144 + t238;
t204 = t161 * t231;
t202 = t161 * t219;
t112 = -t158 * t257 + t226;
t71 = t112 * t160 - t162 * t260;
t199 = g(1) * t212 + g(2) * t71;
t195 = -g(1) * t109 + g(2) * t111;
t194 = g(1) * t112 + g(2) * t110;
t191 = -t12 * t159 + t163 * t190;
t107 = -t158 * t164 + t160 * t262;
t108 = t158 * t160 + t161 * t260;
t94 = t145 + (-pkin(2) - t283) * t158;
t37 = pkin(3) * t107 - pkin(8) * t108 + t94;
t39 = -pkin(8) * t259 + t272;
t15 = -t159 * t39 + t163 * t37;
t16 = t159 * t37 + t163 * t39;
t52 = -t160 * t95 + t164 * t96;
t55 = -t160 * t97 + t164 * t98;
t186 = t112 * pkin(2) + pkin(7) * t111 + t252;
t184 = g(1) * t285 + g(2) * t162;
t29 = pkin(3) * t192 - t33;
t183 = -pkin(8) * t32 + t29 * t80;
t65 = t108 * t159 + t157 * t255;
t23 = t164 * t101 + t160 * t99 + t96 * t244 - t246 * t95;
t179 = t160 * t192;
t177 = g(1) * t71 - g(2) * t212 + g(3) * t107;
t72 = t112 * t164 + t160 * t261;
t176 = -g(1) * t72 - g(2) * t68 - g(3) * t108;
t175 = -t110 * pkin(2) - pkin(7) * t109 + t216;
t6 = -pkin(3) * t174 - t10;
t173 = t177 - t6;
t170 = -g(3) * t262 - t194;
t167 = pkin(8) * qJD(4) * t80 - t173;
t105 = t111 * pkin(2);
t103 = t109 * pkin(2);
t100 = t253 * qJD(1);
t93 = t127 * t159 + t163 * t301;
t92 = -pkin(7) * t258 + t127 * t163;
t78 = (t159 * t161 + t164 * t255) * t250;
t77 = t143 * t258 - t163 * t225;
t66 = t108 * t163 - t230;
t64 = -qJD(3) * t107 + t247 * t260;
t63 = qJD(3) * t108 + t157 * t223;
t51 = pkin(3) * t87 + pkin(8) * t85;
t43 = t111 * t159 + t163 * t72;
t42 = t111 * t163 - t159 * t72;
t40 = -pkin(3) * t225 - t55;
t38 = pkin(3) * t259 - t52;
t27 = -qJD(4) * t65 + t159 * t224 + t163 * t64;
t26 = -qJD(4) * t230 + t108 * t242 + t159 * t64 - t163 * t224;
t25 = pkin(3) * t63 - pkin(8) * t64 + t102;
t20 = -pkin(3) * t224 - t24;
t19 = pkin(8) * t224 + t23;
t18 = t159 * t51 + t163 * t33;
t17 = -t159 * t33 + t163 * t51;
t4 = -qJD(4) * t16 - t159 * t19 + t163 * t25;
t3 = qJD(4) * t15 + t159 * t25 + t163 * t19;
t7 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t162 - g(2) * t285, t184, 0, 0, (qJDD(1) * t155 + 0.2e1 * t202) * t154, (t165 * t236 - t239 * t251) * t293, (t161 * t206 + t207 * t247) * t157, (qJDD(1) * t156 - 0.2e1 * t202) * t154, (t165 * t206 - t207 * t248) * t157, t144 * t158, -t102 * t208 + t113 * t144 + t62 * t158 + g(1) * t110 - g(2) * t112 + (-t220 + t237) * t221, -t101 * t208 - t253 * t144 - t61 * t158 + (-t219 - t236) * t221 + t195, ((-t97 * qJD(2) + qJDD(1) * t253 + t61) * t165 + (-qJD(2) * t100 - qJDD(1) * t113 - t62) * t161 - t184) * t157, t154 * qJDD(1) * pkin(1) ^ 2 - g(1) * t216 - g(2) * t252 + t100 * t101 - t97 * t102 + t62 * t113 + t253 * t61, -t108 * t35 + t64 * t87, t107 * t35 - t108 * t36 - t63 * t87 - t64 * t85, -t64 * t192 + t108 * t235 + (t35 * t165 + (qJD(1) * t108 + t87) * t248) * t157, t107 * t36 + t63 * t85, t63 * t192 - t107 * t235 + (t36 * t165 + (-qJD(1) * t107 - t85) * t248) * t157, (-t235 * t165 + (-t143 - t192) * t248) * t157, -t24 * t192 + t52 * t235 + t102 * t85 + t94 * t36 + t48 * t107 + t74 * t63 + g(1) * t68 - g(2) * t72 + (-t10 * t165 + (qJD(1) * t52 + t33) * t248) * t157, t23 * t192 - t272 * t235 + t102 * t87 - t94 * t35 + t48 * t108 + t74 * t64 + (-t182 * t165 + (-qJD(1) * t272 - t34) * t248) * t157 + t199, -t10 * t108 + t107 * t182 - t23 * t85 - t24 * t87 - t272 * t36 - t33 * t64 - t34 * t63 + t35 * t52 - t195, -g(1) * t175 - g(2) * t186 + t10 * t52 + t74 * t102 - t182 * t272 + t34 * t23 + t33 * t24 + t48 * t94, -t13 * t66 + t27 * t60, t13 * t65 - t14 * t66 - t26 * t60 - t27 * t58, -t107 * t13 + t27 * t80 + t32 * t66 + t60 * t63, t14 * t65 + t26 * t58, -t107 * t14 - t26 * t80 - t32 * t65 - t58 * t63, t107 * t32 + t63 * t80, g(1) * t295 - g(2) * t43 + t2 * t107 + t38 * t14 + t15 * t32 - t190 * t63 + t20 * t58 + t29 * t26 + t4 * t80 + t6 * t65, -g(1) * t296 - g(2) * t42 - t1 * t107 - t12 * t63 - t38 * t13 - t16 * t32 + t20 * t60 + t29 * t27 - t3 * t80 + t6 * t66, -t1 * t65 - t12 * t26 + t13 * t15 - t14 * t16 + t190 * t27 - t2 * t66 - t3 * t58 - t4 * t60 - t199, t1 * t16 + t12 * t3 + t2 * t15 - t190 * t4 + t6 * t38 + t29 * t20 - g(1) * (-pkin(3) * t68 + pkin(8) * t212 + t175) - g(2) * (pkin(3) * t72 + pkin(8) * t71 + t186); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t251 * t263, -t165 * t222 + t218, t204, t161 * t222 + t142, t144, t100 * t208 + t263 * t284 - t171 + t62, t97 * t208 + pkin(1) * t231 + (pkin(6) * t239 + g(3)) * t262 + t194 + t229, 0, 0, -t160 * t35 - t164 * t288, (-t35 + t289) * t164 + (-t36 + t288) * t160, -t164 * t209 + t160 * t235 + (t164 * t210 + (-t87 + t240) * t161) * t250, -t164 * t36 - t179 * t85, t160 * t209 + t164 * t235 + (-t160 * t210 + (qJD(2) * t164 + t85) * t161) * t250, t192 * t225, -pkin(2) * t36 - t100 * t85 + t287 * t160 + t299 * t164 + t55 * t192 - t33 * t225, pkin(2) * t35 - t100 * t87 - t299 * t160 + t287 * t164 - t56 * t192 + t34 * t225, t55 * t87 + t56 * t85 + ((qJD(3) * t87 - t36) * pkin(7) + t297) * t164 + (-t10 + t192 * t34 + (qJD(3) * t85 - t35) * pkin(7)) * t160 + t170, -t48 * pkin(2) - t34 * t56 - t33 * t55 - t74 * t100 + g(1) * t105 + g(2) * t103 - g(3) * t254 + (-t10 * t160 - t182 * t164 + (-t160 * t34 - t164 * t33) * qJD(3) - t194) * pkin(7), -t13 * t160 * t163 + (-t160 * t243 + t163 * t244 - t78) * t60, t58 * t78 + t60 * t77 + (-t159 * t60 - t163 * t58) * t244 + (t271 - t270 + (t159 * t58 - t163 * t60) * qJD(4)) * t160, -t78 * t80 + (t245 * t80 + t13) * t164 + (-t192 * t60 - t243 * t80 + t267) * t160, t14 * t159 * t160 + (t159 * t244 + t160 * t242 - t77) * t58, t77 * t80 + (-qJD(3) * t268 + t14) * t164 + (t192 * t58 - t242 * t80 - t269) * t160, -t164 * t32 - t179 * t80, -t29 * t77 + t92 * t32 - t40 * t58 + t273 * t80 + t170 * t159 + (-t2 + (pkin(7) * t58 + t159 * t29) * qJD(3) - t171 * t163) * t164 + (pkin(7) * t14 + t6 * t159 + t190 * t192 + t242 * t29) * t160, -t29 * t78 - t93 * t32 - t40 * t60 - t274 * t80 + t170 * t163 + (t1 + (pkin(7) * t60 + t163 * t29) * qJD(3) + t171 * t159) * t164 + (-pkin(7) * t13 + t12 * t192 + t6 * t163 - t243 * t29) * t160, -t190 * t78 + t12 * t77 + t13 * t92 - t14 * t93 - t273 * t60 - t274 * t58 + t191 * t244 + (-t1 * t159 - t163 * t2 + (-t12 * t163 - t159 * t190) * qJD(4) - t171) * t160, t1 * t93 + t2 * t92 - t29 * t40 - g(1) * (-t197 * t111 - t105) - g(2) * (-t197 * t109 - t103) - g(3) * (t197 * t259 + t254) + t274 * t12 - t273 * t190 + (t160 * t6 + t244 * t29 - t194) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, -t85 ^ 2 + t87 ^ 2, -t35 - t289, -t275, -t288 - t36, t174, -t34 * t143 - t74 * t87 + t177 - t214, t74 * t85 - t176 - t297, 0, 0, t213 * t60 - t271, (-t13 - t278) * t163 + (-t14 - t276) * t159, t213 * t80 - t60 * t87 + t269, t268 * t58 - t270, -t159 * t80 ^ 2 + t58 * t87 + t267, -t80 * t87, -pkin(3) * t14 + t159 * t183 - t163 * t167 - t17 * t80 + t190 * t87 - t34 * t58, pkin(3) * t13 + t12 * t87 + t159 * t167 + t163 * t183 + t18 * t80 - t34 * t60, t17 * t60 + t18 * t58 + ((-t14 + t266) * pkin(8) + t300) * t163 + ((qJD(4) * t58 - t13) * pkin(8) + t291) * t159 + t176, t190 * t17 - t12 * t18 - t29 * t34 + t173 * pkin(3) + (qJD(4) * t191 + t1 * t163 - t2 * t159 + t176) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, -t58 ^ 2 + t60 ^ 2, -t13 + t278, -t277, t276 - t14, t32, -g(1) * t42 + g(2) * t296 + g(3) * t65 - t29 * t60 - t291, g(1) * t43 + g(2) * t295 + g(3) * t66 + t29 * t58 - t300, 0, 0;];
tau_reg = t7;
