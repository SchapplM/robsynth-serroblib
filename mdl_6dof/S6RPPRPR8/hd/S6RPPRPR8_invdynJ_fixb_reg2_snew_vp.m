% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:38:30
% EndTime: 2019-05-05 14:38:43
% DurationCPUTime: 5.48s
% Computational Cost: add. (9089->337), mult. (21250->435), div. (0->0), fcn. (14640->8), ass. (0->210)
t247 = -qJ(3) - pkin(1);
t178 = sin(pkin(9));
t179 = cos(pkin(9));
t187 = qJD(4) ^ 2;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t227 = t179 * t185;
t160 = (-t178 * t182 + t227) * qJD(1);
t254 = t160 ^ 2;
t141 = t254 + t187;
t228 = t179 * t182;
t204 = t178 * t185 + t228;
t158 = t204 * qJD(1);
t229 = t160 * t158;
t273 = qJDD(4) + t229;
t285 = t273 * t182;
t84 = t141 * t185 + t285;
t284 = t273 * t185;
t86 = -t141 * t182 + t284;
t54 = t178 * t86 + t179 * t84;
t300 = t247 * t54;
t299 = pkin(7) * t84;
t298 = pkin(7) * t86;
t142 = t254 - t187;
t274 = qJDD(4) - t229;
t286 = t185 * t274;
t287 = t182 * t274;
t297 = t178 * (t142 * t185 - t287) + t179 * (t142 * t182 + t286);
t255 = t158 ^ 2;
t137 = t255 - t187;
t294 = t178 * (t137 * t182 + t284) + t179 * (-t137 * t185 + t285);
t218 = t179 * qJDD(1);
t219 = t178 * qJDD(1);
t157 = -t182 * t219 + t185 * t218;
t224 = t158 * qJD(4);
t123 = t157 - t224;
t93 = t123 - t224;
t293 = qJ(5) * t93;
t113 = -t187 - t255;
t72 = t113 * t182 + t286;
t75 = -t113 * t185 + t287;
t45 = t178 * t75 - t179 * t72;
t292 = t247 * t45;
t223 = t160 * qJD(4);
t258 = t204 * qJDD(1);
t120 = t258 + 0.2e1 * t223;
t289 = t120 * t182;
t288 = t120 * t185;
t283 = pkin(7) * t72;
t282 = pkin(7) * t75;
t275 = t123 + t224;
t253 = 2 * qJD(5);
t260 = -t254 - t255;
t272 = pkin(3) * t260;
t108 = pkin(4) * t158 - qJ(5) * t160;
t188 = qJD(1) ^ 2;
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t213 = t183 * g(1) - t186 * g(2);
t205 = qJDD(2) - t213;
t200 = -t188 * qJ(2) + t205;
t209 = -0.2e1 * qJD(1) * qJD(3) + t247 * qJDD(1) + t200;
t248 = t178 * g(3);
t251 = pkin(3) * t188;
t95 = t248 + (-pkin(7) * qJDD(1) - t178 * t251 + t209) * t179;
t112 = -g(3) * t179 + t209 * t178;
t173 = t178 ^ 2;
t99 = -pkin(7) * t219 - t173 * t251 + t112;
t65 = t182 * t99 - t185 * t95;
t43 = -qJDD(4) * pkin(4) - t187 * qJ(5) + t160 * t108 + qJDD(5) + t65;
t271 = -pkin(8) * t274 + t43;
t270 = qJ(2) * t260;
t181 = sin(qJ(6));
t109 = qJDD(6) + t123;
t184 = cos(qJ(6));
t130 = qJD(4) * t181 - t184 * t158;
t132 = qJD(4) * t184 + t158 * t181;
t98 = t132 * t130;
t263 = t109 - t98;
t269 = t181 * t263;
t268 = t182 * t258;
t267 = t184 * t263;
t266 = t185 * t258;
t174 = t179 ^ 2;
t226 = t173 + t174;
t265 = -pkin(3) * t219 + t188 * (t226 * pkin(7) - t247);
t262 = t226 * t188;
t261 = -pkin(4) * t223 + t160 * t253;
t259 = t254 - t255;
t128 = t130 ^ 2;
t129 = t132 ^ 2;
t150 = qJD(6) + t160;
t148 = t150 ^ 2;
t252 = pkin(4) + pkin(8);
t250 = pkin(4) * t185;
t121 = t258 + t223;
t249 = t121 * pkin(4);
t66 = t182 * t95 + t185 * t99;
t38 = t182 * t66 - t185 * t65;
t246 = t179 * t38;
t135 = pkin(5) * t160 - qJD(4) * pkin(8);
t198 = -t187 * pkin(4) - t158 * t108 + t66;
t220 = qJDD(4) * qJ(5);
t31 = t220 - t121 * pkin(5) - t255 * pkin(8) + (t253 + t135) * qJD(4) + t198;
t245 = t181 * t31;
t69 = t109 + t98;
t244 = t181 * t69;
t221 = qJD(2) * qJD(1);
t171 = 0.2e1 * t221;
t175 = qJDD(1) * qJ(2);
t206 = t186 * g(1) + t183 * g(2);
t202 = -t175 + t206;
t201 = -qJDD(3) + t202;
t191 = t171 - t201 - t261 - t265 - t293;
t29 = -pkin(5) * t255 + t252 * t121 - t160 * t135 + t191;
t243 = t184 * t29;
t242 = t184 * t31;
t241 = t184 * t69;
t240 = qJDD(1) * pkin(1);
t197 = t201 - 0.2e1 * t221;
t103 = t197 + t265;
t239 = t103 * t182;
t238 = t103 * t185;
t231 = t150 * t181;
t230 = t150 * t184;
t222 = qJD(6) + t150;
t216 = t182 * t98;
t215 = t185 * t98;
t214 = -qJ(5) * t182 - pkin(3);
t192 = pkin(5) * t275 + t271;
t17 = t181 * t29 - t184 * t192;
t39 = t182 * t65 + t185 * t66;
t134 = -t247 * t188 + t197;
t210 = -t134 + t175;
t208 = t181 * qJDD(4) - t184 * t121;
t18 = t181 * t192 + t243;
t7 = -t184 * t17 + t181 * t18;
t8 = t181 * t17 + t184 * t18;
t196 = qJD(4) * t253 + t198;
t42 = t196 + t220;
t11 = t178 * (t182 * t43 + t185 * t42) + t179 * (t182 * t42 - t185 * t43);
t71 = t179 * (t209 * t179 + t248) + t178 * t112;
t203 = t184 * qJDD(4) + t181 * t121;
t199 = (-qJD(6) + t150) * t132 - t208;
t82 = -qJD(6) * t130 + t203;
t195 = t179 * (t185 * t123 - t182 * t223) - t178 * (t182 * t123 + t185 * t223);
t194 = t179 * (t121 * t182 + t185 * t224) - t178 * (-t185 * t121 + t182 * t224);
t193 = (t179 * (-t158 * t185 + t160 * t182) - t178 * (-t158 * t182 - t160 * t185)) * qJD(4);
t164 = t226 * qJDD(1);
t163 = t178 * t262;
t162 = t179 * t262;
t151 = -t200 + t240;
t122 = t157 - 0.2e1 * t224;
t102 = t150 * t130;
t101 = -t129 + t148;
t100 = t128 - t148;
t94 = t129 - t128;
t88 = -t129 - t148;
t83 = -t148 - t128;
t81 = -qJD(6) * t132 - t208;
t80 = -t128 - t129;
t79 = t182 * t275 - t266;
t78 = t157 * t182 - t266;
t77 = -t185 * t275 - t268;
t76 = -t157 * t185 - t268;
t67 = (t130 * t181 + t132 * t184) * t150;
t63 = -t222 * t130 + t203;
t62 = t102 + t82;
t61 = -t102 + t82;
t58 = t222 * t132 + t208;
t57 = -t132 * t230 - t181 * t82;
t56 = -t130 * t231 - t184 * t81;
t53 = -t100 * t181 - t241;
t52 = -t101 * t184 - t269;
t51 = -t181 * t88 - t241;
t50 = t184 * t88 - t244;
t49 = t178 * t79 + t179 * t77;
t48 = t178 * t78 + t179 * t76;
t47 = t184 * t83 - t269;
t46 = t181 * t83 + t267;
t41 = t191 + t249;
t40 = -qJ(5) * t260 + t43;
t37 = -pkin(4) * t260 + t42;
t36 = (t121 + t120) * pkin(4) + t191;
t35 = t103 - t249 + t261 + 0.2e1 * t293;
t34 = t181 * t62 + t184 * t199;
t33 = t181 * t199 - t184 * t62;
t32 = t181 * t58 - t184 * t61;
t28 = t182 * t50 + t185 * t63;
t27 = t182 * t63 - t185 * t50;
t26 = t182 * t46 + t185 * t58;
t25 = t182 * t58 - t185 * t46;
t24 = t182 * t33 + t185 * t80;
t23 = t182 * t80 - t185 * t33;
t20 = t178 * t39 + t246;
t19 = pkin(5) * t33 - qJ(5) * t34;
t16 = t178 * t28 + t179 * t27;
t15 = t178 * t26 + t179 * t25;
t14 = t178 * t24 + t179 * t23;
t13 = pkin(5) * t63 - t252 * t51 - t245;
t12 = pkin(5) * t58 - t252 * t47 + t242;
t10 = -t243 - t181 * t271 - qJ(5) * t51 + (-t181 * t275 + t50) * pkin(5);
t9 = pkin(5) * t46 - qJ(5) * t47 - t17;
t6 = t182 * t7 + t185 * t31;
t5 = t182 * t31 - t185 * t7;
t4 = pkin(5) * t80 - t252 * t34 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t31 - t252 * t8;
t1 = t178 * t6 + t179 * t5;
t21 = [0, 0, 0, 0, 0, qJDD(1), t213, t206, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t205 - 0.2e1 * t240, t171 + 0.2e1 * t175 - t206, pkin(1) * t151 + qJ(2) * (-t188 * pkin(1) + t171 - t202), t174 * qJDD(1), -0.2e1 * t178 * t218, 0, t173 * qJDD(1), 0, 0, -t247 * t163 + t210 * t178, -t247 * t162 + t210 * t179, -qJ(2) * t262 - t247 * t164 - t71, -qJ(2) * t134 + t247 * t71, t195, t179 * (-t122 * t182 - t288) - t178 * (t122 * t185 - t289), t297, t194, -t294, t193, t179 * (-t239 - t283) - t178 * (-pkin(3) * t120 + t238 - t282) + qJ(2) * t120 - t292, t179 * (-t238 + t299) - t178 * (-pkin(3) * t122 - t239 - t298) + qJ(2) * t122 - t300, t179 * (-pkin(7) * t76 - t38) - t178 * (pkin(7) * t78 - t272 + t39) + t270 + t247 * t48, -pkin(7) * t246 - t178 * (pkin(3) * t103 + pkin(7) * t39) - qJ(2) * t103 + t247 * t20, t193, -t297, t294, t195, t179 * (-t182 * t93 - t288) - t178 * (t185 * t93 - t289), t194, t179 * (-pkin(7) * t77 - t182 * t37 + t185 * t40) - t178 * (pkin(7) * t79 + t182 * t40 + t185 * t37 - t272) + t270 + t247 * t49, t179 * (-t182 * t36 + t283) - t178 * (t185 * t36 + t282) - (-qJ(5) * t227 - t178 * t214 + qJ(2)) * t120 + t292, t179 * (t185 * t35 - t299) - t178 * (t182 * t35 + t298) + (-pkin(4) * t228 - t178 * (pkin(3) + t250) - qJ(2)) * t93 + t300, (t179 * (pkin(4) * t182 - qJ(5) * t185) - t178 * (t214 - t250) + qJ(2)) * t41 + (t247 - pkin(7)) * t11, t179 * (-t182 * t57 + t215) - t178 * (t185 * t57 + t216), t179 * (-t182 * t32 + t185 * t94) - t178 * (t182 * t94 + t185 * t32), t179 * (-t182 * t52 + t185 * t62) - t178 * (t182 * t62 + t185 * t52), t179 * (-t182 * t56 - t215) - t178 * (t185 * t56 - t216), t179 * (-t182 * t53 + t185 * t199) - t178 * (t182 * t199 + t185 * t53), t179 * (t109 * t185 - t182 * t67) - t178 * (t109 * t182 + t185 * t67), t179 * (-pkin(7) * t25 - t12 * t182 + t185 * t9) - t178 * (-pkin(3) * t47 + pkin(7) * t26 + t12 * t185 + t182 * t9) + qJ(2) * t47 + t247 * t15, t179 * (-pkin(7) * t27 + t10 * t185 - t13 * t182) - t178 * (-pkin(3) * t51 + pkin(7) * t28 + t10 * t182 + t13 * t185) + qJ(2) * t51 + t247 * t16, t179 * (-pkin(7) * t23 - t182 * t4 + t185 * t19) - t178 * (-pkin(3) * t34 + pkin(7) * t24 + t182 * t19 + t185 * t4) + qJ(2) * t34 + t247 * t14, t179 * (-pkin(7) * t5 - t182 * t2 + t185 * t3) - t178 * (-pkin(3) * t8 + pkin(7) * t6 + t182 * t3 + t185 * t2) + qJ(2) * t8 + t247 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t188, -t151, 0, 0, 0, 0, 0, 0, -t163, -t162, -t164, t71, 0, 0, 0, 0, 0, 0, -t45, -t54, t48, t20, 0, 0, 0, 0, 0, 0, t49, t45, t54, t11, 0, 0, 0, 0, 0, 0, t15, t16, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t218, -t262, -t134, 0, 0, 0, 0, 0, 0, t120, t122, t260, -t103, 0, 0, 0, 0, 0, 0, t260, -t120, -t93, t41, 0, 0, 0, 0, 0, 0, t47, t51, t34, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t259, t157, -t229, -t258, qJDD(4), -t65, -t66, 0, 0, qJDD(4), -t275, t258, t229, t259, -t229, -pkin(4) * t275 - qJ(5) * t258, -pkin(4) * t274 - qJ(5) * t113 + t43, pkin(4) * t141 + (qJDD(4) + t273) * qJ(5) + t196, -pkin(4) * t43 + qJ(5) * t42, -t132 * t231 + t184 * t82, -t181 * t61 - t184 * t58, -t101 * t181 + t267, t130 * t230 - t181 * t81, t100 * t184 - t244, (-t130 * t184 + t132 * t181) * t150, qJ(5) * t58 - t252 * t46 + t245, qJ(5) * t63 - t252 * t50 + t242, qJ(5) * t80 - t252 * t33 - t7, qJ(5) * t31 - t252 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t274, -t141, t43, 0, 0, 0, 0, 0, 0, t46, t50, t33, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t94, t62, -t98, t199, t109, -t17, -t18, 0, 0;];
tauJ_reg  = t21;
