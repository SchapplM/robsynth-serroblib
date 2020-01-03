% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:41
% EndTime: 2019-12-31 19:44:52
% DurationCPUTime: 4.26s
% Computational Cost: add. (8211->348), mult. (18338->445), div. (0->0), fcn. (11938->8), ass. (0->219)
t182 = sin(qJ(2));
t185 = cos(qJ(2));
t170 = t182 * qJDD(1);
t226 = qJD(1) * qJD(2);
t218 = t185 * t226;
t157 = t170 + t218;
t179 = sin(pkin(8));
t180 = cos(pkin(8));
t130 = t179 * qJDD(2) + t180 * t157;
t231 = qJD(1) * t182;
t150 = -t180 * qJD(2) + t179 * t231;
t228 = t185 * qJD(1);
t220 = t150 * t228;
t192 = t130 + t220;
t166 = t182 * t226;
t225 = t185 * qJDD(1);
t158 = -t166 + t225;
t152 = t179 * qJD(2) + t180 * t231;
t249 = t152 * t150;
t197 = t158 - t249;
t243 = t179 * t197;
t148 = t152 ^ 2;
t177 = t185 ^ 2;
t187 = qJD(1) ^ 2;
t172 = t177 * t187;
t272 = -t148 - t172;
t72 = -t180 * t272 - t243;
t238 = t180 * t197;
t81 = -t179 * t272 + t238;
t309 = pkin(6) * (t182 * t192 + t185 * t81) + pkin(1) * t72;
t307 = pkin(2) * t72;
t306 = qJ(3) * t72;
t305 = qJ(3) * t81;
t129 = -t180 * qJDD(2) + t179 * t157;
t138 = t152 * t228;
t100 = t129 + t138;
t193 = -t130 + t220;
t294 = -t180 * t100 - t179 * t193;
t265 = t150 ^ 2;
t95 = -t265 - t148;
t303 = -pkin(2) * t95 + qJ(3) * t294;
t134 = t265 - t172;
t302 = t185 * t100 + t182 * (t180 * t134 + t243);
t301 = pkin(6) * (t182 * t95 + t185 * t294);
t107 = t158 + t249;
t237 = t180 * t107;
t268 = -t265 - t172;
t279 = t179 * t268 - t237;
t299 = pkin(2) * t279;
t297 = qJ(3) * t279;
t135 = -t148 + t172;
t242 = t179 * t107;
t296 = t180 * t135 - t242;
t101 = t129 - t138;
t278 = t180 * t268 + t242;
t295 = -pkin(2) * t101 + qJ(3) * t278;
t63 = -t179 * t100 + t180 * t193;
t293 = t179 * t134 - t238;
t240 = t180 * t101;
t244 = t179 * t192;
t291 = t182 * (t240 + t244) + t185 * (t148 - t265);
t290 = t182 * (-t179 * t135 - t237) + t185 * t193;
t289 = pkin(6) * (t182 * t101 + t185 * t278) - pkin(1) * t279;
t288 = pkin(3) * t101;
t181 = sin(qJ(5));
t153 = qJDD(5) + t158;
t184 = cos(qJ(5));
t114 = -t184 * t150 + t181 * t152;
t116 = t181 * t150 + t184 * t152;
t78 = t116 * t114;
t275 = t153 - t78;
t285 = t181 * t275;
t284 = t184 * t275;
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t217 = t183 * g(1) - t186 * g(2);
t143 = qJDD(1) * pkin(1) + t187 * pkin(6) + t217;
t207 = t157 + t218;
t85 = -t207 * qJ(3) + (-t158 + t166) * pkin(2) - t143;
t209 = t186 * g(1) + t183 * g(2);
t250 = qJDD(1) * pkin(6);
t144 = -t187 * pkin(1) - t209 + t250;
t208 = -t185 * pkin(2) - t182 * qJ(3);
t155 = t208 * qJD(1);
t212 = qJD(1) * t155 + t144;
t260 = t182 * g(3);
t264 = qJD(2) ^ 2;
t90 = -t264 * pkin(2) + qJDD(2) * qJ(3) + t212 * t185 - t260;
t216 = t179 * t90 - t180 * t85;
t196 = t158 * pkin(3) - qJ(4) * t172 + qJDD(4) + t216;
t188 = t193 * pkin(7) + t196;
t119 = t150 * pkin(3) - t152 * qJ(4);
t263 = 2 * qJD(3);
t227 = t263 + t119;
t281 = t188 + (pkin(4) * t150 + t227) * t152 + t158 * pkin(4);
t68 = -t114 * qJD(5) + t181 * t129 + t184 * t130;
t164 = qJD(5) + t228;
t98 = t164 * t114;
t276 = t68 - t98;
t273 = t212 * t182;
t258 = t179 * t85 + t180 * t90;
t269 = -t158 * qJ(4) - 0.2e1 * qJD(4) * t228 - t150 * t119 + t258;
t246 = t179 * t101;
t267 = t180 * t192 - t246;
t259 = t185 * g(3);
t203 = -qJDD(2) * pkin(2) - t264 * qJ(3) + qJDD(3) + t259;
t194 = t130 * qJ(4) - t203 - t288;
t235 = t182 * t144;
t266 = -(qJ(4) * t150 * t185 - t155 * t182) * qJD(1) - t194 + t235;
t112 = t114 ^ 2;
t113 = t116 ^ 2;
t162 = t164 ^ 2;
t262 = pkin(3) + pkin(4);
t89 = t203 + t273;
t257 = t179 * t89;
t256 = t180 * t89;
t202 = pkin(4) * t228 - t152 * pkin(7);
t229 = qJD(4) * t152;
t40 = -0.2e1 * t229 + t266;
t33 = t129 * pkin(4) + pkin(7) * t265 - t152 * t202 + t40;
t255 = t181 * t33;
t70 = t153 + t78;
t254 = t181 * t70;
t230 = qJD(3) * t150;
t141 = -0.2e1 * t230;
t205 = t141 + t269;
t34 = -pkin(3) * t172 + t205;
t30 = -pkin(4) * t265 + t129 * pkin(7) - t202 * t228 + t34;
t253 = t184 * t30;
t252 = t184 * t33;
t251 = t184 * t70;
t248 = t164 * t181;
t247 = t164 * t184;
t163 = t185 * t187 * t182;
t234 = t182 * (qJDD(2) + t163);
t232 = t185 * (-t163 + qJDD(2));
t224 = t152 * t263;
t59 = t141 + t258;
t223 = t185 * t78;
t222 = t185 * t249;
t221 = pkin(3) * t180 + pkin(2);
t58 = t216 + t224;
t29 = t179 * t58 + t180 * t59;
t11 = t181 * t30 - t184 * t281;
t124 = t235 + t259;
t125 = t185 * t144 - t260;
t214 = t182 * t124 + t185 * t125;
t213 = -t184 * t129 + t181 * t130;
t211 = t179 * t220;
t210 = t182 * (t180 * t130 + t179 * t138) - t222;
t12 = t281 * t181 + t253;
t7 = -t184 * t11 + t181 * t12;
t8 = t181 * t11 + t184 * t12;
t206 = t179 * t59 - t180 * t58;
t204 = -pkin(1) + t208;
t201 = -t180 * t129 - t211;
t131 = t180 * t138;
t200 = t131 + t211;
t195 = (-qJD(5) + t164) * t116 - t213;
t145 = t185 * t158;
t191 = t145 + t182 * (t150 * t180 - t152 * t179) * t228;
t35 = t227 * t152 + t196;
t189 = t182 * (t179 * t129 - t180 * t220) + t222;
t176 = t182 ^ 2;
t171 = t176 * t187;
t159 = -0.2e1 * t166 + t225;
t156 = t170 + 0.2e1 * t218;
t140 = 0.2e1 * t229;
t94 = -t113 + t162;
t93 = t112 - t162;
t91 = t179 * t130 - t131;
t88 = -t113 - t162;
t77 = t113 - t112;
t76 = -t162 - t112;
t67 = -t116 * qJD(5) - t213;
t62 = (-t114 * t184 + t116 * t181) * t164;
t61 = (t114 * t181 + t116 * t184) * t164;
t60 = -t112 - t113;
t57 = t68 + t98;
t52 = (qJD(5) + t164) * t116 + t213;
t50 = t184 * t93 - t254;
t49 = -t181 * t94 + t284;
t48 = -t181 * t93 - t251;
t47 = -t184 * t94 - t285;
t46 = -t116 * t248 + t184 * t68;
t45 = -t116 * t247 - t181 * t68;
t44 = t114 * t247 - t181 * t67;
t43 = -t114 * t248 - t184 * t67;
t42 = -t181 * t88 - t251;
t41 = t184 * t88 - t254;
t39 = t184 * t76 - t285;
t38 = t181 * t76 + t284;
t37 = t140 - t266 - t288;
t36 = t140 - t273 + (t192 + t220) * qJ(4) + t194;
t32 = -qJ(4) * t95 + t35;
t31 = (-t95 - t172) * pkin(3) + t205;
t27 = t181 * t57 + t184 * t195;
t26 = -t181 * t276 - t184 * t52;
t25 = t181 * t195 - t184 * t57;
t24 = t181 * t52 - t184 * t276;
t22 = t179 * t41 + t180 * t42;
t21 = t179 * t42 - t180 * t41;
t20 = t179 * t38 + t180 * t39;
t19 = t179 * t39 - t180 * t38;
t18 = t179 * t35 + t180 * t34;
t17 = t179 * t34 - t180 * t35;
t16 = -pkin(7) * t41 + qJ(4) * t276 - t252;
t15 = -pkin(7) * t38 + qJ(4) * t52 - t255;
t14 = t179 * t25 + t180 * t27;
t13 = t179 * t27 - t180 * t25;
t10 = -pkin(7) * t42 + t262 * t276 + t255;
t9 = -pkin(7) * t39 + t262 * t52 - t252;
t6 = -pkin(7) * t7 - qJ(4) * t33;
t5 = -pkin(7) * t25 + qJ(4) * t60 - t7;
t4 = -pkin(7) * t8 - t262 * t33;
t3 = -pkin(7) * t27 + t262 * t60 - t8;
t2 = t179 * t7 + t180 * t8;
t1 = t179 * t8 - t180 * t7;
t23 = [0, 0, 0, 0, 0, qJDD(1), t217, t209, 0, 0, t207 * t182, t185 * t156 + t182 * t159, t234 + t185 * (-t171 + t264), -t182 * t218 + t145, t182 * (t172 - t264) + t232, 0, t185 * t143 + pkin(1) * t159 + pkin(6) * (t185 * (-t172 - t264) - t234), -t182 * t143 - pkin(1) * t156 + pkin(6) * (-t232 - t182 * (-t171 - t264)), pkin(1) * (t171 + t172) + (t176 + t177) * t250 + t214, pkin(1) * t143 + pkin(6) * t214, t210, -t291, t290, t189, t302, t191, t182 * (t257 - t297) + t185 * (t58 - t299) + t289, t182 * (t256 + t306) + t185 * (t59 + t307) + t309, -t182 * t206 + t204 * t63 + t301, pkin(6) * (t182 * t89 + t185 * t29) + t204 * t206, t210, t290, t291, t191, -t302, t189, t182 * (-qJ(4) * t240 - t179 * t37 - t297) + t185 * (pkin(3) * t107 - qJ(4) * t268 - t299 + t35) + t289, t182 * (-qJ(3) * t63 - t179 * t31 + t180 * t32) + t185 * (-pkin(2) * t63 - pkin(3) * t193 + qJ(4) * t100) - pkin(1) * t63 + t301, t182 * (-pkin(3) * t244 + t180 * t36 - t306) + t185 * (-t307 + qJ(4) * t197 + 0.2e1 * t230 + (t272 + t172) * pkin(3) - t269) - t309, t182 * (-qJ(3) * t17 + (pkin(3) * t179 - qJ(4) * t180) * t40) + t185 * (-pkin(2) * t17 + pkin(3) * t35 - qJ(4) * t34) - pkin(1) * t17 + pkin(6) * (t185 * t18 + t182 * t40), t182 * (-t179 * t45 + t180 * t46) + t223, t182 * (-t179 * t24 + t180 * t26) + t185 * t77, t182 * (-t179 * t47 + t180 * t49) + t185 * t57, t182 * (-t179 * t43 + t180 * t44) - t223, t182 * (-t179 * t48 + t180 * t50) + t185 * t195, t182 * (-t179 * t61 + t180 * t62) + t185 * t153, t182 * (-qJ(3) * t19 + t180 * t15 - t179 * t9) + t185 * (-pkin(2) * t19 - qJ(4) * t39 + t262 * t38 - t11) - pkin(1) * t19 + pkin(6) * (-t182 * t52 + t185 * t20), t182 * (-qJ(3) * t21 - t179 * t10 + t180 * t16) + t185 * (-pkin(2) * t21 + pkin(3) * t41 - qJ(4) * t42 - t253 - t181 * (t152 * t119 + t188 + t224) + (-t107 * t181 + t41) * pkin(4)) - pkin(1) * t21 + pkin(6) * (-t182 * t276 + t185 * t22), t182 * (-qJ(3) * t13 - t179 * t3 + t180 * t5) + t185 * (-pkin(2) * t13 - qJ(4) * t27 + t262 * t25) - pkin(1) * t13 + pkin(6) * (t185 * t14 - t182 * t60), t182 * (-qJ(3) * t1 - t179 * t4 + t180 * t6) + t185 * (-pkin(2) * t1 - qJ(4) * t8 + t262 * t7) - pkin(1) * t1 + pkin(6) * (t182 * t33 + t185 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t171 - t172, t170, t163, t225, qJDD(2), -t124, -t125, 0, 0, t91, t267, t296, t201, t293, t200, -t256 + t295, -pkin(2) * t192 + t257 + t305, t29 + t303, -pkin(2) * t89 + qJ(3) * t29, t91, t296, -t267, t200, -t293, t201, -qJ(4) * t246 + t180 * t37 + t295, t179 * t32 + t180 * t31 + t303, t179 * t36 + t192 * t221 - t305, qJ(3) * t18 + (-qJ(4) * t179 - t221) * t40, t179 * t46 + t180 * t45, t179 * t26 + t180 * t24, t179 * t49 + t180 * t47, t179 * t44 + t180 * t43, t179 * t50 + t180 * t48, t179 * t62 + t180 * t61, pkin(2) * t52 + qJ(3) * t20 + t179 * t15 + t180 * t9, pkin(2) * t276 + qJ(3) * t22 + t180 * t10 + t179 * t16, pkin(2) * t60 + qJ(3) * t14 + t179 * t5 + t180 * t3, -pkin(2) * t33 + qJ(3) * t2 + t179 * t6 + t180 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t192, t95, t89, 0, 0, 0, 0, 0, 0, t101, t95, -t192, t40, 0, 0, 0, 0, 0, 0, -t52, -t276, -t60, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t193, t272, t35, 0, 0, 0, 0, 0, 0, t38, t41, t25, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t77, t57, -t78, t195, t153, -t11, -t12, 0, 0;];
tauJ_reg = t23;
