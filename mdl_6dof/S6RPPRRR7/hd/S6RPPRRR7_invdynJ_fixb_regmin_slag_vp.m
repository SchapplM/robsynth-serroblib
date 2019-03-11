% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:51
% EndTime: 2019-03-09 02:34:00
% DurationCPUTime: 3.30s
% Computational Cost: add. (4723->341), mult. (9937->431), div. (0->0), fcn. (7823->14), ass. (0->197)
t163 = cos(qJ(6));
t225 = qJD(6) * t163;
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t156 = sin(pkin(10));
t157 = cos(pkin(10));
t165 = cos(qJ(4));
t231 = qJD(1) * t165;
t161 = sin(qJ(4));
t232 = qJD(1) * t161;
t96 = -t156 * t231 - t157 * t232;
t214 = t156 * t232;
t216 = t157 * t231;
t97 = -t214 + t216;
t61 = t160 * t97 - t164 * t96;
t291 = t163 * t61;
t302 = t225 + t291;
t152 = qJD(4) + qJD(5);
t255 = t152 * t61;
t227 = qJD(5) * t164;
t228 = qJD(5) * t160;
t222 = t157 * qJDD(1);
t223 = t156 * qJDD(1);
t189 = -t161 * t223 + t165 * t222;
t105 = t156 * t165 + t157 * t161;
t98 = t105 * qJD(4);
t68 = -qJD(1) * t98 + t189;
t176 = -qJD(4) * t214 + t105 * qJDD(1);
t229 = qJD(4) * t165;
t215 = t157 * t229;
t69 = qJD(1) * t215 + t176;
t31 = -t160 * t69 + t164 * t68 + t96 * t227 - t97 * t228;
t301 = t31 + t255;
t148 = qJDD(4) + qJDD(5);
t159 = sin(qJ(6));
t188 = t160 * t96 + t164 * t97;
t226 = qJD(6) * t159;
t13 = t159 * t148 + t152 * t225 + t163 * t31 - t188 * t226;
t106 = -t156 * t161 + t157 * t165;
t277 = -t105 * t160 + t164 * t106;
t185 = t164 * t105 + t106 * t160;
t230 = qJD(4) * t161;
t99 = -t156 * t230 + t215;
t42 = qJD(5) * t185 + t160 * t99 + t164 * t98;
t55 = t152 * t159 + t163 * t188;
t300 = t13 * t277 - t42 * t55;
t299 = t148 * t277 - t152 * t42;
t14 = qJD(6) * t55 - t163 * t148 + t159 * t31;
t53 = -t163 * t152 + t159 * t188;
t298 = t13 * t163 - t159 * t14 - t302 * t53;
t11 = t13 * t159;
t297 = t302 * t55 + t11;
t260 = t55 * t188;
t32 = qJD(5) * t188 + t160 * t68 + t164 * t69;
t30 = qJDD(6) + t32;
t27 = t159 * t30;
t289 = qJD(6) + t61;
t58 = t289 * t225;
t296 = t289 * t291 - t260 + t27 + t58;
t158 = -pkin(1) - qJ(3);
t118 = t158 * qJD(1) + qJD(2);
t204 = -pkin(7) * qJD(1) + t118;
t90 = t204 * t156;
t91 = t204 * t157;
t186 = -t161 * t91 - t165 * t90;
t49 = pkin(8) * t96 - t186;
t253 = t160 * t49;
t278 = -t161 * t90 + t165 * t91;
t48 = -pkin(8) * t97 + t278;
t47 = qJD(4) * pkin(4) + t48;
t23 = t164 * t47 - t253;
t19 = -pkin(5) * t152 - t23;
t294 = t19 * t61;
t151 = pkin(10) + qJ(4);
t141 = qJ(5) + t151;
t133 = cos(t141);
t162 = sin(qJ(1));
t147 = g(1) * t162;
t219 = t133 * t147;
t267 = -qJD(1) * qJD(3) + qJDD(1) * t158;
t109 = qJDD(2) + t267;
t199 = -pkin(7) * qJDD(1) + t109;
t76 = t199 * t156;
t77 = t199 * t157;
t205 = -t161 * t76 + t165 * t77;
t21 = qJDD(4) * pkin(4) - pkin(8) * t68 + qJD(4) * t186 + t205;
t187 = t161 * t77 + t165 * t76;
t22 = -pkin(8) * t69 + t278 * qJD(4) + t187;
t249 = t164 * t49;
t24 = t160 * t47 + t249;
t269 = qJD(5) * t24 + t160 * t22 - t164 * t21;
t3 = -pkin(5) * t148 + t269;
t293 = t3 + t219;
t290 = t188 * t61;
t132 = sin(t141);
t166 = cos(qJ(1));
t146 = g(2) * t166;
t288 = g(3) * t132 + t133 * t146;
t274 = t147 - t146;
t256 = t152 * t188;
t287 = -t32 + t256;
t285 = t188 ^ 2 - t61 ^ 2;
t127 = g(3) * t133;
t268 = (qJD(5) * t47 + t22) * t164 + t160 * t21 - t49 * t228;
t138 = qJD(1) * qJ(2) + qJD(3);
t142 = t156 * pkin(3);
t112 = qJD(1) * t142 + t138;
t74 = -pkin(4) * t96 + t112;
t284 = t274 * t132 + t61 * t74 + t127 - t268;
t282 = pkin(5) * t188;
t136 = pkin(4) * t160 + pkin(9);
t281 = (pkin(4) * t97 + pkin(9) * t61 + qJD(6) * t136 + t282) * t289;
t280 = (t289 * pkin(9) + t282) * t289;
t261 = t53 * t188;
t279 = t289 * t188;
t233 = t156 ^ 2 + t157 ^ 2;
t276 = t118 * t233;
t259 = -pkin(7) + t158;
t110 = t259 * t156;
t111 = t259 * t157;
t236 = t165 * t110 + t161 * t111;
t153 = qJDD(1) * qJ(2);
t154 = qJD(1) * qJD(2);
t273 = t153 + t154;
t116 = qJDD(3) + t273;
t192 = g(1) * t166 + g(2) * t162;
t275 = t116 - t192;
t20 = pkin(9) * t152 + t24;
t33 = pkin(5) * t61 - pkin(9) * t188 + t74;
t6 = -t159 * t20 + t163 * t33;
t272 = t288 * t163 - t6 * t188 + t19 * t226;
t7 = t159 * t33 + t163 * t20;
t271 = t293 * t159 + t7 * t188 + t19 * t225;
t270 = -t188 * t74 - t219 - t269 + t288;
t170 = t277 * qJD(5) - t160 * t98 + t164 * t99;
t266 = -t148 * t185 - t152 * t170;
t212 = pkin(9) * t148 + qJD(6) * t33 + t268;
t197 = -t110 * t161 + t165 * t111;
t56 = -pkin(8) * t106 + t197;
t57 = -pkin(8) * t105 + t236;
t36 = t160 * t56 + t164 * t57;
t131 = qJ(2) + t142;
t79 = pkin(4) * t105 + t131;
t37 = pkin(5) * t185 - pkin(9) * t277 + t79;
t35 = t160 * t57 - t164 * t56;
t172 = -qJD(3) * t105 - t110 * t230 + t111 * t229;
t45 = -pkin(8) * t99 + t172;
t171 = -t106 * qJD(3) - t236 * qJD(4);
t46 = pkin(8) * t98 + t171;
t8 = -qJD(5) * t35 + t160 * t46 + t164 * t45;
t265 = -t185 * t212 - t19 * t42 - (qJD(6) * t37 + t8) * t289 + t3 * t277 - t36 * t30;
t264 = 0.2e1 * t154;
t263 = t19 * t277;
t262 = t37 * t30;
t254 = t159 * t55;
t28 = t163 * t30;
t248 = -t98 * qJD(4) + t106 * qJDD(4);
t247 = pkin(1) * qJDD(1);
t242 = t159 * t162;
t241 = t159 * t166;
t240 = t162 * t163;
t239 = t163 * t166;
t235 = t166 * pkin(1) + t162 * qJ(2);
t217 = t277 * t226;
t81 = pkin(4) * t99 + qJD(2);
t101 = pkin(3) * t223 + t116;
t52 = pkin(4) * t69 + t101;
t5 = pkin(5) * t32 - pkin(9) * t31 + t52;
t211 = qJD(6) * t20 - t5;
t201 = t159 * t289;
t200 = t233 * t109;
t196 = qJD(6) * t185 + qJD(1);
t195 = qJDD(2) - t247;
t25 = t160 * t48 + t249;
t193 = pkin(4) * t228 - t25;
t190 = t277 * t30 - t289 * t42;
t184 = t28 - (t159 * t61 + t226) * t289;
t183 = -qJD(4) * t99 - qJDD(4) * t105;
t181 = -t212 + t127;
t179 = -pkin(9) * t30 + t23 * t289 + t294;
t177 = t105 * qJD(1);
t26 = t164 * t48 - t253;
t174 = -t136 * t30 + t294 + (-pkin(4) * t227 + t26) * t289;
t173 = t273 + t275;
t167 = qJD(1) ^ 2;
t144 = t166 * qJ(2);
t140 = cos(t151);
t139 = sin(t151);
t137 = -pkin(4) * t164 - pkin(5);
t89 = t132 * t239 - t242;
t88 = t132 * t241 + t240;
t87 = t132 * t240 + t241;
t86 = -t132 * t242 + t239;
t15 = pkin(5) * t170 + pkin(9) * t42 + t81;
t9 = qJD(5) * t36 + t160 * t45 - t164 * t46;
t4 = t163 * t5;
t1 = [qJDD(1), t274, t192, qJDD(2) - t274 - 0.2e1 * t247, 0.2e1 * t153 + t264 - t192, -t195 * pkin(1) - g(1) * (-pkin(1) * t162 + t144) - g(2) * t235 + (t153 + t264) * qJ(2), t173 * t156, t173 * t157, t274 + t233 * (-t109 - t267) t116 * qJ(2) + t138 * qJD(2) - g(1) * (t158 * t162 + t144) - g(2) * (qJ(3) * t166 + t235) + t158 * t200 - qJD(3) * t276, t106 * t68 - t97 * t98, -t105 * t68 - t106 * t69 - t96 * t98 - t97 * t99, t248, t183, 0, -qJD(2) * t96 + qJD(4) * t171 + qJDD(4) * t197 + t101 * t105 + t112 * t99 + t131 * t69 - t139 * t192, qJD(2) * t97 - qJD(4) * t172 - qJDD(4) * t236 + t101 * t106 - t112 * t98 + t131 * t68 - t140 * t192, -t188 * t42 + t277 * t31, -t170 * t188 - t185 * t31 - t277 * t32 + t42 * t61, t299, t266, 0, -t132 * t192 - t148 * t35 - t152 * t9 + t170 * t74 + t185 * t52 + t32 * t79 + t61 * t81, -t133 * t192 - t148 * t36 - t152 * t8 + t188 * t81 + t277 * t52 + t31 * t79 - t42 * t74, t163 * t300 - t55 * t217 -(-t163 * t53 - t254) * t42 + (-t11 - t14 * t163 + (t159 * t53 - t163 * t55) * qJD(6)) * t277, t13 * t185 + t163 * t190 + t170 * t55 - t217 * t289, -t14 * t185 - t159 * t190 - t170 * t53 - t277 * t58, t170 * t289 + t185 * t30, -g(1) * t89 - g(2) * t87 + t35 * t14 + t4 * t185 + t6 * t170 + t9 * t53 + (t15 * t289 + t262 + (-t185 * t20 - t289 * t36 + t263) * qJD(6)) * t163 + t265 * t159, g(1) * t88 - g(2) * t86 + t35 * t13 - t7 * t170 + t9 * t55 + (-(-qJD(6) * t36 + t15) * t289 - t262 + t211 * t185 - qJD(6) * t263) * t159 + t265 * t163; 0, 0, 0, qJDD(1), -t167, -qJ(2) * t167 + t195 - t274, -t167 * t156, -t167 * t157, -t233 * qJDD(1), -qJD(1) * t138 + t200 - t274, 0, 0, 0, 0, 0, qJD(1) * t96 + t248, -qJD(1) * t97 + t183, 0, 0, 0, 0, 0, -qJD(1) * t61 + t299, -qJD(1) * t188 + t266, 0, 0, 0, 0, 0, -t185 * t27 - t14 * t277 + t42 * t53 + (-t159 * t170 - t163 * t196) * t289, -t185 * t28 + (t159 * t196 - t163 * t170) * t289 - t300; 0, 0, 0, 0, 0, 0, t223, t222, -t233 * t167, qJD(1) * t276 + t275, 0, 0, 0, 0, 0 (t97 + t216) * qJD(4) + t176 (t96 - t177) * qJD(4) + t189, 0, 0, 0, 0, 0, t32 + t256, t31 - t255, 0, 0, 0, 0, 0, t184 - t261, -t163 * t289 ^ 2 - t260 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t96, -t96 ^ 2 + t97 ^ 2 (-t96 - t177) * qJD(4) + t189 (t97 - t216) * qJD(4) - t176, qJDD(4), g(3) * t139 - t112 * t97 - t140 * t274 + t205, g(3) * t140 - t112 * t96 + t139 * t274 - t187, t290, t285, t301, t287, t148, t152 * t25 + (t148 * t164 - t152 * t228 - t61 * t97) * pkin(4) + t270, t152 * t26 + (-t148 * t160 - t152 * t227 - t188 * t97) * pkin(4) + t284, t297, -t254 * t289 + t298, t296, t184 + t261, -t279, t137 * t14 + t193 * t53 + (-t293 - t281) * t163 + t174 * t159 + t272, t137 * t13 + t193 * t55 + t174 * t163 + (-t288 + t281) * t159 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t285, t301, t287, t148, t152 * t24 + t270, t152 * t23 + t284, t297, -t201 * t55 + t298, t296, -t201 * t289 + t261 + t28, -t279, -pkin(5) * t14 - t24 * t53 + t179 * t159 + (-t293 - t280) * t163 + t272, -pkin(5) * t13 - t24 * t55 + t179 * t163 + (-t288 + t280) * t159 + t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t53 ^ 2 + t55 ^ 2, t289 * t53 + t13, t289 * t55 - t14, t30, -g(1) * t86 - g(2) * t88 + t159 * t181 - t19 * t55 - t20 * t225 + t289 * t7 + t4, g(1) * t87 - g(2) * t89 + t159 * t211 + t163 * t181 + t19 * t53 + t289 * t6;];
tau_reg  = t1;
