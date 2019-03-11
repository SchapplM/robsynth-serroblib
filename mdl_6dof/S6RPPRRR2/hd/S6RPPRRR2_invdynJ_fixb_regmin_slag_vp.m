% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:25
% EndTime: 2019-03-09 02:21:34
% DurationCPUTime: 3.56s
% Computational Cost: add. (4428->380), mult. (9759->496), div. (0->0), fcn. (7655->18), ass. (0->211)
t185 = cos(qJ(4));
t177 = cos(pkin(11));
t243 = qJD(1) * t177;
t150 = t185 * t243;
t175 = sin(pkin(11));
t181 = sin(qJ(4));
t247 = t175 * t181;
t228 = qJD(1) * t247;
t125 = t150 - t228;
t284 = qJD(5) + qJD(6);
t297 = t125 - t284;
t176 = sin(pkin(10));
t151 = pkin(1) * t176 + qJ(3);
t143 = t151 * qJD(1);
t160 = t177 * qJD(2);
t105 = t160 + (-pkin(7) * qJD(1) - t143) * t175;
t119 = qJD(2) * t175 + t143 * t177;
t106 = pkin(7) * t243 + t119;
t47 = t105 * t181 + t106 * t185;
t296 = qJD(4) * t47;
t122 = qJD(5) - t125;
t172 = pkin(11) + qJ(4);
t161 = sin(t172);
t163 = cos(t172);
t173 = qJ(1) + pkin(10);
t162 = sin(t173);
t164 = cos(t173);
t217 = g(1) * t164 + g(2) * t162;
t193 = -g(3) * t163 + t161 * t217;
t131 = qJD(1) * qJD(3) + qJDD(1) * t151;
t158 = t177 * qJDD(2);
t93 = t158 + (-pkin(7) * qJDD(1) - t131) * t175;
t111 = qJDD(2) * t175 + t131 * t177;
t234 = t177 * qJDD(1);
t94 = pkin(7) * t234 + t111;
t207 = -t181 * t94 + t185 * t93 - t296;
t20 = -qJDD(4) * pkin(4) - t207;
t295 = -qJD(5) * pkin(8) * t122 + t193 - t20;
t180 = sin(qJ(5));
t183 = cos(qJ(6));
t179 = sin(qJ(6));
t184 = cos(qJ(5));
t246 = t179 * t184;
t136 = t180 * t183 + t246;
t271 = t297 * t136;
t134 = t175 * t185 + t177 * t181;
t126 = t134 * qJD(1);
t238 = t184 * qJD(4);
t107 = t126 * t180 - t238;
t109 = qJD(4) * t180 + t126 * t184;
t205 = t107 * t179 - t109 * t183;
t49 = t107 * t183 + t109 * t179;
t294 = t205 * t49;
t242 = qJD(5) * t180;
t259 = t125 * t180;
t293 = t242 - t259;
t292 = t205 ^ 2 - t49 ^ 2;
t120 = qJD(6) + t122;
t239 = qJD(6) * t183;
t240 = qJD(6) * t179;
t235 = t175 * qJDD(1);
t231 = qJD(4) * t150 + t181 * t234 + t185 * t235;
t84 = -qJD(4) * t228 + t231;
t40 = qJD(5) * t238 + qJDD(4) * t180 - t126 * t242 + t184 * t84;
t41 = qJD(5) * t109 - qJDD(4) * t184 + t180 * t84;
t9 = -t107 * t239 - t109 * t240 - t179 * t41 + t183 * t40;
t291 = t120 * t49 + t9;
t174 = qJ(5) + qJ(6);
t168 = sin(t174);
t251 = t164 * t168;
t169 = cos(t174);
t254 = t162 * t169;
t102 = -t163 * t254 + t251;
t250 = t164 * t169;
t255 = t162 * t168;
t104 = t163 * t250 + t255;
t44 = qJD(4) * pkin(8) + t47;
t178 = cos(pkin(10));
t153 = -pkin(1) * t178 - pkin(2);
t142 = -pkin(3) * t177 + t153;
t124 = qJD(1) * t142 + qJD(3);
t58 = -t125 * pkin(4) - t126 * pkin(8) + t124;
t25 = t180 * t58 + t184 * t44;
t15 = -pkin(9) * t107 + t25;
t13 = t15 * t240;
t279 = g(3) * t161;
t286 = t105 * t185 - t106 * t181;
t43 = -qJD(4) * pkin(4) - t286;
t32 = pkin(5) * t107 + t43;
t290 = g(1) * t104 - g(2) * t102 + t169 * t279 + t32 * t49 + t13;
t101 = t163 * t255 + t250;
t103 = -t163 * t251 + t254;
t209 = t181 * t93 + t185 * t94;
t19 = qJDD(4) * pkin(8) + qJD(4) * t286 + t209;
t121 = qJDD(1) * t142 + qJDD(3);
t128 = t134 * qJD(4);
t211 = t181 * t235 - t185 * t234;
t85 = qJD(1) * t128 + t211;
t31 = t85 * pkin(4) - t84 * pkin(8) + t121;
t30 = t184 * t31;
t81 = qJDD(5) + t85;
t2 = t81 * pkin(5) - t40 * pkin(9) - qJD(5) * t25 - t180 * t19 + t30;
t241 = qJD(5) * t184;
t200 = t180 * t31 + t184 * t19 + t241 * t58 - t242 * t44;
t3 = -pkin(9) * t41 + t200;
t230 = -t179 * t3 + t183 * t2;
t24 = -t180 * t44 + t184 * t58;
t14 = -pkin(9) * t109 + t24;
t12 = pkin(5) * t122 + t14;
t268 = t15 * t183;
t5 = t12 * t179 + t268;
t289 = -g(1) * t103 + g(2) * t101 - qJD(6) * t5 + t168 * t279 + t205 * t32 + t230;
t191 = qJD(6) * t205 - t179 * t40 - t183 * t41;
t288 = -t120 * t205 + t191;
t74 = t136 * t134;
t236 = qJDD(1) * t153;
t139 = qJDD(3) + t236;
t216 = g(1) * t162 - g(2) * t164;
t285 = t216 - t139;
t135 = t179 * t180 - t183 * t184;
t272 = t297 * t135;
t79 = qJDD(6) + t81;
t283 = -t120 * t272 - t136 * t79;
t227 = t134 * t242;
t133 = -t177 * t185 + t247;
t127 = t133 * qJD(4);
t258 = t127 * t184;
t196 = t227 + t258;
t256 = t134 * t184;
t282 = -t122 * t196 + t256 * t81;
t281 = pkin(8) + pkin(9);
t280 = -t128 * t205 + t133 * t9;
t277 = pkin(7) + t151;
t245 = t180 * t127;
t257 = t134 * t180;
t22 = -t127 * t246 - t179 * t227 - t240 * t257 + (t256 * t284 - t245) * t183;
t276 = -t120 * t22 - t74 * t79;
t275 = t109 * t128 + t133 * t40;
t82 = pkin(4) * t126 - pkin(8) * t125;
t274 = t180 * t82 + t184 * t286;
t129 = t277 * t175;
t130 = t277 * t177;
t77 = -t129 * t181 + t130 * t185;
t67 = t184 * t77;
t68 = pkin(4) * t133 - pkin(8) * t134 + t142;
t273 = t180 * t68 + t67;
t270 = t126 * t49;
t267 = t180 * t81;
t266 = t40 * t180;
t265 = t205 * t126;
t263 = t107 * t122;
t262 = t107 * t126;
t261 = t109 * t122;
t260 = t109 * t126;
t253 = t162 * t180;
t252 = t162 * t184;
t249 = t164 * t180;
t248 = t164 * t184;
t244 = t175 ^ 2 + t177 ^ 2;
t229 = qJD(5) * t281;
t226 = t134 * t241;
t225 = qJD(6) * t12 + t3;
t223 = -qJD(5) * t58 - t19;
t221 = t122 * t184;
t220 = t120 * t271 - t135 * t79;
t219 = pkin(5) * t293 - t47;
t218 = -t241 * t44 + t30;
t182 = sin(qJ(1));
t186 = cos(qJ(1));
t215 = g(1) * t182 - g(2) * t186;
t144 = t281 * t180;
t214 = -pkin(9) * t259 + qJD(6) * t144 + t180 * t229 + t274;
t145 = t281 * t184;
t72 = t184 * t82;
t213 = pkin(5) * t126 + qJD(6) * t145 - t180 * t286 + t72 + (-pkin(9) * t125 + t229) * t184;
t21 = t135 * t127 - t284 * t74;
t75 = t135 * t134;
t212 = -t120 * t21 + t75 * t79;
t210 = -t128 * t49 + t133 * t191;
t208 = -t128 * t107 - t133 * t41;
t110 = -t131 * t175 + t158;
t204 = -t110 * t175 + t111 * t177;
t203 = (-t143 * t175 + t160) * t175 - t119 * t177;
t202 = -t129 * t185 - t130 * t181;
t201 = -t122 * t293 + t184 * t81;
t54 = -qJD(3) * t133 + qJD(4) * t202;
t83 = pkin(4) * t128 + pkin(8) * t127;
t199 = t180 * t83 + t184 * t54 + t241 * t68 - t242 * t77;
t198 = -pkin(8) * t81 + t122 * t43;
t197 = t226 - t245;
t194 = -t236 + t285;
t190 = -t122 * t197 - t257 * t81;
t55 = qJD(3) * t134 + qJD(4) * t77;
t156 = -pkin(5) * t184 - pkin(4);
t115 = t163 * t248 + t253;
t114 = -t163 * t249 + t252;
t113 = -t163 * t252 + t249;
t112 = t163 * t253 + t248;
t87 = -qJD(4) * t127 + qJDD(4) * t134;
t86 = -qJD(4) * t128 - qJDD(4) * t133;
t73 = t184 * t83;
t65 = t184 * t68;
t56 = pkin(5) * t257 - t202;
t28 = pkin(5) * t197 + t55;
t27 = -pkin(9) * t257 + t273;
t26 = pkin(5) * t133 - pkin(9) * t256 - t180 * t77 + t65;
t11 = pkin(5) * t41 + t20;
t7 = -pkin(9) * t197 + t199;
t6 = pkin(9) * t258 + t128 * pkin(5) - t180 * t54 + t73 + (-t67 + (pkin(9) * t134 - t68) * t180) * qJD(5);
t4 = t12 * t183 - t15 * t179;
t1 = [qJDD(1), t215, g(1) * t186 + g(2) * t182 (t215 + (t176 ^ 2 + t178 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t194 * t177, -t194 * t175, t131 * t244 + t204 - t217, t139 * t153 - g(1) * (-pkin(1) * t182 - pkin(2) * t162 + qJ(3) * t164) - g(2) * (pkin(1) * t186 + pkin(2) * t164 + qJ(3) * t162) + t204 * t151 - t203 * qJD(3), -t126 * t127 + t134 * t84, -t125 * t127 - t126 * t128 - t133 * t84 - t134 * t85, t87, t86, 0, -t55 * qJD(4) + qJDD(4) * t202 + t121 * t133 + t124 * t128 + t142 * t85 + t163 * t216, -t54 * qJD(4) - t77 * qJDD(4) + t121 * t134 - t124 * t127 + t142 * t84 - t161 * t216, -t109 * t196 + t256 * t40 -(-t107 * t184 - t109 * t180) * t127 + (-t266 - t184 * t41 + (t107 * t180 - t109 * t184) * qJD(5)) * t134, t275 + t282, t190 + t208, t122 * t128 + t133 * t81 (-t77 * t241 + t73) * t122 + t65 * t81 + t218 * t133 + t24 * t128 + t55 * t107 - t202 * t41 + t43 * t226 - g(1) * t113 - g(2) * t115 + ((-qJD(5) * t68 - t54) * t122 - t77 * t81 + t223 * t133 + t20 * t134 - t43 * t127) * t180, -t199 * t122 - t273 * t81 - t200 * t133 - t25 * t128 + t55 * t109 - t202 * t40 - t43 * t258 - g(1) * t112 - g(2) * t114 + (t20 * t184 - t242 * t43) * t134, -t205 * t21 - t75 * t9, -t191 * t75 + t205 * t22 - t21 * t49 - t74 * t9, -t212 + t280, t210 + t276, t120 * t128 + t133 * t79 (-t179 * t7 + t183 * t6) * t120 + (-t179 * t27 + t183 * t26) * t79 + t230 * t133 + t4 * t128 + t28 * t49 - t56 * t191 + t11 * t74 + t32 * t22 - g(1) * t102 - g(2) * t104 + ((-t179 * t26 - t183 * t27) * t120 - t5 * t133) * qJD(6), -g(1) * t101 - g(2) * t103 - t11 * t75 - t5 * t128 + t13 * t133 + t32 * t21 - t28 * t205 + t56 * t9 + (-(-qJD(6) * t27 + t6) * t120 - t26 * t79 - t2 * t133) * t179 + (-(qJD(6) * t26 + t7) * t120 - t27 * t79 - t225 * t133) * t183; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t110 * t177 + t111 * t175 - g(3), 0, 0, 0, 0, 0, t86, -t87, 0, 0, 0, 0, 0, t190 - t208, t275 - t282, 0, 0, 0, 0, 0, -t210 + t276, t212 + t280; 0, 0, 0, 0, -t234, t235, -t244 * qJD(1) ^ 2, qJD(1) * t203 - t285, 0, 0, 0, 0, 0, 0.2e1 * t126 * qJD(4) + t211 (t125 - t228) * qJD(4) + t231, 0, 0, 0, 0, 0, t201 - t262, -t122 ^ 2 * t184 - t260 - t267, 0, 0, 0, 0, 0, t220 - t270, t265 + t283; 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t126, -t125 ^ 2 + t126 ^ 2 (-t125 - t228) * qJD(4) + t231, -t211, qJDD(4), -t124 * t126 + t193 + t207 + t296, -t124 * t125 + t163 * t217 - t209 + t279, t109 * t221 + t266 (t40 - t263) * t184 + (-t41 - t261) * t180, t122 * t221 - t260 + t267, t201 + t262, -t122 * t126, -pkin(4) * t41 - t47 * t107 - t72 * t122 - t24 * t126 + (t122 * t286 + t198) * t180 + t295 * t184, -pkin(4) * t40 - t47 * t109 + t122 * t274 + t25 * t126 - t180 * t295 + t184 * t198, t136 * t9 - t205 * t272, -t9 * t135 + t136 * t191 - t205 * t271 - t272 * t49, t265 - t283, t220 + t270, -t120 * t126 (-t144 * t183 - t145 * t179) * t79 - t156 * t191 + t11 * t135 - t4 * t126 + t219 * t49 - t271 * t32 + (t179 * t214 - t183 * t213) * t120 + t193 * t169 -(-t144 * t179 + t145 * t183) * t79 + t156 * t9 + t11 * t136 + t5 * t126 - t219 * t205 + t272 * t32 + (t179 * t213 + t183 * t214) * t120 - t193 * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t107, -t107 ^ 2 + t109 ^ 2, t40 + t263, t261 - t41, t81, -g(1) * t114 + g(2) * t112 - t109 * t43 + t122 * t25 + (t223 + t279) * t180 + t218, g(1) * t115 - g(2) * t113 + t107 * t43 + t122 * t24 + t184 * t279 - t200, -t294, t292, t291, t288, t79 -(-t14 * t179 - t268) * t120 + (-t109 * t49 - t120 * t240 + t183 * t79) * pkin(5) + t289 (-t120 * t15 - t2) * t179 + (t120 * t14 - t225) * t183 + (t109 * t205 - t120 * t239 - t179 * t79) * pkin(5) + t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, t292, t291, t288, t79, t120 * t5 + t289, t120 * t4 - t179 * t2 - t183 * t225 + t290;];
tau_reg  = t1;
