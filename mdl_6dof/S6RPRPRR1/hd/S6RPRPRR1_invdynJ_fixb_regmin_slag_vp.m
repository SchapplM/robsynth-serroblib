% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:23
% EndTime: 2019-03-09 03:35:31
% DurationCPUTime: 3.12s
% Computational Cost: add. (5131->337), mult. (11493->442), div. (0->0), fcn. (8855->16), ass. (0->204)
t182 = cos(qJ(6));
t239 = qJD(6) * t182;
t173 = sin(pkin(11));
t175 = cos(pkin(11));
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t133 = -t173 * t180 + t175 * t184;
t124 = t133 * qJD(1);
t183 = cos(qJ(5));
t114 = t183 * t124;
t134 = t173 * t184 + t175 * t180;
t126 = t134 * qJD(1);
t179 = sin(qJ(5));
t78 = -t179 * t126 + t114;
t295 = t182 * t78;
t303 = t239 - t295;
t166 = qJ(3) + pkin(11) + qJ(5);
t155 = sin(t166);
t170 = qJ(1) + pkin(10);
t163 = sin(t170);
t164 = cos(t170);
t214 = g(1) * t164 + g(2) * t163;
t302 = t214 * t155;
t169 = qJD(3) + qJD(5);
t255 = t78 * t169;
t241 = qJD(5) * t179;
t125 = t134 * qJD(3);
t85 = -qJD(1) * t125 + qJDD(1) * t133;
t238 = qJD(1) * qJD(3);
t229 = t184 * t238;
t230 = t180 * t238;
t86 = qJDD(1) * t134 - t173 * t230 + t175 * t229;
t36 = qJD(5) * t114 - t126 * t241 + t179 * t85 + t183 * t86;
t301 = t36 - t255;
t168 = qJDD(3) + qJDD(5);
t178 = sin(qJ(6));
t208 = t179 * t124 + t183 * t126;
t240 = qJD(6) * t178;
t14 = t178 * t168 + t169 * t239 + t182 * t36 - t208 * t240;
t72 = t178 * t169 + t182 * t208;
t15 = qJD(6) * t72 - t182 * t168 + t178 * t36;
t70 = -t182 * t169 + t178 * t208;
t300 = t14 * t182 - t178 * t15 - t303 * t70;
t12 = t14 * t178;
t299 = t303 * t72 + t12;
t265 = t72 * t208;
t293 = qJD(6) - t78;
t37 = qJD(5) * t208 + t179 * t86 - t183 * t85;
t34 = qJDD(6) + t37;
t31 = t178 * t34;
t73 = t293 * t239;
t298 = -t293 * t295 - t265 + t31 + t73;
t269 = t126 * pkin(8);
t174 = sin(pkin(10));
t157 = t174 * pkin(1) + pkin(7);
t246 = qJ(4) + t157;
t216 = t246 * qJD(1);
t102 = t180 * qJD(2) + t184 * t216;
t92 = t173 * t102;
t101 = t184 * qJD(2) - t180 * t216;
t260 = qJD(3) * pkin(3);
t96 = t101 + t260;
t59 = t175 * t96 - t92;
t46 = qJD(3) * pkin(4) - t269 + t59;
t270 = t124 * pkin(8);
t247 = t175 * t102;
t60 = t173 * t96 + t247;
t47 = t60 + t270;
t24 = -t179 * t47 + t183 * t46;
t20 = -t169 * pkin(5) - t24;
t297 = t20 * t78;
t156 = cos(t166);
t272 = g(3) * t156;
t165 = t184 * qJDD(2);
t142 = t157 * qJDD(1);
t195 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t142;
t204 = t216 * qJD(3);
t57 = qJDD(3) * pkin(3) - t180 * t195 - t184 * t204 + t165;
t61 = (qJDD(2) - t204) * t180 + t195 * t184;
t28 = -t173 * t61 + t175 * t57;
t22 = qJDD(3) * pkin(4) - t86 * pkin(8) + t28;
t29 = t173 * t57 + t175 * t61;
t23 = t85 * pkin(8) + t29;
t25 = t179 * t46 + t183 * t47;
t279 = qJD(5) * t25 + t179 * t23 - t183 * t22;
t3 = -t168 * pkin(5) + t279;
t296 = t3 + t272;
t294 = t208 * t78;
t256 = t208 * t169;
t291 = -t37 + t256;
t289 = t208 ^ 2 - t78 ^ 2;
t149 = g(3) * t155;
t277 = (qJD(5) * t46 + t23) * t183 + t179 * t22 - t47 * t241;
t160 = t184 * pkin(3) + pkin(2);
t176 = cos(pkin(10));
t268 = t176 * pkin(1);
t206 = -t160 - t268;
t122 = qJD(1) * t206 + qJD(4);
t87 = -t124 * pkin(4) + t122;
t288 = t156 * t214 - t87 * t78 + t149 - t277;
t286 = pkin(5) * t208;
t264 = t208 * t70;
t103 = t180 * qJD(1) * pkin(3) + t126 * pkin(4);
t158 = t175 * pkin(3) + pkin(4);
t275 = pkin(3) * t173;
t243 = t179 * t158 + t183 * t275;
t121 = pkin(9) + t243;
t285 = (-t78 * pkin(9) + qJD(6) * t121 + t103 + t286) * t293;
t284 = (pkin(9) * t293 + t286) * t293;
t283 = t293 * t208;
t21 = t169 * pkin(9) + t25;
t38 = -pkin(5) * t78 - pkin(9) * t208 + t87;
t6 = -t178 * t21 + t182 * t38;
t282 = t182 * t302 + t20 * t240 - t6 * t208;
t7 = t178 * t38 + t182 * t21;
t281 = t178 * t296 + t20 * t239 + t7 * t208;
t280 = -t208 * t87 - t272 - t279 + t302;
t128 = t133 * qJD(3);
t207 = t183 * t133 - t179 * t134;
t48 = qJD(5) * t207 - t179 * t125 + t183 * t128;
t89 = t179 * t133 + t183 * t134;
t210 = t293 * t48 + t34 * t89;
t232 = t89 * t240;
t278 = -t182 * t210 + t232 * t293;
t194 = pkin(3) * t230 + qJDD(1) * t206 + qJDD(4);
t227 = t168 * pkin(9) + qJD(6) * t38 + t277;
t131 = t246 * t180;
t132 = t246 * t184;
t83 = -t175 * t131 - t173 * t132;
t68 = -t134 * pkin(8) + t83;
t84 = -t173 * t131 + t175 * t132;
t69 = t133 * pkin(8) + t84;
t40 = t179 * t68 + t183 * t69;
t100 = -t133 * pkin(4) + t206;
t42 = -pkin(5) * t207 - t89 * pkin(9) + t100;
t39 = t179 * t69 - t183 * t68;
t217 = qJD(3) * t246;
t109 = t184 * qJD(4) - t180 * t217;
t110 = -t180 * qJD(4) - t184 * t217;
t66 = -t173 * t109 + t175 * t110;
t52 = -t128 * pkin(8) + t66;
t67 = t175 * t109 + t173 * t110;
t53 = -t125 * pkin(8) + t67;
t8 = -qJD(5) * t39 + t179 * t52 + t183 * t53;
t276 = t20 * t48 - (qJD(6) * t42 + t8) * t293 + t227 * t207 + t3 * t89 - t40 * t34;
t271 = g(3) * t184;
t267 = t20 * t89;
t266 = t42 * t34;
t49 = qJD(5) * t89 + t183 * t125 + t179 * t128;
t263 = -t14 * t207 + t72 * t49;
t63 = t175 * t101 - t92;
t259 = t178 * t72;
t202 = t183 * t158 - t179 * t275;
t62 = -t173 * t101 - t247;
t50 = t62 - t270;
t51 = t63 - t269;
t254 = -t202 * qJD(5) + t179 * t50 + t183 * t51;
t253 = t243 * qJD(5) - t179 * t51 + t183 * t50;
t251 = t163 * t178;
t250 = t163 * t182;
t249 = t164 * t178;
t248 = t164 * t182;
t244 = qJDD(2) - g(3);
t171 = t180 ^ 2;
t242 = -t184 ^ 2 + t171;
t159 = -pkin(2) - t268;
t145 = qJD(1) * t159;
t237 = t184 * qJDD(1);
t162 = t180 * t260;
t104 = t125 * pkin(4) + t162;
t58 = -t85 * pkin(4) + t194;
t5 = t37 * pkin(5) - t36 * pkin(9) + t58;
t226 = qJD(6) * t21 - t5;
t219 = t178 * t293;
t213 = g(1) * t163 - g(2) * t164;
t181 = sin(qJ(1));
t185 = cos(qJ(1));
t212 = g(1) * t181 - g(2) * t185;
t211 = t15 * t207 - t49 * t70;
t209 = t89 * t168 + t48 * t169;
t32 = t182 * t34;
t205 = t32 - (-t178 * t78 + t240) * t293;
t203 = -t227 + t149;
t200 = -pkin(9) * t34 + t24 * t293 - t297;
t198 = -t145 * qJD(1) - t142 + t214;
t197 = -t121 * t34 + t254 * t293 - t297;
t196 = 0.2e1 * qJD(3) * t145 - qJDD(3) * t157;
t186 = qJD(3) ^ 2;
t192 = -0.2e1 * qJDD(1) * t159 - t157 * t186 + t213;
t191 = -t178 * t210 - t73 * t89;
t187 = qJD(1) ^ 2;
t177 = -qJ(4) - pkin(7);
t140 = qJDD(3) * t184 - t186 * t180;
t139 = qJDD(3) * t180 + t186 * t184;
t120 = -pkin(5) - t202;
t108 = t156 * t248 + t251;
t107 = -t156 * t249 + t250;
t106 = -t156 * t250 + t249;
t105 = t156 * t251 + t248;
t35 = t168 * t207 - t49 * t169;
t16 = t49 * pkin(5) - t48 * pkin(9) + t104;
t9 = qJD(5) * t40 + t179 * t53 - t183 * t52;
t4 = t182 * t5;
t1 = [qJDD(1), t212, g(1) * t185 + g(2) * t181 (t212 + (t174 ^ 2 + t176 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t171 * qJDD(1) + 0.2e1 * t180 * t229, 0.2e1 * t180 * t237 - 0.2e1 * t238 * t242, t139, t140, 0, t180 * t196 + t184 * t192, -t180 * t192 + t184 * t196, t67 * t124 - t60 * t125 - t66 * t126 - t59 * t128 + t29 * t133 - t28 * t134 - t83 * t86 + t84 * t85 - t214, t29 * t84 + t60 * t67 + t28 * t83 + t59 * t66 + t122 * t162 - g(1) * (-t181 * pkin(1) - t163 * t160 - t164 * t177) - g(2) * (t185 * pkin(1) + t164 * t160 - t163 * t177) + t194 * t206, t208 * t48 + t36 * t89, t207 * t36 - t208 * t49 - t89 * t37 + t48 * t78, t209, t35, 0, t100 * t37 - t104 * t78 + t156 * t213 - t39 * t168 - t9 * t169 - t207 * t58 + t87 * t49, t100 * t36 + t104 * t208 - t155 * t213 - t40 * t168 - t8 * t169 + t87 * t48 + t58 * t89, -t72 * t232 + (t14 * t89 + t48 * t72) * t182 (-t182 * t70 - t259) * t48 + (-t12 - t15 * t182 + (t178 * t70 - t182 * t72) * qJD(6)) * t89, t263 - t278, t191 + t211, -t207 * t34 + t293 * t49, -g(1) * t106 - g(2) * t108 + t39 * t15 - t4 * t207 + t6 * t49 + t9 * t70 + (t16 * t293 + t266 + (t207 * t21 - t293 * t40 + t267) * qJD(6)) * t182 + t276 * t178, -g(1) * t105 - g(2) * t107 + t39 * t14 - t7 * t49 + t9 * t72 + (-(-qJD(6) * t40 + t16) * t293 - t266 - t226 * t207 - qJD(6) * t267) * t178 + t276 * t182; 0, 0, 0, t244, 0, 0, 0, 0, 0, t140, -t139, t128 * t124 + t125 * t126 - t133 * t86 + t134 * t85, -t59 * t125 + t60 * t128 + t28 * t133 + t29 * t134 - g(3), 0, 0, 0, 0, 0, t35, -t209, 0, 0, 0, 0, 0, t191 - t211, t263 + t278; 0, 0, 0, 0, -t180 * t187 * t184, t242 * t187, t180 * qJDD(1), t237, qJDD(3), t180 * t198 + t165 - t271, -t180 * t244 + t184 * t198 (t60 + t62) * t126 + (t59 - t63) * t124 + (t173 * t85 - t175 * t86) * pkin(3), -t59 * t62 - t60 * t63 + (-t271 + t173 * t29 + t175 * t28 + (-qJD(1) * t122 + t214) * t180) * pkin(3), -t294, t289, t301, t291, t168, t103 * t78 + t168 * t202 - t169 * t253 + t280, -t103 * t208 - t168 * t243 + t169 * t254 + t288, t299, -t259 * t293 + t300, t298, t205 + t264, -t283, t120 * t15 + t253 * t70 + (-t296 - t285) * t182 + t197 * t178 + t282, t120 * t14 + t253 * t72 + t197 * t182 + (-t302 + t285) * t178 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 ^ 2 - t126 ^ 2, -t60 * t124 + t59 * t126 + t194 - t213, 0, 0, 0, 0, 0, t37 + t256, t36 + t255, 0, 0, 0, 0, 0, t205 - t264, -t182 * t293 ^ 2 - t265 - t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, t289, t301, t291, t168, t25 * t169 + t280, t24 * t169 + t288, t299, -t219 * t72 + t300, t298, -t219 * t293 + t264 + t32, -t283, -pkin(5) * t15 - t25 * t70 + t200 * t178 + (-t296 - t284) * t182 + t282, -pkin(5) * t14 - t25 * t72 + t200 * t182 + (-t302 + t284) * t178 + t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t70, -t70 ^ 2 + t72 ^ 2, t293 * t70 + t14, t293 * t72 - t15, t34, -g(1) * t107 + g(2) * t105 + t178 * t203 - t20 * t72 - t21 * t239 + t293 * t7 + t4, g(1) * t108 - g(2) * t106 + t178 * t226 + t182 * t203 + t20 * t70 + t293 * t6;];
tau_reg  = t1;
