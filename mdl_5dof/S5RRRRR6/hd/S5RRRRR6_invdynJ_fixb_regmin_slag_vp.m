% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:21
% EndTime: 2020-01-03 12:15:28
% DurationCPUTime: 2.44s
% Computational Cost: add. (3550->296), mult. (5394->387), div. (0->0), fcn. (3946->16), ass. (0->209)
t176 = qJ(3) + qJ(4);
t167 = qJ(5) + t176;
t153 = cos(t167);
t178 = sin(qJ(5));
t183 = cos(qJ(5));
t173 = qJD(1) + qJD(2);
t184 = cos(qJ(4));
t185 = cos(qJ(3));
t259 = t184 * t185;
t238 = t173 * t259;
t179 = sin(qJ(4));
t180 = sin(qJ(3));
t264 = t179 * t180;
t239 = t173 * t264;
t92 = -t238 + t239;
t260 = t180 * t184;
t114 = t179 * t185 + t260;
t94 = t114 * t173;
t210 = t178 * t92 - t183 * t94;
t152 = sin(t167);
t177 = qJ(1) + qJ(2);
t166 = cos(t177);
t270 = t152 * t166;
t164 = sin(t177);
t271 = t152 * t164;
t170 = qJDD(3) + qJDD(4);
t181 = sin(qJ(2));
t287 = pkin(1) * t181;
t242 = qJD(1) * t287;
t125 = pkin(7) * t173 + t242;
t233 = pkin(8) * t173 + t125;
t84 = t233 * t185;
t77 = t184 * t84;
t83 = t233 * t180;
t78 = qJD(3) * pkin(3) - t83;
t209 = -t179 * t78 - t77;
t249 = qJD(3) * t185;
t234 = t173 * t249;
t171 = qJDD(1) + qJDD(2);
t261 = t180 * t171;
t244 = qJDD(1) * t181;
t186 = cos(qJ(2));
t251 = qJD(2) * t186;
t97 = t171 * pkin(7) + (qJD(1) * t251 + t244) * pkin(1);
t41 = -t125 * t249 + qJDD(3) * pkin(3) - t180 * t97 + (-t234 - t261) * pkin(8);
t250 = qJD(3) * t180;
t235 = t173 * t250;
t258 = t185 * t171;
t42 = -t125 * t250 + t185 * t97 + (-t235 + t258) * pkin(8);
t198 = qJD(4) * t209 - t179 * t42 + t184 * t41;
t172 = qJD(3) + qJD(4);
t37 = qJD(4) * t238 + t171 * t260 - t172 * t239 + t179 * t258 + t184 * t234;
t5 = pkin(4) * t170 - pkin(9) * t37 + t198;
t157 = -pkin(3) * t185 - pkin(2);
t253 = qJD(1) * t186;
t240 = pkin(1) * t253;
t95 = t157 * t173 - t240;
t57 = t92 * pkin(4) + t95;
t248 = qJD(4) * t179;
t291 = (qJD(4) * t78 + t42) * t184 + t179 * t41 - t84 * t248;
t214 = t179 * t261 - t184 * t258;
t67 = t172 * t114;
t38 = t173 * t67 + t214;
t6 = -t38 * pkin(9) + t291;
t299 = -g(1) * t153 + g(2) * t271 - g(3) * t270 - t178 * t6 + t183 * t5 + t57 * t210;
t232 = g(2) * t164 - g(3) * t166;
t246 = qJD(5) * t178;
t289 = pkin(9) * t92;
t31 = -t209 - t289;
t54 = -t178 * t94 - t183 * t92;
t298 = g(1) * t152 + t232 * t153 + t31 * t246 - t57 * t54;
t215 = -g(2) * t166 - g(3) * t164;
t252 = qJD(2) * t181;
t160 = pkin(1) * t252;
t286 = pkin(1) * t186;
t255 = -qJD(1) * t160 + qJDD(1) * t286;
t285 = pkin(2) * t171;
t96 = -t255 - t285;
t297 = t96 - t215;
t272 = t183 * t31;
t75 = t179 * t84;
t228 = t184 * t78 - t75;
t90 = t94 * pkin(9);
t30 = t228 - t90;
t28 = pkin(4) * t172 + t30;
t213 = -t178 * t28 - t272;
t296 = qJD(5) * t213 + t299;
t243 = -qJD(4) - qJD(5);
t162 = qJD(3) - t243;
t295 = (-t162 * t31 - t5) * t178 + t298;
t278 = t210 * t54;
t9 = t210 ^ 2 - t54 ^ 2;
t245 = qJD(5) * t183;
t10 = -t178 * t38 + t183 * t37 - t92 * t245 - t94 * t246;
t7 = -t162 * t54 + t10;
t197 = qJD(5) * t210 - t178 * t37 - t183 * t38;
t8 = -t162 * t210 + t197;
t113 = -t259 + t264;
t290 = -pkin(7) - pkin(8);
t237 = qJD(3) * t290;
t118 = t180 * t237;
t119 = t185 * t237;
t142 = t290 * t180;
t168 = t185 * pkin(8);
t143 = pkin(7) * t185 + t168;
t247 = qJD(4) * t184;
t294 = -t113 * t240 - t184 * t118 - t179 * t119 - t142 * t247 + t143 * t248;
t256 = t179 * t142 + t184 * t143;
t293 = -t256 * qJD(4) + t114 * t240 - t179 * t118 + t184 * t119;
t154 = pkin(7) + t287;
t274 = -pkin(8) - t154;
t110 = t274 * t180;
t111 = t154 * t185 + t168;
t257 = t179 * t110 + t184 * t111;
t188 = qJD(3) ^ 2;
t292 = pkin(7) * t188 - t285;
t66 = t172 * t113;
t288 = t66 * pkin(9);
t284 = pkin(2) * t173;
t282 = pkin(9) * t114;
t275 = t94 * t92;
t273 = -t184 * t83 - t75;
t161 = qJDD(5) + t170;
t269 = t161 * t179;
t268 = t161 * t183;
t163 = sin(t176);
t267 = t163 * t164;
t266 = t163 * t166;
t265 = t173 * t180;
t263 = t179 * t183;
t174 = t180 ^ 2;
t254 = -t185 ^ 2 + t174;
t241 = pkin(1) * t251;
t159 = pkin(3) * t250;
t236 = t173 * t252;
t60 = pkin(4) * t67 + t159;
t231 = -qJD(5) * t28 - t6;
t227 = t179 * t83 - t77;
t226 = qJD(3) * t274;
t224 = t184 * t110 - t111 * t179;
t223 = t184 * t142 - t143 * t179;
t63 = t183 * t113 + t114 * t178;
t21 = -qJD(5) * t63 - t178 * t67 - t183 * t66;
t61 = pkin(3) * t235 + t157 * t171 - t255;
t23 = pkin(4) * t38 + t61;
t64 = -t113 * t178 + t114 * t183;
t222 = g(2) * t270 + g(3) * t271 + t57 * t21 + t23 * t64;
t221 = g(2) * t266 + g(3) * t267 + t61 * t114 - t95 * t66;
t220 = t173 * t242;
t126 = -t240 - t284;
t219 = t126 * t249 + t297 * t180;
t218 = t60 - t242;
t58 = t223 - t282;
t65 = t67 * pkin(9);
t217 = -qJD(5) * t58 + t294 + t65;
t109 = t113 * pkin(9);
t59 = -t109 + t256;
t216 = qJD(5) * t59 - t288 - t293;
t47 = t224 - t282;
t48 = -t109 + t257;
t212 = -t178 * t48 + t183 * t47;
t211 = t178 * t47 + t183 * t48;
t91 = pkin(4) * t113 + t157;
t79 = t180 * t226 + t185 * t241;
t80 = -t180 * t241 + t185 * t226;
t207 = t110 * t247 - t111 * t248 + t179 * t80 + t184 * t79;
t205 = -t159 + t242;
t204 = -t126 * t173 + t232 - t97;
t22 = qJD(5) * t64 - t178 * t66 + t183 * t67;
t203 = t153 * t215 + t57 * t22 + t23 * t63;
t165 = cos(t176);
t202 = t61 * t113 + t165 * t215 + t95 * t67;
t201 = t215 + t220;
t156 = -pkin(2) - t286;
t200 = pkin(1) * t236 + t154 * t188 + t156 * t171;
t196 = -pkin(7) * qJDD(3) + (t240 - t284) * qJD(3);
t195 = -t257 * qJD(4) - t179 * t79 + t184 * t80;
t193 = -qJDD(3) * t154 + (t156 * t173 - t241) * qJD(3);
t191 = g(1) * t163 + t232 * t165 + t95 * t92 - t291;
t189 = -g(1) * t165 + g(2) * t267 - g(3) * t266 - t95 * t94 + t198;
t187 = cos(qJ(1));
t182 = sin(qJ(1));
t169 = t173 ^ 2;
t155 = pkin(3) * t184 + pkin(4);
t135 = t157 - t286;
t134 = qJDD(3) * t185 - t180 * t188;
t133 = qJDD(3) * t180 + t185 * t188;
t120 = t160 + t159;
t103 = t126 * t250;
t98 = t171 * t174 + 0.2e1 * t180 * t234;
t82 = t91 - t286;
t69 = -0.2e1 * t254 * t173 * qJD(3) + 0.2e1 * t180 * t258;
t68 = pkin(3) * t265 + pkin(4) * t94;
t56 = t160 + t60;
t44 = -t113 * t170 - t172 * t67;
t43 = t114 * t170 - t172 * t66;
t40 = -t92 ^ 2 + t94 ^ 2;
t33 = -t90 + t273;
t32 = t227 + t289;
t26 = t172 * t92 + t37;
t20 = t195 + t288;
t19 = t207 - t65;
t18 = t114 * t37 - t66 * t94;
t13 = -t161 * t63 - t162 * t22;
t12 = t161 * t64 + t162 * t21;
t3 = -t113 * t37 - t114 * t38 + t66 * t92 - t67 * t94;
t2 = t10 * t64 - t21 * t210;
t1 = -t10 * t63 + t197 * t64 + t21 * t54 + t210 * t22;
t4 = [qJDD(1), -g(2) * t187 - g(3) * t182, g(2) * t182 - g(3) * t187, t171, (t171 * t186 - t236) * pkin(1) + t215 + t255, ((-qJDD(1) - t171) * t181 + (-qJD(1) - t173) * t251) * pkin(1) + t232, t98, t69, t133, t134, 0, t103 + t193 * t180 + (-t200 - t297) * t185, t180 * t200 + t185 * t193 + t219, t18, t3, t43, t44, 0, t120 * t92 + t135 * t38 + t170 * t224 + t172 * t195 + t202, t120 * t94 + t135 * t37 - t170 * t257 - t172 * t207 + t221, t2, t1, t12, t13, 0, -t56 * t54 - t82 * t197 + (-qJD(5) * t211 - t178 * t19 + t183 * t20) * t162 + t212 * t161 + t203, -t56 * t210 + t82 * t10 - (qJD(5) * t212 + t178 * t20 + t183 * t19) * t162 - t211 * t161 + t222; 0, 0, 0, t171, t201 + t255, (-t244 + (-qJD(2) + t173) * t253) * pkin(1) + t232, t98, t69, t133, t134, 0, t103 + t196 * t180 + (t201 - t96 - t292) * t185, t196 * t185 + (-t220 + t292) * t180 + t219, t18, t3, t43, t44, 0, t157 * t38 + t223 * t170 + t293 * t172 - t205 * t92 + t202, t157 * t37 - t256 * t170 + t294 * t172 - t205 * t94 + t221, t2, t1, t12, t13, 0, -t91 * t197 + (-t178 * t59 + t183 * t58) * t161 - t218 * t54 + (t178 * t217 - t183 * t216) * t162 + t203, t91 * t10 - (t178 * t58 + t183 * t59) * t161 - t218 * t210 + (t178 * t216 + t183 * t217) * t162 + t222; 0, 0, 0, 0, 0, 0, -t180 * t169 * t185, t254 * t169, t261, t258, qJDD(3), -g(1) * t185 + t180 * t204, g(1) * t180 + t185 * t204, t275, t40, t26, -t214, t170, -t227 * t172 + (t170 * t184 - t172 * t248 - t265 * t92) * pkin(3) + t189, t273 * t172 + (-t170 * t179 - t172 * t247 - t265 * t94) * pkin(3) + t191, t278, t9, t7, t8, t161, t155 * t268 + t68 * t54 - (-t178 * t33 + t183 * t32) * t162 + (-t178 * t269 + (-t178 * t184 - t263) * t162 * qJD(4)) * pkin(3) + ((-pkin(3) * t263 - t155 * t178) * t162 + t213) * qJD(5) + t299, t68 * t210 + (-t155 * t161 - t5 + (-pkin(3) * t179 * t243 + t32) * t162) * t178 + (-pkin(3) * t269 + (-pkin(3) * t247 - qJD(5) * t155 + t33) * t162 + t231) * t183 + t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t40, t26, -t214, t170, -t172 * t209 + t189, t172 * t228 + t191, t278, t9, t7, t8, t161, -(-t178 * t30 - t272) * t162 + (-t162 * t246 + t54 * t94 + t268) * pkin(4) + t296, (t162 * t30 + t231) * t183 + (-t161 * t178 - t162 * t245 + t210 * t94) * pkin(4) + t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, t9, t7, t8, t161, -t162 * t213 + t296, (-t6 + (-qJD(5) + t162) * t28) * t183 + t295;];
tau_reg = t4;
