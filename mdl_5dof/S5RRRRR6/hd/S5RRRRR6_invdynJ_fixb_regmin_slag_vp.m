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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:08:22
% EndTime: 2022-01-20 12:08:31
% DurationCPUTime: 2.53s
% Computational Cost: add. (3550->300), mult. (5394->395), div. (0->0), fcn. (3946->16), ass. (0->214)
t182 = qJ(3) + qJ(4);
t173 = qJ(5) + t182;
t159 = cos(t173);
t184 = sin(qJ(5));
t189 = cos(qJ(5));
t179 = qJD(1) + qJD(2);
t190 = cos(qJ(4));
t191 = cos(qJ(3));
t267 = t190 * t191;
t243 = t179 * t267;
t185 = sin(qJ(4));
t186 = sin(qJ(3));
t273 = t185 * t186;
t244 = t179 * t273;
t92 = -t243 + t244;
t269 = t186 * t190;
t114 = t185 * t191 + t269;
t94 = t114 * t179;
t218 = t184 * t92 - t189 * t94;
t158 = sin(t173);
t183 = qJ(1) + qJ(2);
t172 = cos(t183);
t282 = t158 * t172;
t170 = sin(t183);
t283 = t158 * t170;
t176 = qJDD(3) + qJDD(4);
t187 = sin(qJ(2));
t297 = pkin(1) * t187;
t248 = qJD(1) * t297;
t127 = pkin(7) * t179 + t248;
t237 = pkin(8) * t179 + t127;
t84 = t237 * t191;
t77 = t190 * t84;
t83 = t237 * t186;
t78 = qJD(3) * pkin(3) - t83;
t217 = -t185 * t78 - t77;
t255 = qJD(3) * t191;
t239 = t179 * t255;
t177 = qJDD(1) + qJDD(2);
t270 = t186 * t177;
t250 = qJDD(1) * t187;
t192 = cos(qJ(2));
t257 = qJD(2) * t192;
t97 = pkin(7) * t177 + (qJD(1) * t257 + t250) * pkin(1);
t41 = -t127 * t255 + qJDD(3) * pkin(3) - t186 * t97 + (-t239 - t270) * pkin(8);
t256 = qJD(3) * t186;
t240 = t179 * t256;
t266 = t191 * t177;
t42 = -t127 * t256 + t191 * t97 + (-t240 + t266) * pkin(8);
t203 = t217 * qJD(4) - t185 * t42 + t190 * t41;
t178 = qJD(3) + qJD(4);
t37 = qJD(4) * t243 + t177 * t269 - t178 * t244 + t185 * t266 + t190 * t239;
t5 = pkin(4) * t176 - pkin(9) * t37 + t203;
t163 = -pkin(3) * t191 - pkin(2);
t259 = qJD(1) * t192;
t247 = pkin(1) * t259;
t95 = t163 * t179 - t247;
t57 = pkin(4) * t92 + t95;
t254 = qJD(4) * t185;
t301 = (qJD(4) * t78 + t42) * t190 + t185 * t41 - t84 * t254;
t222 = t185 * t270 - t190 * t266;
t67 = t178 * t114;
t38 = t67 * t179 + t222;
t6 = -pkin(9) * t38 + t301;
t309 = g(1) * t282 + g(2) * t283 - g(3) * t159 - t184 * t6 + t189 * t5 + t57 * t218;
t252 = qJD(5) * t184;
t280 = t159 * t172;
t281 = t159 * t170;
t298 = pkin(9) * t92;
t31 = -t217 - t298;
t54 = -t184 * t94 - t189 * t92;
t308 = g(1) * t280 + g(2) * t281 + g(3) * t158 + t31 * t252 - t57 * t54;
t75 = t185 * t84;
t233 = t190 * t78 - t75;
t90 = t94 * pkin(9);
t30 = t233 - t90;
t28 = pkin(4) * t178 + t30;
t284 = t189 * t31;
t221 = -t184 * t28 - t284;
t307 = t221 * qJD(5) + t309;
t249 = -qJD(4) - qJD(5);
t168 = qJD(3) - t249;
t306 = (-t31 * t168 - t5) * t184 + t308;
t258 = qJD(2) * t187;
t166 = pkin(1) * t258;
t296 = pkin(1) * t192;
t261 = -qJD(1) * t166 + qJDD(1) * t296;
t305 = -g(2) * t172 + t261;
t290 = t218 * t54;
t295 = pkin(2) * t177;
t304 = -t295 - t305;
t9 = t218 ^ 2 - t54 ^ 2;
t251 = qJD(5) * t189;
t10 = -t184 * t38 + t189 * t37 - t92 * t251 - t94 * t252;
t7 = -t168 * t54 + t10;
t202 = t218 * qJD(5) - t184 * t37 - t189 * t38;
t8 = -t168 * t218 + t202;
t113 = -t267 + t273;
t300 = -pkin(7) - pkin(8);
t242 = qJD(3) * t300;
t118 = t186 * t242;
t119 = t191 * t242;
t146 = t300 * t186;
t174 = t191 * pkin(8);
t147 = pkin(7) * t191 + t174;
t253 = qJD(4) * t190;
t303 = -t113 * t247 - t190 * t118 - t185 * t119 - t146 * t253 + t147 * t254;
t263 = t185 * t146 + t190 * t147;
t302 = -t263 * qJD(4) + t114 * t247 - t185 * t118 + t190 * t119;
t160 = pkin(7) + t297;
t286 = -pkin(8) - t160;
t110 = t286 * t186;
t111 = t160 * t191 + t174;
t265 = t185 * t110 + t190 * t111;
t262 = g(1) * t172 + g(2) * t170;
t66 = t178 * t113;
t299 = pkin(9) * t66;
t294 = pkin(2) * t179;
t293 = pkin(3) * t185;
t292 = pkin(9) * t114;
t155 = g(1) * t170;
t287 = t94 * t92;
t285 = -t190 * t83 - t75;
t169 = sin(t182);
t279 = t169 * t170;
t278 = t169 * t172;
t171 = cos(t182);
t277 = t170 * t171;
t276 = t171 * t172;
t275 = t179 * t186;
t167 = qJDD(5) + t176;
t274 = t184 * t167;
t272 = t185 * t189;
t268 = t189 * t167;
t128 = -t247 - t294;
t264 = t128 * t256 + t191 * t155;
t180 = t186 ^ 2;
t260 = -t191 ^ 2 + t180;
t246 = pkin(1) * t257;
t165 = pkin(3) * t256;
t245 = t128 * t255 + t304 * t186;
t241 = t179 * t258;
t60 = pkin(4) * t67 + t165;
t236 = -qJD(5) * t28 - t6;
t232 = t185 * t83 - t77;
t231 = qJD(3) * t286;
t229 = t190 * t110 - t111 * t185;
t228 = t190 * t146 - t147 * t185;
t227 = t179 * t248;
t225 = t60 - t248;
t58 = t228 - t292;
t65 = t67 * pkin(9);
t224 = -qJD(5) * t58 + t303 + t65;
t109 = t113 * pkin(9);
t59 = -t109 + t263;
t223 = qJD(5) * t59 - t299 - t302;
t47 = t229 - t292;
t48 = -t109 + t265;
t220 = -t184 * t48 + t189 * t47;
t219 = t184 * t47 + t189 * t48;
t63 = t189 * t113 + t114 * t184;
t64 = -t113 * t184 + t114 * t189;
t91 = pkin(4) * t113 + t163;
t216 = t155 + t305;
t21 = -t63 * qJD(5) - t184 * t67 - t189 * t66;
t61 = pkin(3) * t240 + t163 * t177 - t261;
t23 = pkin(4) * t38 + t61;
t215 = -g(1) * t283 + g(2) * t282 + t57 * t21 + t23 * t64;
t214 = -g(1) * t279 + g(2) * t278 + t61 * t114 - t95 * t66;
t22 = t64 * qJD(5) - t184 * t66 + t189 * t67;
t213 = g(1) * t281 - g(2) * t280 + t57 * t22 + t23 * t63;
t212 = g(1) * t277 - g(2) * t276 + t61 * t113 + t95 * t67;
t79 = t186 * t231 + t191 * t246;
t80 = -t186 * t246 + t191 * t231;
t211 = t110 * t253 - t111 * t254 + t185 * t80 + t190 * t79;
t209 = -t165 + t248;
t207 = -t128 * t179 + t262 - t97;
t194 = qJD(3) ^ 2;
t206 = pkin(7) * t194 - t227 - t295;
t162 = -pkin(2) - t296;
t204 = pkin(1) * t241 + t160 * t194 + t162 * t177;
t201 = -pkin(7) * qJDD(3) + (t247 - t294) * qJD(3);
t200 = -t265 * qJD(4) - t185 * t79 + t190 * t80;
t198 = -qJDD(3) * t160 + (t162 * t179 - t246) * qJD(3);
t197 = g(1) * t276 + g(2) * t277 + g(3) * t169 + t95 * t92 - t301;
t195 = g(1) * t278 + g(2) * t279 - g(3) * t171 - t95 * t94 + t203;
t193 = cos(qJ(1));
t188 = sin(qJ(1));
t175 = t179 ^ 2;
t161 = pkin(3) * t190 + pkin(4);
t139 = t163 - t296;
t138 = qJDD(3) * t191 - t186 * t194;
t137 = qJDD(3) * t186 + t191 * t194;
t120 = t166 + t165;
t98 = t177 * t180 + 0.2e1 * t186 * t239;
t82 = t91 - t296;
t69 = -0.2e1 * t260 * t179 * qJD(3) + 0.2e1 * t186 * t266;
t68 = pkin(3) * t275 + pkin(4) * t94;
t56 = t166 + t60;
t44 = -t113 * t176 - t178 * t67;
t43 = t114 * t176 - t178 * t66;
t40 = -t92 ^ 2 + t94 ^ 2;
t33 = -t90 + t285;
t32 = t232 + t298;
t26 = t178 * t92 + t37;
t20 = t200 + t299;
t19 = t211 - t65;
t18 = t114 * t37 - t66 * t94;
t13 = -t167 * t63 - t168 * t22;
t12 = t167 * t64 + t168 * t21;
t3 = -t113 * t37 - t114 * t38 + t66 * t92 - t67 * t94;
t2 = t10 * t64 - t21 * t218;
t1 = -t10 * t63 + t202 * t64 + t21 * t54 + t218 * t22;
t4 = [qJDD(1), g(1) * t188 - g(2) * t193, g(1) * t193 + g(2) * t188, t177, (t177 * t192 - t241) * pkin(1) + t216, ((-qJDD(1) - t177) * t187 + (-qJD(1) - t179) * t257) * pkin(1) + t262, t98, t69, t137, t138, 0, t198 * t186 + (-t204 - t304) * t191 + t264, t198 * t191 + (t204 - t155) * t186 + t245, t18, t3, t43, t44, 0, t120 * t92 + t139 * t38 + t176 * t229 + t178 * t200 + t212, t120 * t94 + t139 * t37 - t176 * t265 - t178 * t211 + t214, t2, t1, t12, t13, 0, -t56 * t54 - t82 * t202 + (-qJD(5) * t219 - t184 * t19 + t189 * t20) * t168 + t220 * t167 + t213, -t56 * t218 + t82 * t10 - (qJD(5) * t220 + t184 * t20 + t189 * t19) * t168 - t219 * t167 + t215; 0, 0, 0, t177, t216 + t227, (-t250 + (-qJD(2) + t179) * t259) * pkin(1) + t262, t98, t69, t137, t138, 0, t201 * t186 + (-t206 - t304) * t191 + t264, t201 * t191 + (t206 - t155) * t186 + t245, t18, t3, t43, t44, 0, t163 * t38 + t228 * t176 + t302 * t178 - t209 * t92 + t212, t163 * t37 - t263 * t176 + t303 * t178 - t209 * t94 + t214, t2, t1, t12, t13, 0, -t91 * t202 + (-t184 * t59 + t189 * t58) * t167 - t225 * t54 + (t184 * t224 - t189 * t223) * t168 + t213, t91 * t10 - (t184 * t58 + t189 * t59) * t167 - t225 * t218 + (t184 * t223 + t189 * t224) * t168 + t215; 0, 0, 0, 0, 0, 0, -t186 * t175 * t191, t260 * t175, t270, t266, qJDD(3), -g(3) * t191 + t186 * t207, g(3) * t186 + t191 * t207, t287, t40, t26, -t222, t176, -t232 * t178 + (t176 * t190 - t178 * t254 - t275 * t92) * pkin(3) + t195, t285 * t178 + (-t176 * t185 - t178 * t253 - t275 * t94) * pkin(3) + t197, t290, t9, t7, t8, t167, t161 * t268 + t68 * t54 - (-t184 * t33 + t189 * t32) * t168 + (-t185 * t274 + (-t184 * t190 - t272) * t168 * qJD(4)) * pkin(3) + ((-pkin(3) * t272 - t161 * t184) * t168 + t221) * qJD(5) + t309, t68 * t218 + (-t161 * t167 - t5 + (-t249 * t293 + t32) * t168) * t184 + (-t167 * t293 + (-pkin(3) * t253 - qJD(5) * t161 + t33) * t168 + t236) * t189 + t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, t40, t26, -t222, t176, -t178 * t217 + t195, t178 * t233 + t197, t290, t9, t7, t8, t167, -(-t184 * t30 - t284) * t168 + (-t168 * t252 + t54 * t94 + t268) * pkin(4) + t307, (t30 * t168 + t236) * t189 + (-t168 * t251 + t218 * t94 - t274) * pkin(4) + t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t9, t7, t8, t167, -t168 * t221 + t307, (-t6 + (-qJD(5) + t168) * t28) * t189 + t306;];
tau_reg = t4;
