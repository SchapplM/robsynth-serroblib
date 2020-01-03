% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:05
% EndTime: 2019-12-31 19:08:15
% DurationCPUTime: 5.57s
% Computational Cost: add. (8580->453), mult. (21235->582), div. (0->0), fcn. (16536->14), ass. (0->235)
t178 = cos(qJ(5));
t250 = qJD(5) * t178;
t172 = cos(pkin(9));
t171 = sin(pkin(9));
t176 = sin(qJ(3));
t261 = t176 * t171;
t294 = cos(qJ(3));
t200 = -t294 * t172 + t261;
t119 = t200 * qJD(1);
t175 = sin(qJ(4));
t293 = cos(qJ(4));
t127 = t294 * t171 + t176 * t172;
t300 = t127 * qJD(1);
t84 = -t293 * t119 - t175 * t300;
t316 = t178 * t84;
t326 = t250 - t316;
t174 = sin(qJ(5));
t228 = qJD(4) * t293;
t249 = t175 * qJD(4);
t231 = qJD(1) * t261;
t224 = qJDD(1) * t294;
t229 = qJD(3) * t294;
t246 = t172 * qJDD(1);
t233 = t172 * qJD(1) * t229 + t171 * t224 + t176 * t246;
t89 = qJD(3) * t231 - t233;
t247 = t171 * qJDD(1);
t207 = -t172 * t224 + t176 * t247;
t310 = t300 * qJD(3);
t90 = t207 + t310;
t197 = -t119 * t228 - t175 * t90 - t249 * t300 - t293 * t89;
t199 = t175 * t119 - t293 * t300;
t245 = qJD(3) + qJD(4);
t217 = t178 * t245;
t244 = qJDD(3) + qJDD(4);
t251 = qJD(5) * t174;
t27 = -qJD(5) * t217 - t174 * t244 - t178 * t197 - t199 * t251;
t72 = t174 * t245 - t178 * t199;
t272 = qJD(5) * t72;
t28 = t174 * t197 - t178 * t244 + t272;
t70 = -t174 * t199 - t217;
t281 = -t174 * t28 - t70 * t250;
t314 = qJD(5) - t84;
t323 = t174 * t314;
t304 = t72 * t323;
t325 = -t178 * t27 + t316 * t70 + t281 - t304;
t282 = pkin(6) + qJ(2);
t142 = t282 * t171;
t128 = qJD(1) * t142;
t143 = t282 * t172;
t129 = qJD(1) * t143;
t94 = -t176 * t128 + t294 * t129;
t69 = -t119 * pkin(7) + t94;
t238 = t293 * t69;
t268 = t129 * t176;
t93 = -t294 * t128 - t268;
t68 = -pkin(7) * t300 + t93;
t67 = qJD(3) * pkin(3) + t68;
t46 = t175 * t67 + t238;
t44 = t245 * pkin(8) + t46;
t159 = t172 * pkin(2) + pkin(1);
t134 = -t159 * qJD(1) + qJD(2);
t95 = pkin(3) * t119 + t134;
t47 = -pkin(4) * t84 + pkin(8) * t199 + t95;
t14 = t174 * t47 + t178 * t44;
t289 = t14 * t314;
t225 = t175 * t89 - t293 * t90;
t311 = t199 * qJD(4);
t42 = -t225 - t311;
t133 = -t159 * qJDD(1) + qJDD(2);
t295 = pkin(3) * t90;
t73 = t133 + t295;
t12 = pkin(4) * t42 - pkin(8) * t197 + t73;
t11 = t178 * t12;
t248 = qJD(1) * qJD(2);
t298 = t282 * qJDD(1) + t248;
t103 = t298 * t171;
t104 = t298 * t172;
t223 = -t294 * t103 - t176 * t104;
t57 = -t94 * qJD(3) + t223;
t37 = qJDD(3) * pkin(3) + t89 * pkin(7) + t57;
t240 = t176 * t103 - t294 * t104 + t128 * t229;
t252 = qJD(3) * t176;
t56 = -t129 * t252 - t240;
t40 = -pkin(7) * t90 + t56;
t9 = t175 * t37 + t67 * t228 - t69 * t249 + t293 * t40;
t7 = t244 * pkin(8) + t9;
t3 = -qJD(5) * t14 - t174 * t7 + t11;
t324 = t3 + t289;
t13 = -t174 * t44 + t178 * t47;
t290 = t13 * t314;
t270 = qJDD(1) * pkin(1);
t177 = sin(qJ(1));
t179 = cos(qJ(1));
t302 = -g(1) * t177 + g(2) * t179;
t204 = -qJDD(2) + t270 - t302;
t305 = t245 * t84;
t322 = t197 - t305;
t284 = t84 ^ 2;
t285 = t199 ^ 2;
t321 = -t284 + t285;
t24 = t27 * t174;
t319 = t326 * t72 - t24;
t286 = t72 * t199;
t39 = qJDD(5) + t42;
t32 = t174 * t39;
t318 = t314 * t326 + t286 + t32;
t275 = t175 * t69;
t45 = t293 * t67 - t275;
t43 = -t245 * pkin(4) - t45;
t317 = t43 * t84;
t170 = pkin(9) + qJ(3);
t165 = qJ(4) + t170;
t158 = cos(t165);
t291 = g(3) * t158;
t226 = t175 * t40 - t293 * t37;
t10 = -t46 * qJD(4) - t226;
t8 = -t244 * pkin(4) - t10;
t232 = -t8 - t291;
t288 = t70 * t199;
t315 = t314 * t199;
t283 = t84 * t199;
t157 = sin(t165);
t266 = t157 * t179;
t267 = t157 * t177;
t313 = g(1) * t266 + g(2) * t267;
t258 = t199 * qJD(3);
t312 = -t258 + t225;
t152 = g(3) * t157;
t264 = t158 * t179;
t265 = t158 * t177;
t234 = -g(1) * t264 - g(2) * t265 - t152;
t309 = -t95 * t84 - t234 - t9;
t34 = t43 * t251;
t308 = t13 * t199 + t313 * t178 + t34;
t35 = t43 * t250;
t307 = -t14 * t199 - t232 * t174 + t35;
t306 = t95 * t199 - t226 - t291 + t313;
t55 = -pkin(4) * t199 - pkin(8) * t84;
t97 = -t176 * t142 + t294 * t143;
t303 = t158 * pkin(4) + t157 * pkin(8);
t301 = -t13 * t174 + t14 * t178;
t299 = qJ(2) * qJDD(1);
t297 = t300 ^ 2;
t296 = qJD(3) ^ 2;
t2 = qJD(5) * t13 + t12 * t174 + t178 * t7;
t1 = t2 * t178;
t287 = t72 * t70;
t276 = t174 * t70;
t26 = t28 * t178;
t271 = qJD(5) * t314;
t269 = t300 * t119;
t263 = t174 * t177;
t262 = t174 * t179;
t260 = t177 * t178;
t259 = t178 * t179;
t230 = qJD(2) * t294;
t256 = -t142 * t229 + t172 * t230;
t168 = t171 ^ 2;
t169 = t172 ^ 2;
t254 = t168 + t169;
t253 = qJD(2) * t176;
t241 = pkin(8) * t271;
t160 = pkin(3) * t175 + pkin(8);
t237 = t160 * t271;
t195 = t175 * t200;
t92 = t293 * t127 - t195;
t236 = t92 * t251;
t235 = t92 * t250;
t219 = t254 * qJD(1) ^ 2;
t96 = -t294 * t142 - t176 * t143;
t218 = t1 + t234;
t216 = pkin(3) * t228;
t215 = 0.2e1 * t254;
t48 = t175 * t68 + t238;
t214 = pkin(3) * t249 - t48;
t213 = -pkin(8) * t39 - t317;
t189 = t293 * t200;
t193 = qJD(3) * t127;
t58 = t127 * t249 + t175 * t193 + t245 * t189;
t212 = -t43 * t58 + t8 * t92;
t210 = g(1) * t179 + g(2) * t177;
t208 = -t314 * t58 + t39 * t92;
t206 = t13 * t178 + t14 * t174;
t78 = -pkin(7) * t127 + t96;
t79 = -t200 * pkin(7) + t97;
t52 = t175 * t78 + t293 * t79;
t101 = t200 * pkin(3) - t159;
t91 = t127 * t175 + t189;
t53 = t91 * pkin(4) - t92 * pkin(8) + t101;
t22 = t174 * t53 + t178 * t52;
t21 = -t174 * t52 + t178 * t53;
t33 = t178 * t39;
t205 = -t251 * t314 + t323 * t84 + t33;
t203 = -qJD(5) * t47 + t152 - t7;
t202 = t210 * t157;
t163 = sin(t170);
t201 = t210 * t163;
t194 = qJD(3) * t200;
t191 = t204 + t270;
t190 = pkin(3) * t193;
t164 = cos(t170);
t188 = -g(3) * t164 + t201;
t187 = t133 + t302;
t186 = -t160 * t39 - t216 * t314 - t317;
t185 = t215 * t248 - t210;
t184 = -t206 * qJD(5) - t3 * t174 + t1;
t183 = -g(1) * (-pkin(4) * t266 + pkin(8) * t264) - g(2) * (-pkin(4) * t267 + pkin(8) * t265);
t182 = pkin(7) * t194 + t142 * t252 - t143 * t229 - t171 * t230 - t172 * t253;
t167 = -pkin(7) - t282;
t161 = -t293 * pkin(3) - pkin(4);
t156 = pkin(3) * t164;
t132 = t156 + t159;
t122 = t179 * t132;
t116 = t119 ^ 2;
t113 = t158 * t259 + t263;
t112 = -t158 * t262 + t260;
t111 = -t158 * t260 + t262;
t110 = t158 * t263 + t259;
t75 = -t127 * qJD(2) - t97 * qJD(3);
t74 = (-qJD(2) * t171 - qJD(3) * t143) * t176 + t256;
t61 = -pkin(7) * t193 - t143 * t252 - t171 * t253 + t256;
t59 = -qJD(4) * t195 + t127 * t228 - t175 * t194 + t293 * t193;
t51 = t175 * t79 - t293 * t78;
t50 = pkin(3) * t300 + t55;
t49 = t293 * t68 - t275;
t29 = t59 * pkin(4) + t58 * pkin(8) + t190;
t20 = t174 * t55 + t178 * t45;
t19 = -t174 * t45 + t178 * t55;
t18 = t174 * t50 + t178 * t49;
t17 = -t174 * t49 + t178 * t50;
t16 = t52 * qJD(4) + t175 * t61 - t293 * t182;
t15 = t175 * t182 + t78 * t228 - t79 * t249 + t293 * t61;
t5 = -t22 * qJD(5) - t15 * t174 + t178 * t29;
t4 = t21 * qJD(5) + t15 * t178 + t174 * t29;
t6 = [0, 0, 0, 0, 0, qJDD(1), -t302, t210, 0, 0, t168 * qJDD(1), 0.2e1 * t171 * t246, 0, t169 * qJDD(1), 0, 0, t191 * t172, -t191 * t171, t215 * t299 + t185, t204 * pkin(1) + (t254 * t299 + t185) * qJ(2), -t89 * t127 - t194 * t300, t119 * t194 - t127 * t90 - t193 * t300 + t89 * t200, t127 * qJDD(3) - t200 * t296, t119 * t193 + t200 * t90, -qJDD(3) * t200 - t296 * t127, 0, t96 * qJDD(3) + t133 * t200 - t159 * t90 - t302 * t164 + (t127 * t134 + t75) * qJD(3), -t74 * qJD(3) - t97 * qJDD(3) + t133 * t127 - t134 * t194 + t159 * t89 + t302 * t163, -t74 * t119 - t75 * t300 - t57 * t127 - t193 * t94 + t96 * t89 - t97 * t90 - t210 + (qJD(3) * t93 - t56) * t200, t56 * t97 + t94 * t74 + t57 * t96 + t93 * t75 - t133 * t159 - g(1) * (-t159 * t177 + t179 * t282) - g(2) * (t159 * t179 + t177 * t282), t197 * t92 + t199 * t58, -t197 * t91 + t199 * t59 - t42 * t92 - t58 * t84, t244 * t92 - t245 * t58, t42 * t91 - t59 * t84, -t244 * t91 - t245 * t59, 0, g(1) * t265 - g(2) * t264 + t101 * t42 - t16 * t245 - t190 * t84 - t244 * t51 + t95 * t59 + t73 * t91, -g(1) * t267 + g(2) * t266 + t101 * t197 - t15 * t245 - t190 * t199 - t244 * t52 - t95 * t58 + t73 * t92, -t10 * t92 + t15 * t84 - t16 * t199 + t197 * t51 - t42 * t52 + t45 * t58 - t46 * t59 - t9 * t91 - t210, t9 * t52 + t46 * t15 - t10 * t51 - t45 * t16 + t73 * t101 + t95 * t190 - g(1) * (-t132 * t177 - t167 * t179) - g(2) * (-t167 * t177 + t122), -t72 * t236 + (-t27 * t92 - t58 * t72) * t178, (t174 * t72 + t178 * t70) * t58 + (t24 - t26 + (-t178 * t72 + t276) * qJD(5)) * t92, t208 * t178 - t236 * t314 - t27 * t91 + t59 * t72, t70 * t235 + (t28 * t92 - t58 * t70) * t174, -t174 * t208 - t235 * t314 - t28 * t91 - t59 * t70, t314 * t59 + t39 * t91, -g(1) * t111 - g(2) * t113 + t13 * t59 + t16 * t70 + t174 * t212 + t21 * t39 + t28 * t51 + t3 * t91 + t314 * t5 + t35 * t92, -g(1) * t110 - g(2) * t112 - t14 * t59 + t16 * t72 + t178 * t212 - t2 * t91 - t22 * t39 - t27 * t51 - t314 * t4 - t34 * t92, t21 * t27 - t22 * t28 - t4 * t70 - t5 * t72 + t206 * t58 - t302 * t157 + (-t301 * qJD(5) - t174 * t2 - t178 * t3) * t92, -g(2) * t122 + t13 * t5 + t14 * t4 + t43 * t16 + t2 * t22 + t3 * t21 + t8 * t51 + (g(1) * t167 - g(2) * t303) * t179 + (-g(1) * (-t132 - t303) + g(2) * t167) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t247, -t219, -qJ(2) * t219 - t204, 0, 0, 0, 0, 0, 0, t207 + 0.2e1 * t310, (-t119 - t231) * qJD(3) + t233, -t116 - t297, t119 * t94 + t300 * t93 + t187, 0, 0, 0, 0, 0, 0, -t225 - t258 - 0.2e1 * t311, t197 + t305, -t284 - t285, -t199 * t45 - t46 * t84 + t187 + t295, 0, 0, 0, 0, 0, 0, t205 + t288, -t178 * t314 ^ 2 + t286 - t32, (t70 * t84 + t27) * t178 + t304 + t281, t43 * t199 + t324 * t178 + (t2 - t290) * t174 + t302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, -t116 + t297, (t119 - t231) * qJD(3) + t233, -t269, -t207, qJDD(3), -t134 * t300 + t188 + t223, g(3) * t163 + t119 * t134 + t210 * t164 + (t93 + t268) * qJD(3) + t240, 0, 0, t283, t321, t322, -t283, t312, t244, t48 * qJD(3) + (t48 - t46) * qJD(4) + (t293 * t244 - t245 * t249 + t300 * t84) * pkin(3) + t306, t49 * t245 + (-t175 * t244 + t199 * t300 - t228 * t245) * pkin(3) + t309, t45 * t84 - t46 * t199 + t48 * t199 - t49 * t84 + (-t293 * t197 - t175 * t42 + (-t175 * t199 + t293 * t84) * qJD(4)) * pkin(3), t45 * t48 - t46 * t49 + (t293 * t10 - t300 * t95 + t175 * t9 + (-t175 * t45 + t293 * t46) * qJD(4) + t188) * pkin(3), t319, t325, t318, t323 * t70 - t26, t205 - t288, t315, t161 * t28 - t17 * t314 + t214 * t70 + (t232 - t237) * t178 + t186 * t174 + t308, -t161 * t27 + t18 * t314 + t214 * t72 + t186 * t178 + (-t202 + t237) * t174 + t307, t17 * t72 + t18 * t70 + (-t70 * t216 + t13 * t84 - t160 * t28 + (t160 * t72 - t13) * qJD(5)) * t178 + (t72 * t216 + t14 * t84 - t160 * t27 - t3 + (t160 * t70 - t14) * qJD(5)) * t174 + t218, t8 * t161 - t14 * t18 - t13 * t17 - t43 * t48 - g(3) * (t156 + t303) + (t201 + (t175 * t43 + t301 * t293) * qJD(4)) * pkin(3) + t184 * t160 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, t321, t322, -t283, t312, t244, t46 * qJD(3) + t306, t245 * t45 + t309, 0, 0, t319, t325, t318, t276 * t314 - t26, -t314 * t323 - t288 + t33, t315, -pkin(4) * t28 - t19 * t314 - t46 * t70 + t213 * t174 + (t232 - t241) * t178 + t308, pkin(4) * t27 + t20 * t314 - t46 * t72 + t213 * t178 + (-t202 + t241) * t174 + t307, t19 * t72 + t20 * t70 + (-t290 + (-t28 + t272) * pkin(8)) * t178 + ((qJD(5) * t70 - t27) * pkin(8) - t324) * t174 + t218, -t8 * pkin(4) + pkin(8) * t184 - g(3) * t303 - t13 * t19 - t14 * t20 - t43 * t46 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t70 ^ 2 + t72 ^ 2, t314 * t70 - t27, -t287, t314 * t72 - t28, t39, -g(1) * t112 + g(2) * t110 + t174 * t203 - t250 * t44 - t43 * t72 + t11 + t289, g(1) * t113 - g(2) * t111 + t290 + t43 * t70 + (qJD(5) * t44 - t12) * t174 + t203 * t178, 0, 0;];
tau_reg = t6;
