% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:15
% EndTime: 2019-03-08 20:34:27
% DurationCPUTime: 4.81s
% Computational Cost: add. (5122->397), mult. (12311->550), div. (0->0), fcn. (10655->18), ass. (0->223)
t175 = cos(qJ(6));
t250 = qJD(6) * t175;
t166 = sin(pkin(12));
t173 = sin(qJ(4));
t271 = t166 * t173;
t232 = qJD(2) * t271;
t177 = cos(qJ(4));
t169 = cos(pkin(12));
t255 = qJD(2) * t169;
t234 = t177 * t255;
t119 = -t232 + t234;
t129 = t166 * t177 + t169 * t173;
t120 = t129 * qJD(2);
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t74 = -t176 * t119 + t120 * t172;
t320 = t175 * t74;
t331 = t250 + t320;
t289 = pkin(8) + qJ(3);
t133 = t289 * t166;
t126 = t177 * t133;
t134 = t289 * t169;
t265 = t177 * t169;
t128 = -t265 + t271;
t168 = sin(pkin(6));
t178 = cos(qJ(2));
t266 = t168 * t178;
t190 = t128 * t266;
t330 = -qJD(1) * t190 + t173 * (qJD(3) * t166 + qJD(4) * t134) - qJD(3) * t265 + qJD(4) * t126;
t165 = qJD(4) + qJD(5);
t283 = t165 * t74;
t252 = qJD(5) * t176;
t253 = qJD(5) * t172;
t245 = t169 * qJDD(2);
t246 = t166 * qJDD(2);
t237 = qJD(4) * t234 + t173 * t245 + t177 * t246;
t88 = -qJD(4) * t232 + t237;
t122 = t129 * qJD(4);
t209 = t173 * t246 - t177 * t245;
t89 = qJD(2) * t122 + t209;
t34 = t119 * t252 - t120 * t253 - t172 * t89 + t176 * t88;
t329 = t34 + t283;
t191 = t129 * t266;
t259 = -t173 * t133 + t177 * t134;
t328 = qJD(1) * t191 - t129 * qJD(3) - t259 * qJD(4);
t317 = -qJD(6) - t74;
t327 = t317 + qJD(6);
t171 = sin(qJ(6));
t161 = qJDD(4) + qJDD(5);
t204 = t119 * t172 + t176 * t120;
t251 = qJD(6) * t171;
t22 = t171 * t161 + t165 * t250 + t175 * t34 - t204 * t251;
t61 = t165 * t171 + t175 * t204;
t23 = qJD(6) * t61 - t175 * t161 + t171 * t34;
t59 = -t175 * t165 + t171 * t204;
t326 = -t171 * t23 + t22 * t175 - t331 * t59;
t20 = t22 * t171;
t325 = t331 * t61 + t20;
t35 = t204 * qJD(5) + t172 * t88 + t176 * t89;
t32 = qJDD(6) + t35;
t29 = t171 * t32;
t290 = t61 * t204;
t63 = t317 * t250;
t324 = -t317 * t320 + t29 - t290 - t63;
t174 = sin(qJ(2));
t257 = qJD(1) * t168;
t233 = t174 * t257;
t132 = qJD(2) * qJ(3) + t233;
t170 = cos(pkin(6));
t256 = qJD(1) * t170;
t143 = t169 * t256;
t93 = t143 + (-pkin(8) * qJD(2) - t132) * t166;
t102 = t169 * t132 + t166 * t256;
t94 = pkin(8) * t255 + t102;
t206 = -t173 * t93 - t177 * t94;
t45 = pkin(9) * t119 - t206;
t280 = t172 * t45;
t307 = -t173 * t94 + t177 * t93;
t44 = -pkin(9) * t120 + t307;
t43 = qJD(4) * pkin(4) + t44;
t24 = t176 * t43 - t280;
t17 = -pkin(5) * t165 - t24;
t323 = t17 * t74;
t275 = cos(pkin(11));
t224 = t275 * t174;
t167 = sin(pkin(11));
t268 = t167 * t178;
t116 = t170 * t224 + t268;
t223 = t275 * t178;
t269 = t167 * t174;
t118 = -t170 * t269 + t223;
t164 = pkin(12) + qJ(4);
t160 = qJ(5) + t164;
t153 = sin(t160);
t154 = cos(t160);
t225 = t168 * t275;
t267 = t168 * t174;
t270 = t167 * t168;
t195 = -g(3) * (-t153 * t267 + t154 * t170) - g(2) * (-t116 * t153 - t154 * t225) - g(1) * (-t118 * t153 + t154 * t270);
t236 = t178 * t257;
t247 = qJDD(2) * qJ(3);
t249 = qJDD(1) * t168;
t100 = t174 * t249 + t247 + (qJD(3) + t236) * qJD(2);
t248 = qJDD(1) * t170;
t140 = t169 * t248;
t64 = t140 + (-pkin(8) * qJDD(2) - t100) * t166;
t82 = t169 * t100 + t166 * t248;
t65 = pkin(8) * t245 + t82;
t226 = -t173 * t65 + t177 * t64;
t13 = qJDD(4) * pkin(4) - pkin(9) * t88 + t206 * qJD(4) + t226;
t207 = t173 * t64 + t177 * t65;
t14 = -pkin(9) * t89 + t307 * qJD(4) + t207;
t276 = t176 * t45;
t25 = t172 * t43 + t276;
t302 = -t25 * qJD(5) + t176 * t13 - t14 * t172;
t3 = -pkin(5) * t161 - t302;
t192 = t195 - t3;
t322 = -pkin(9) * t122 - t330;
t121 = t128 * qJD(4);
t321 = -pkin(9) * t121 - t328;
t319 = t204 * t74;
t318 = t171 * t317;
t284 = t165 * t204;
t316 = -t35 + t284;
t115 = -t170 * t223 + t269;
t117 = t170 * t268 + t224;
t214 = g(1) * t117 + g(2) * t115;
t189 = g(3) * t266 - t214;
t187 = t189 * t154;
t155 = -pkin(3) * t169 - pkin(2);
t106 = pkin(4) * t128 + t155;
t90 = t176 * t128 + t129 * t172;
t91 = -t128 * t172 + t129 * t176;
t39 = pkin(5) * t90 - pkin(10) * t91 + t106;
t314 = t39 * t32 - t187;
t313 = t204 ^ 2 - t74 ^ 2;
t47 = pkin(5) * t204 + pkin(10) * t74;
t104 = t153 * t170 + t154 * t267;
t299 = (qJD(5) * t43 + t14) * t176 + t13 * t172 - t45 * t253;
t84 = t116 * t154 - t153 * t225;
t86 = t118 * t154 + t153 * t270;
t210 = qJD(3) - t236;
t111 = t155 * qJD(2) + t210;
t87 = -pkin(4) * t119 + t111;
t312 = g(1) * t86 + g(2) * t84 + g(3) * t104 + t74 * t87 - t299;
t291 = t59 * t204;
t309 = t317 * t204;
t30 = t175 * t32;
t308 = -t251 * t317 - t30;
t179 = qJD(2) ^ 2;
t194 = (qJDD(2) * t178 - t174 * t179) * t168;
t200 = pkin(4) * t122 - t233;
t18 = pkin(10) * t165 + t25;
t33 = pkin(5) * t74 - pkin(10) * t204 + t87;
t208 = t171 * t18 - t175 * t33;
t306 = t17 * t251 + t204 * t208;
t5 = t171 * t33 + t175 * t18;
t305 = t17 * t250 - t192 * t171 + t5 * t204;
t2 = pkin(10) * t161 + t299;
t213 = g(1) * t118 + g(2) * t116;
t242 = g(3) * t267;
t217 = -t134 * t173 - t126;
t66 = -pkin(9) * t129 + t217;
t67 = -pkin(9) * t128 + t259;
t37 = t172 * t67 - t176 * t66;
t288 = t37 * qJD(5) + t321 * t172 - t322 * t176;
t38 = t172 * t66 + t176 * t67;
t48 = -t90 * qJD(5) - t121 * t176 - t122 * t172;
t304 = -(qJD(6) * t33 + t2) * t90 + t17 * t48 + t3 * t91 - (-qJD(6) * t39 + t288) * t317 - t38 * t32 - t242 - t213;
t303 = -t204 * t87 + t195 + t302;
t254 = qJD(2) * t174;
t231 = qJD(1) * t254;
t244 = t168 * t231 + qJDD(3);
t196 = -t178 * t249 + t244;
t273 = qJDD(2) * pkin(2);
t105 = t196 - t273;
t193 = t214 + t273;
t301 = t168 * (-g(3) * t178 + t231) - t105 + t193;
t205 = (-t132 * t166 + t143) * t166 - t102 * t169;
t298 = t205 * t178 - (-qJD(2) * pkin(2) + t210) * t174;
t293 = t17 * t91;
t287 = t38 * qJD(5) + t322 * t172 + t321 * t176;
t282 = t171 * t61;
t264 = t178 * t179;
t261 = qJDD(1) - g(3);
t258 = t166 ^ 2 + t169 ^ 2;
t240 = t91 * t251;
t239 = t171 * t266;
t238 = t175 * t266;
t235 = t168 * t254;
t220 = t168 * t261;
t156 = pkin(4) * t172 + pkin(10);
t218 = pkin(4) * t120 + qJD(6) * t156 + t47;
t26 = t172 * t44 + t276;
t215 = pkin(4) * t253 - t26;
t49 = t91 * qJD(5) - t121 * t172 + t176 * t122;
t212 = pkin(5) * t49 - pkin(10) * t48 + t200;
t211 = -t317 * t48 + t32 * t91;
t112 = -t166 * t267 + t169 * t170;
t113 = t166 * t170 + t169 * t267;
t71 = t112 * t177 - t113 * t173;
t72 = t112 * t173 + t113 * t177;
t41 = t172 * t72 - t176 * t71;
t42 = t172 * t71 + t176 * t72;
t203 = t318 * t74 - t308;
t199 = -t171 * t42 - t238;
t198 = -t175 * t42 + t239;
t81 = -t100 * t166 + t140;
t186 = -t166 * t81 + t169 * t82 - t213;
t27 = t176 * t44 - t280;
t184 = -t156 * t32 + t323 - (-pkin(4) * t252 + t27) * t317;
t95 = t155 * qJDD(2) + t196;
t52 = pkin(4) * t89 + t95;
t159 = cos(t164);
t158 = sin(t164);
t157 = -pkin(4) * t176 - pkin(5);
t51 = -qJD(2) * t191 - t72 * qJD(4);
t50 = -qJD(2) * t190 + t71 * qJD(4);
t9 = t42 * qJD(5) + t172 * t50 - t176 * t51;
t8 = -t41 * qJD(5) + t172 * t51 + t176 * t50;
t7 = pkin(5) * t35 - pkin(10) * t34 + t52;
t6 = t175 * t7;
t1 = [t261, 0, t194 (-qJDD(2) * t174 - t264) * t168, t169 * t194, -t166 * t194, t258 * t168 * t264 + (-t112 * t166 + t113 * t169) * qJDD(2), t112 * t81 + t113 * t82 - g(3) + (-t298 * qJD(2) - t105 * t178) * t168, 0, 0, 0, 0, 0, qJD(4) * t51 + qJDD(4) * t71 + (-t119 * t254 - t178 * t89) * t168, -qJD(4) * t50 - qJDD(4) * t72 + (t120 * t254 - t178 * t88) * t168, 0, 0, 0, 0, 0, -t161 * t41 - t165 * t9 + (-t178 * t35 + t74 * t254) * t168, -t161 * t42 - t165 * t8 + (-t178 * t34 + t204 * t254) * t168, 0, 0, 0, 0, 0 -(qJD(6) * t198 - t171 * t8 + t175 * t235) * t317 + t199 * t32 + t9 * t59 + t41 * t23 (qJD(6) * t199 + t171 * t235 + t175 * t8) * t317 + t198 * t32 + t9 * t61 + t41 * t22; 0, qJDD(2), t261 * t266 + t214, -t174 * t220 + t213, t301 * t169, -t301 * t166, -t242 + t186 + (t210 * qJD(2) + t247) * t258, -t205 * qJD(3) + (-t105 + t214) * pkin(2) + t186 * qJ(3) + (-g(3) * (pkin(2) * t178 + qJ(3) * t174) + t298 * qJD(1)) * t168, -t120 * t121 + t129 * t88, -t119 * t121 - t120 * t122 - t128 * t88 - t129 * t89, -qJD(4) * t121 + qJDD(4) * t129, -qJD(4) * t122 - qJDD(4) * t128, 0, qJD(4) * t328 + t217 * qJDD(4) + t111 * t122 + t119 * t233 + t95 * t128 + t155 * t89 - t189 * t159, qJD(4) * t330 - t259 * qJDD(4) - t111 * t121 - t120 * t233 + t95 * t129 + t155 * t88 + t189 * t158, t204 * t48 + t34 * t91, -t204 * t49 - t34 * t90 - t35 * t91 - t48 * t74, t161 * t91 + t165 * t48, -t161 * t90 - t165 * t49, 0, t106 * t35 - t161 * t37 - t287 * t165 + t200 * t74 + t49 * t87 + t52 * t90 - t187, t106 * t34 + t189 * t153 - t161 * t38 + t288 * t165 + t200 * t204 + t48 * t87 + t52 * t91, -t61 * t240 + (t22 * t91 + t48 * t61) * t175 (-t175 * t59 - t282) * t48 + (-t20 - t175 * t23 + (t171 * t59 - t175 * t61) * qJD(6)) * t91, t175 * t211 + t22 * t90 + t240 * t317 + t49 * t61, -t171 * t211 - t23 * t90 - t49 * t59 + t63 * t91, -t317 * t49 + t32 * t90, t37 * t23 - t208 * t49 + t6 * t90 + t287 * t59 + (-t212 * t317 + (-t18 * t90 + t317 * t38 + t293) * qJD(6) + t314) * t175 + t304 * t171, t37 * t22 - t5 * t49 + t287 * t61 + (-(-qJD(6) * t18 + t7) * t90 - qJD(6) * t293 - (qJD(6) * t38 - t212) * t317 - t314) * t171 + t304 * t175; 0, 0, 0, 0, -t245, t246, -t258 * t179, t205 * qJD(2) - t178 * t220 - t193 + t244, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t120 + t209 (t119 - t232) * qJD(4) + t237, 0, 0, 0, 0, 0, t35 + t284, t34 - t283, 0, 0, 0, 0, 0, t203 - t291, -t175 * t317 ^ 2 - t29 - t290; 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t119, -t119 ^ 2 + t120 ^ 2 (-t119 - t232) * qJD(4) + t237, -t209, qJDD(4), -t111 * t120 - g(1) * (-t118 * t158 + t159 * t270) - g(2) * (-t116 * t158 - t159 * t225) - g(3) * (-t158 * t267 + t159 * t170) + t226, -t111 * t119 - g(1) * (-t118 * t159 - t158 * t270) - g(2) * (-t116 * t159 + t158 * t225) - g(3) * (-t158 * t170 - t159 * t267) - t207, t319, t313, t329, t316, t161, t165 * t26 + (-t120 * t74 + t161 * t176 - t165 * t253) * pkin(4) + t303, t165 * t27 + (-t120 * t204 - t161 * t172 - t165 * t252) * pkin(4) + t312, t325, t282 * t317 + t326, t324, t203 + t291, t309, t157 * t23 + t215 * t59 + t184 * t171 + (t218 * t317 + t192) * t175 + t306, t157 * t22 + t175 * t184 + t215 * t61 - t218 * t318 + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t319, t313, t329, t316, t161, t165 * t25 + t303, t165 * t24 + t312, t325, t318 * t61 + t326, t324, -t317 * t318 + t291 + t30, t309, -pkin(5) * t23 - t25 * t59 + (-pkin(10) * t32 - t24 * t317 + t323) * t171 + (-(-pkin(10) * qJD(6) - t47) * t317 + t192) * t175 + t306, -pkin(5) * t22 - (t171 * t47 + t175 * t24) * t317 - t25 * t61 + t17 * t320 + t308 * pkin(10) + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t59, -t59 ^ 2 + t61 ^ 2, -t317 * t59 + t22, -t317 * t61 - t23, t32, -t171 * t2 + t6 - t17 * t61 - g(1) * (t117 * t175 - t171 * t86) - g(2) * (t115 * t175 - t171 * t84) - g(3) * (-t104 * t171 - t238) - t327 * t5, -t175 * t2 - t171 * t7 + t17 * t59 - g(1) * (-t117 * t171 - t175 * t86) - g(2) * (-t115 * t171 - t175 * t84) - g(3) * (-t104 * t175 + t239) + t327 * t208;];
tau_reg  = t1;
