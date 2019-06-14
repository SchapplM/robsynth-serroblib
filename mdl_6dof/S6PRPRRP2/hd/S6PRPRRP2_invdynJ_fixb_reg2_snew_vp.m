% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:38:17
% EndTime: 2019-05-04 23:38:24
% DurationCPUTime: 3.09s
% Computational Cost: add. (7201->301), mult. (13463->404), div. (0->0), fcn. (9710->12), ass. (0->205)
t179 = cos(qJ(5));
t176 = sin(qJ(5));
t177 = sin(qJ(4));
t220 = qJD(2) * qJD(4);
t158 = t177 * t220;
t180 = cos(qJ(4));
t218 = t180 * qJDD(2);
t139 = -t158 + t218;
t132 = -qJDD(5) + t139;
t222 = qJD(2) * t177;
t133 = -t179 * qJD(4) + t176 * t222;
t135 = qJD(4) * t176 + t179 * t222;
t229 = t135 * t133;
t189 = t132 - t229;
t240 = t176 * t189;
t131 = t135 ^ 2;
t155 = qJD(2) * t180 - qJD(5);
t247 = t155 ^ 2;
t256 = -t131 - t247;
t60 = -t179 * t256 - t240;
t302 = pkin(3) * t60;
t301 = pkin(4) * t60;
t300 = pkin(9) * t60;
t235 = t179 * t189;
t66 = -t176 * t256 + t235;
t299 = pkin(9) * t66;
t168 = sin(pkin(11));
t298 = t168 * t60;
t171 = cos(pkin(11));
t297 = t171 * t60;
t296 = t177 * t66;
t295 = t180 * t66;
t89 = t132 + t229;
t234 = t179 * t89;
t248 = t133 ^ 2;
t251 = -t247 - t248;
t261 = t176 * t251 - t234;
t121 = t135 * t155;
t214 = t180 * t220;
t219 = t177 * qJDD(2);
t138 = t214 + t219;
t209 = t179 * qJDD(4) - t176 * t138;
t193 = qJD(5) * t135 - t209;
t257 = -t121 + t193;
t239 = t176 * t89;
t260 = t179 * t251 + t239;
t274 = t177 * t257 + t180 * t260;
t290 = t168 * t274 - t171 * t261;
t294 = pkin(2) * t290 - pkin(3) * t261 + pkin(8) * t274;
t170 = sin(pkin(6));
t173 = cos(pkin(6));
t178 = sin(qJ(2));
t181 = cos(qJ(2));
t275 = t177 * t260 - t180 * t257;
t293 = t173 * t275 + (t181 * t290 + t178 * (t168 * t261 + t171 * t274)) * t170;
t115 = t248 - t247;
t72 = t121 + t193;
t289 = t177 * (t115 * t179 + t240) + t180 * t72;
t197 = -t176 * qJDD(4) - t179 * t138;
t188 = -qJD(5) * t133 - t197;
t230 = t133 * t155;
t252 = t230 + t188;
t242 = t176 * t252;
t255 = t131 - t248;
t288 = t177 * (t179 * t257 + t242) + t180 * t255;
t285 = pkin(4) * t261;
t284 = pkin(9) * t260;
t283 = pkin(9) * t261;
t282 = qJ(6) * t252;
t116 = -t131 + t247;
t279 = t179 * t116 - t239;
t277 = t115 * t176 - t235;
t253 = -t230 + t188;
t273 = t177 * (-t116 * t176 - t234) - t180 * t253;
t254 = t131 + t248;
t272 = pkin(4) * t254;
t169 = sin(pkin(10));
t172 = cos(pkin(10));
t201 = g(1) * t169 - g(2) * t172;
t165 = -g(3) + qJDD(1);
t208 = t173 * t165 + qJDD(3);
t187 = -t170 * t201 + t208;
t108 = t180 * t187;
t182 = qJD(2) ^ 2;
t204 = -pkin(4) * t180 - pkin(9) * t177;
t145 = -g(1) * t172 - g(2) * t169;
t258 = t173 * t201;
t262 = t170 * t165 + t258;
t93 = -t145 * t178 + t181 * t262;
t191 = qJDD(2) * pkin(2) + t93;
t94 = t181 * t145 + t178 * t262;
t85 = -t182 * pkin(2) + t94;
t59 = t168 * t191 + t171 * t85;
t57 = -pkin(3) * t182 + qJDD(2) * pkin(8) + t59;
t210 = t182 * t204 + t57;
t246 = qJD(4) ^ 2;
t36 = -qJDD(4) * pkin(4) - t246 * pkin(9) + t210 * t177 - t108;
t271 = pkin(5) * t193 - t282 + t36;
t269 = t177 * t254;
t265 = t180 * t254;
t259 = -t176 * t257 + t179 * t252;
t101 = pkin(5) * t133 - qJ(6) * t135;
t185 = t177 * t187;
t37 = -t246 * pkin(4) + qJDD(4) * pkin(9) + t210 * t180 + t185;
t199 = -t139 + t158;
t200 = t138 + t214;
t212 = t168 * t85 - t171 * t191;
t56 = -qJDD(2) * pkin(3) - t182 * pkin(8) + t212;
t44 = t199 * pkin(4) - t200 * pkin(9) + t56;
t20 = t176 * t44 + t179 * t37;
t211 = -t132 * qJ(6) - t133 * t101 + t20;
t249 = -(t256 + t247) * pkin(5) - qJ(6) * t189 + t211;
t245 = pkin(5) * t179;
t244 = t176 * t36;
t241 = t176 * t253;
t238 = t179 * t36;
t236 = t179 * t253;
t232 = qJ(6) * t179;
t228 = t155 * t176;
t227 = t155 * t179;
t154 = t177 * t182 * t180;
t146 = qJDD(4) + t154;
t225 = t177 * t146;
t147 = qJDD(4) - t154;
t224 = t180 * t147;
t221 = qJD(6) * t155;
t217 = t133 * t227;
t216 = t180 * t229;
t213 = -qJ(6) * t176 - pkin(4);
t19 = t176 * t37 - t179 * t44;
t10 = t176 * t19 + t179 * t20;
t50 = t177 * t57 - t108;
t51 = t180 * t57 + t185;
t25 = t177 * t50 + t180 * t51;
t114 = t135 * t228;
t207 = t177 * (t179 * t188 + t114) - t216;
t206 = -t133 * t228 - t179 * t193;
t148 = -0.2e1 * t221;
t205 = t148 + t211;
t13 = -pkin(5) * t247 + t205;
t15 = t132 * pkin(5) - qJ(6) * t247 + t101 * t135 + qJDD(6) + t19;
t203 = -pkin(5) * t15 + qJ(6) * t13;
t202 = -pkin(5) * t253 - qJ(6) * t72;
t9 = t176 * t20 - t179 * t19;
t198 = -pkin(3) + t204;
t192 = (t133 * t176 + t135 * t179) * t155;
t190 = t177 * (-t114 + t217) + t180 * t132;
t186 = t177 * (t176 * t193 - t217) + t216;
t184 = 0.2e1 * qJD(6) * t135 - t271;
t183 = -pkin(5) * t89 + qJ(6) * t251 - t15;
t164 = t180 ^ 2;
t163 = t177 ^ 2;
t162 = t164 * t182;
t160 = t163 * t182;
t152 = -t162 - t246;
t151 = -t160 - t246;
t144 = t160 + t162;
t143 = (t163 + t164) * qJDD(2);
t142 = -qJDD(2) * t168 - t171 * t182;
t141 = qJDD(2) * t171 - t168 * t182;
t140 = -0.2e1 * t158 + t218;
t137 = 0.2e1 * t214 + t219;
t113 = -t151 * t177 - t224;
t112 = t152 * t180 - t225;
t111 = -t147 * t177 + t151 * t180;
t110 = t146 * t180 + t152 * t177;
t102 = t143 * t168 + t144 * t171;
t83 = t113 * t168 - t137 * t171;
t82 = t112 * t168 + t140 * t171;
t78 = (qJD(5) - t155) * t133 + t197;
t73 = (-qJD(5) - t155) * t135 + t209;
t69 = -t135 * t227 + t176 * t188;
t55 = t179 * t73 + t241;
t54 = -t179 * t72 + t241;
t53 = t176 * t73 - t236;
t52 = -t176 * t72 - t236;
t47 = -t177 * t78 + t295;
t45 = t180 * t78 + t296;
t40 = -t177 * t252 - t295;
t38 = t180 * t252 - t296;
t34 = t180 * t55 - t269;
t33 = t180 * t54 - t269;
t32 = t177 * t55 + t265;
t31 = t177 * t54 + t265;
t30 = t168 * t59 - t171 * t212;
t28 = t168 * t47 + t297;
t26 = t168 * t40 - t297;
t24 = t177 * t51 - t180 * t50;
t23 = (-pkin(5) * t155 - 0.2e1 * qJD(6)) * t135 + t271;
t22 = t168 * t34 - t171 * t53;
t21 = t168 * t33 - t171 * t52;
t17 = (-t257 + t121) * pkin(5) + t184;
t16 = pkin(5) * t121 + t184 + t282;
t14 = t168 * t25 - t171 * t56;
t12 = qJ(6) * t254 + t15;
t11 = (-t247 + t254) * pkin(5) + t205;
t8 = t13 * t179 + t15 * t176;
t7 = t13 * t176 - t15 * t179;
t6 = t10 * t180 + t177 * t36;
t5 = t10 * t177 - t180 * t36;
t4 = t177 * t23 + t180 * t8;
t3 = t177 * t8 - t180 * t23;
t2 = t168 * t6 - t171 * t9;
t1 = t168 * t4 - t171 * t7;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, 0, 0, 0, 0, 0, (qJDD(2) * t181 - t178 * t182) * t170, (-qJDD(2) * t178 - t181 * t182) * t170, 0, t173 ^ 2 * t165 + (t178 * t94 + t181 * t93 - t258) * t170, 0, 0, 0, 0, 0, 0, (t141 * t181 + t142 * t178) * t170, (-t141 * t178 + t142 * t181) * t170, 0, t173 * t208 + (t178 * (t168 * t212 + t171 * t59) + t181 * t30 - t258) * t170, 0, 0, 0, 0, 0, 0, t173 * t110 + (t178 * (t112 * t171 - t140 * t168) + t181 * t82) * t170, t173 * t111 + (t178 * (t113 * t171 + t137 * t168) + t181 * t83) * t170, (t178 * (t143 * t171 - t144 * t168) + t181 * t102) * t170, t173 * t24 + (t178 * (t168 * t56 + t171 * t25) + t181 * t14) * t170, 0, 0, 0, 0, 0, 0, t293, t173 * t45 + (t178 * (t171 * t47 - t298) + t181 * t28) * t170, t173 * t32 + (t178 * (t168 * t53 + t171 * t34) + t181 * t22) * t170, t173 * t5 + (t178 * (t168 * t9 + t171 * t6) + t181 * t2) * t170, 0, 0, 0, 0, 0, 0, t293, t173 * t31 + (t178 * (t168 * t52 + t171 * t33) + t181 * t21) * t170, t173 * t38 + (t178 * (t171 * t40 + t298) + t181 * t26) * t170, t173 * t3 + (t178 * (t168 * t7 + t171 * t4) + t181 * t1) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t93, -t94, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t141 - t212, pkin(2) * t142 - t59, 0, pkin(2) * t30, t200 * t177, t137 * t180 + t140 * t177, t225 + t180 * (-t160 + t246), -t199 * t180, t177 * (t162 - t246) + t224, 0, pkin(2) * t82 + pkin(3) * t140 + pkin(8) * t112 - t180 * t56, pkin(2) * t83 - pkin(3) * t137 + pkin(8) * t113 + t177 * t56, pkin(2) * t102 + pkin(3) * t144 + pkin(8) * t143 + t25, pkin(2) * t14 - pkin(3) * t56 + pkin(8) * t25, t207, -t288, t273, t186, t289, t190, t177 * (t244 - t283) + t180 * (t19 - t285) + t294, t177 * (t238 + t300) + t180 * (t20 + t301) + t302 + pkin(8) * t47 + pkin(2) * t28, pkin(2) * t22 + pkin(8) * t34 - t177 * t9 + t198 * t53, pkin(2) * t2 + pkin(8) * t6 + t198 * t9, t207, t273, t288, t190, -t289, t186, t177 * (-t17 * t176 - t232 * t257 - t283) + t180 * (-t183 - t285) + t294, t177 * (-pkin(9) * t52 - t11 * t176 + t12 * t179) + t180 * (-pkin(4) * t52 - t202) - pkin(3) * t52 + pkin(8) * t33 + pkin(2) * t21, t177 * (-pkin(5) * t242 + t16 * t179 - t300) + t180 * (0.2e1 * t221 - t249 - t301) - t302 + pkin(8) * t40 + pkin(2) * t26, t177 * (-pkin(9) * t7 + (pkin(5) * t176 - t232) * t23) + t180 * (-pkin(4) * t7 - t203) - pkin(3) * t7 + pkin(8) * t4 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, 0, 0, 0, 0, 0, 0, t110, t111, 0, t24, 0, 0, 0, 0, 0, 0, t275, t45, t32, t5, 0, 0, 0, 0, 0, 0, t275, t31, t38, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t160 - t162, t219, t154, t218, qJDD(4), -t50, -t51, 0, 0, t69, t259, t279, t206, t277, t192, -pkin(4) * t257 - t238 + t284, pkin(4) * t78 + t244 + t299, pkin(9) * t55 + t10 + t272, -pkin(4) * t36 + pkin(9) * t10, t69, t279, -t259, t192, -t277, t206, t179 * t17 + t213 * t257 + t284, pkin(9) * t54 + t11 * t179 + t12 * t176 + t272, -t299 + t176 * t16 + (pkin(4) + t245) * t252, pkin(9) * t8 + (t213 - t245) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t255, t253, -t229, -t72, -t132, -t19, -t20, 0, 0, t229, t253, -t255, -t132, t72, -t229, t183, t202, t148 + t249, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t253, t256, t15;];
tauJ_reg  = t18;
