% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP1
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
% tau_reg [6x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:32
% EndTime: 2019-03-08 19:58:39
% DurationCPUTime: 2.89s
% Computational Cost: add. (2763->358), mult. (6454->496), div. (0->0), fcn. (5276->12), ass. (0->195)
t146 = sin(qJ(4));
t149 = cos(qJ(4));
t179 = pkin(4) * t146 - pkin(9) * t149;
t102 = t179 * qJD(4);
t141 = cos(pkin(11));
t147 = sin(qJ(2));
t140 = sin(pkin(6));
t216 = qJD(1) * t140;
t194 = t147 * t216;
t108 = t141 * t194;
t138 = sin(pkin(11));
t150 = cos(qJ(2));
t193 = t150 * t216;
t68 = t138 * t193 + t108;
t271 = t102 - t68;
t215 = qJD(2) * t140;
t190 = qJD(1) * t215;
t201 = qJDD(1) * t140;
t270 = t147 * t201 + t150 * t190;
t206 = qJD(5) * t146;
t269 = qJD(2) * t206 - qJDD(4);
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t200 = t146 * qJDD(2);
t213 = qJD(2) * t149;
t55 = ((qJD(5) + t213) * qJD(4) + t200) * t145 + t269 * t148;
t130 = pkin(5) * t148 + pkin(4);
t268 = t130 * t149 + pkin(3);
t252 = pkin(5) * t145;
t267 = pkin(8) + t252;
t210 = qJD(4) * t146;
t222 = t145 * t149;
t127 = pkin(2) * t138 + pkin(8);
t232 = t127 * t145;
t107 = t138 * t194;
t71 = t141 * t193 - t107;
t266 = -t148 * t271 - t210 * t232 - t71 * t222;
t205 = qJD(5) * t148;
t220 = t148 * t149;
t173 = -pkin(4) * t149 - pkin(9) * t146 - pkin(3);
t253 = pkin(2) * t141;
t88 = t173 - t253;
t265 = t145 * t271 + t88 * t205 - t71 * t220;
t143 = cos(pkin(6));
t123 = qJD(1) * t143 + qJD(3);
t104 = qJD(2) * pkin(2) + t193;
t64 = t138 * t104 + t108;
t62 = qJD(2) * pkin(8) + t64;
t264 = t123 * t149 - t146 * t62;
t174 = t138 * t150 + t141 * t147;
t75 = t174 * t140;
t60 = t143 * t146 + t149 * t75;
t219 = t150 * t141;
t228 = t140 * t147;
t74 = t138 * t228 - t140 * t219;
t24 = -t145 * t60 + t148 * t74;
t142 = cos(pkin(10));
t229 = t140 * t146;
t139 = sin(pkin(10));
t76 = t174 * t143;
t92 = t138 * t147 - t219;
t47 = t139 * t92 - t142 * t76;
t30 = -t142 * t229 - t149 * t47;
t50 = t139 * t76 + t142 * t92;
t32 = t139 * t229 - t149 * t50;
t162 = t92 * t143;
t48 = -t139 * t174 - t142 * t162;
t51 = t139 * t162 - t142 * t174;
t263 = -g(1) * (-t145 * t32 - t148 * t51) - g(2) * (-t145 * t30 - t148 * t48) - g(3) * t24;
t223 = t143 * t150;
t262 = -t139 * t223 - t142 * t147;
t261 = pkin(5) * t55 + qJDD(6);
t124 = -qJD(5) + t213;
t203 = t148 * qJD(4);
t191 = t149 * t203;
t221 = t146 * t148;
t134 = t149 * qJDD(2);
t202 = qJD(2) * qJD(4);
t91 = t146 * t202 + qJDD(5) - t134;
t260 = -t124 * (-t145 * t206 + t191) + t91 * t221;
t164 = g(1) * t51 + g(2) * t48 - g(3) * t74;
t42 = t146 * t123 + t149 * t62;
t39 = qJD(4) * pkin(9) + t42;
t258 = (t124 * t127 + t39) * qJD(5) - t164;
t211 = qJD(4) * t145;
t214 = qJD(2) * t146;
t96 = t148 * t214 + t211;
t257 = t96 ^ 2;
t63 = t104 * t141 - t107;
t46 = t173 * qJD(2) - t63;
t12 = -t145 * t39 + t148 * t46;
t10 = -qJ(6) * t96 + t12;
t7 = -pkin(5) * t124 + t10;
t254 = t10 - t7;
t144 = -qJ(6) - pkin(9);
t100 = t127 * t220;
t171 = pkin(5) * t146 - qJ(6) * t220;
t204 = qJD(6) * t148;
t239 = qJ(6) * t146;
t251 = -t146 * t204 + t171 * qJD(4) + (-t100 + (-t88 + t239) * t145) * qJD(5) - t266;
t212 = qJD(4) * t127;
t250 = (-qJ(6) * qJD(5) - t212) * t221 + (-qJD(6) * t146 + (-qJ(6) * qJD(4) - qJD(5) * t127) * t149) * t145 + t265;
t101 = t179 * qJD(2);
t249 = t145 * t101 + t148 * t264;
t94 = t145 * t214 - t203;
t248 = -t94 * t191 - t55 * t221;
t184 = qJD(5) * t144;
t246 = t204 - t249 + (qJ(6) * t213 + t184) * t145;
t85 = t148 * t101;
t245 = -t171 * qJD(2) + t148 * t184 - t85 + (-qJD(6) + t264) * t145;
t243 = t124 * t94;
t242 = t124 * t96;
t241 = t146 * t96;
t240 = t145 * t88 + t100;
t238 = qJD(4) * t94;
t237 = qJD(5) * t94;
t236 = qJDD(4) * pkin(4);
t120 = t143 * qJDD(1) + qJDD(3);
t235 = t120 * t146;
t233 = t124 * t148;
t230 = t139 * t147;
t227 = t140 * t149;
t226 = t140 * t150;
t224 = t143 * t147;
t218 = qJDD(1) - g(3);
t136 = t146 ^ 2;
t217 = -t149 ^ 2 + t136;
t209 = qJD(4) * t149;
t208 = qJD(5) * t124;
t207 = qJD(5) * t145;
t119 = t150 * t201;
t72 = qJDD(2) * pkin(2) - t147 * t190 + t119;
t35 = -t138 * t270 + t141 * t72;
t19 = qJD(2) * t102 + t173 * qJDD(2) - t35;
t36 = t138 * t72 + t141 * t270;
t34 = qJDD(2) * pkin(8) + t36;
t8 = qJDD(4) * pkin(9) + qJD(4) * t264 + t149 * t34 + t235;
t199 = t145 * t19 + t148 * t8 + t46 * t205;
t197 = t96 * t209;
t195 = t142 * t223;
t192 = t124 * t211;
t189 = t149 * t202;
t186 = t127 + t252;
t54 = -qJD(5) * t203 + (-t189 - t200) * t148 + t269 * t145;
t185 = t149 * t54 + t96 * t210;
t181 = t205 * t241;
t13 = t145 * t46 + t148 * t39;
t11 = -qJ(6) * t94 + t13;
t178 = t11 * t148 - t145 * t7;
t177 = -t11 * t145 - t148 * t7;
t25 = t145 * t74 + t148 * t60;
t175 = t143 * t149 - t146 * t75;
t172 = -t149 * t120 + t123 * t210 + t146 * t34 + t62 * t209;
t38 = -qJD(4) * pkin(4) - t264;
t169 = t39 * t207 - t199;
t168 = t124 * t205 - t145 * t91;
t29 = t142 * t227 - t146 * t47;
t31 = -t139 * t227 - t146 * t50;
t167 = g(1) * t31 + g(2) * t29 - g(3) * t175;
t166 = g(1) * t32 + g(2) * t30 + g(3) * t60;
t165 = g(1) * t50 + g(2) * t47 - g(3) * t75;
t163 = -g(3) * t143 + (-g(1) * t139 + g(2) * t142) * t140;
t9 = t172 - t236;
t160 = -pkin(9) * t91 - t124 * t38;
t128 = -pkin(3) - t253;
t61 = -qJD(2) * pkin(3) - t63;
t159 = -qJDD(4) * t127 + (qJD(2) * t128 + t61 + t71) * qJD(4);
t17 = t148 * t19;
t158 = -t13 * qJD(5) - t145 * t8 + t17;
t157 = pkin(9) * t208 + t167 - t9;
t156 = -g(1) * t262 - g(3) * t226;
t155 = t167 - t172;
t151 = qJD(4) ^ 2;
t153 = -qJD(2) * t68 + t127 * t151 + t164 - t35 + (-pkin(3) + t128) * qJDD(2);
t152 = qJD(2) ^ 2;
t114 = t144 * t148;
t113 = t144 * t145;
t112 = pkin(2) * t195;
t111 = qJDD(4) * t149 - t146 * t151;
t110 = qJDD(4) * t146 + t149 * t151;
t90 = t94 ^ 2;
t79 = t148 * t88;
t70 = t92 * t215;
t69 = qJD(2) * t75;
t53 = -t145 * t239 + t240;
t45 = -qJ(6) * t221 + t79 + (-pkin(5) - t232) * t149;
t23 = t175 * qJD(4) - t149 * t70;
t22 = t60 * qJD(4) - t146 * t70;
t21 = pkin(5) * t94 + qJD(6) + t38;
t5 = t24 * qJD(5) + t145 * t69 + t148 * t23;
t4 = -t25 * qJD(5) - t145 * t23 + t148 * t69;
t3 = t9 + t261;
t2 = -qJ(6) * t55 - qJD(6) * t94 - t169;
t1 = pkin(5) * t91 + qJ(6) * t54 - qJD(6) * t96 + t158;
t6 = [t218, 0 (qJDD(2) * t150 - t147 * t152) * t140 (-qJDD(2) * t147 - t150 * t152) * t140, t120 * t143 - t35 * t74 + t36 * t75 - t63 * t69 - t64 * t70 - g(3), 0, 0, 0, 0, 0, -t74 * t134 - qJD(4) * t22 + qJDD(4) * t175 + (-t149 * t69 + t210 * t74) * qJD(2), t74 * t200 - qJD(4) * t23 - qJDD(4) * t60 + (t146 * t69 + t209 * t74) * qJD(2), 0, 0, 0, 0, 0, -t124 * t4 - t175 * t55 + t22 * t94 + t24 * t91, t124 * t5 + t175 * t54 + t22 * t96 - t25 * t91, t24 * t54 - t25 * t55 - t4 * t96 - t5 * t94, t1 * t24 + t11 * t5 - t175 * t3 + t2 * t25 + t21 * t22 + t4 * t7 - g(3); 0, qJDD(2), t119 - g(2) * (t195 - t230) + t156, -g(1) * (t139 * t224 - t142 * t150) - g(2) * (-t139 * t150 - t142 * t224) - t218 * t228, -g(2) * t112 + t63 * t68 - t64 * t71 + (g(2) * t230 + t36 * t138 + t35 * t141 + t156) * pkin(2), qJDD(2) * t136 + 0.2e1 * t146 * t189, 0.2e1 * t146 * t134 - 0.2e1 * t217 * t202, t110, t111, 0, t146 * t159 - t149 * t153, t146 * t153 + t149 * t159, t96 * t191 + (-t148 * t54 - t207 * t96) * t146, -t181 + (-t197 + (t54 + t237) * t146) * t145 + t248, t185 + t260 (t55 + t192) * t149 + (t168 - t238) * t146, -t124 * t210 - t149 * t91, t79 * t91 + t266 * t124 + (t208 * t88 + t165) * t145 + (t94 * t212 - t17 + (qJD(4) * t38 + qJD(5) * t46 - t127 * t91 + t8) * t145 + t258 * t148) * t149 + (t12 * qJD(4) + t127 * t55 + t9 * t145 + t205 * t38 - t71 * t94) * t146, -t240 * t91 + t265 * t124 + t165 * t148 + ((t127 * t96 + t38 * t148) * qJD(4) - t258 * t145 + t199) * t149 + (-t38 * t207 - t127 * t54 + t9 * t148 - t71 * t96 + (-t127 * t233 - t13) * qJD(4)) * t146, t45 * t54 - t53 * t55 - t251 * t96 - t250 * t94 + t177 * t209 + (-qJD(5) * t178 - t1 * t148 - t145 * t2 - t164) * t146, t2 * t53 + t1 * t45 - g(1) * (t262 * pkin(2) - t267 * t50 + t268 * t51) - g(2) * (-pkin(2) * t230 - t267 * t47 + t268 * t48 + t112) - g(3) * (pkin(2) * t226 + t267 * t75 - t268 * t74) + t251 * t7 + t250 * t11 + t21 * t186 * t209 + (t3 * t186 + (pkin(5) * t205 - t71) * t21 + t164 * t144) * t146; 0, 0, 0, 0, t163 + t120, 0, 0, 0, 0, 0, t111, -t110, 0, 0, 0, 0, 0 (-t55 + t192) * t149 + (t168 + t238) * t146, t185 - t260, t181 + (t197 + (-t54 + t237) * t146) * t145 + t248 (qJD(4) * t178 - t3) * t149 + (qJD(4) * t21 + qJD(5) * t177 - t1 * t145 + t148 * t2) * t146 + t163; 0, 0, 0, 0, 0, -t146 * t152 * t149, t217 * t152, t200, t134, qJDD(4), qJD(4) * t42 - t214 * t61 + t155, -t235 + (-qJD(2) * t61 - t34) * t149 + t166, -t145 * t54 - t96 * t233 (-t54 + t243) * t148 + (-t55 + t242) * t145 (t124 * t220 - t241) * qJD(2) - t168, t124 * t207 + t148 * t91 + (-t124 * t222 + t146 * t94) * qJD(2), t124 * t214, -t12 * t214 - pkin(4) * t55 + t85 * t124 - t42 * t94 + (-t124 * t264 + t160) * t145 + t157 * t148, pkin(4) * t54 - t249 * t124 + t13 * t214 - t157 * t145 + t160 * t148 - t42 * t96, t113 * t54 + t114 * t55 - t245 * t96 - t246 * t94 + (t124 * t7 + t2) * t148 + (t11 * t124 - t1) * t145 - t166, -t2 * t114 + t1 * t113 - t3 * t130 - g(1) * (-t130 * t31 - t32 * t144) - g(2) * (-t130 * t29 - t144 * t30) - g(3) * (t130 * t175 - t144 * t60) + t245 * t7 + (-t124 * t252 - t42) * t21 + t246 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t94, -t90 + t257, -t54 - t243, -t242 - t55, t91, -t13 * t124 - t38 * t96 + t158 + t263, -t12 * t124 + t38 * t94 - g(1) * (t145 * t51 - t148 * t32) - g(2) * (t145 * t48 - t148 * t30) + g(3) * t25 + t169, pkin(5) * t54 + t254 * t94, -t254 * t11 + (-t21 * t96 + t1 + t263) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 - t257, t11 * t94 + t7 * t96 - t155 - t236 + t261;];
tau_reg  = t6;
