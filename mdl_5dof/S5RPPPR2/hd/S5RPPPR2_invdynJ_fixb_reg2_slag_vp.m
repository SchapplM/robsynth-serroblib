% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:23:05
% EndTime: 2020-01-03 11:23:13
% DurationCPUTime: 2.94s
% Computational Cost: add. (2911->337), mult. (7297->508), div. (0->0), fcn. (5685->10), ass. (0->186)
t148 = sin(pkin(9));
t151 = cos(pkin(9));
t153 = cos(pkin(7));
t202 = qJD(1) * qJD(4);
t149 = sin(pkin(8));
t152 = cos(pkin(8));
t197 = t153 * qJDD(1);
t183 = t152 * t197;
t203 = qJD(1) * qJD(2);
t188 = t153 * t203;
t150 = sin(pkin(7));
t241 = t153 * pkin(2);
t102 = -t150 * qJ(3) - pkin(1) - t241;
t206 = qJD(3) * t150;
t66 = -qJD(1) * t206 + qJDD(1) * t102 + qJDD(2);
t33 = qJ(2) * t183 + t149 * t66 + t152 * t188;
t26 = (-qJ(4) * qJDD(1) - t202) * t153 + t33;
t172 = pkin(3) * t149 - qJ(4) * t152;
t198 = t150 * qJDD(1);
t99 = qJ(2) * t198 + t150 * t203 + qJDD(3);
t37 = (qJDD(1) * t172 - t152 * t202) * t150 + t99;
t12 = t148 * t37 + t151 * t26;
t186 = t149 * t198;
t10 = pkin(6) * t186 + t12;
t184 = t149 * t197;
t32 = -qJ(2) * t184 - t149 * t188 + t152 * t66;
t29 = pkin(3) * t197 + qJDD(4) - t32;
t223 = t153 * t151;
t226 = t150 * t152;
t89 = t148 * t226 + t223;
t75 = qJDD(1) * t89;
t117 = t148 * t197;
t200 = qJDD(1) * t152;
t185 = t151 * t200;
t76 = t150 * t185 - t117;
t14 = t75 * pkin(4) - t76 * pkin(6) + t29;
t154 = sin(qJ(5));
t156 = cos(qJ(5));
t208 = qJD(1) * t153;
t194 = qJ(2) * t208;
t88 = qJD(1) * t102 + qJD(2);
t49 = -t149 * t194 + t152 * t88;
t39 = pkin(3) * t208 + qJD(4) - t49;
t77 = t89 * qJD(1);
t210 = qJD(1) * t150;
t189 = t151 * t210;
t190 = t148 * t208;
t80 = t152 * t189 - t190;
t15 = t77 * pkin(4) - t80 * pkin(6) + t39;
t50 = t149 * t88 + t152 * t194;
t40 = -qJ(4) * t208 + t50;
t116 = qJ(2) * t210 + qJD(3);
t62 = t172 * t210 + t116;
t19 = t148 * t62 + t151 * t40;
t193 = t149 * t210;
t17 = pkin(6) * t193 + t19;
t6 = t154 * t15 + t156 * t17;
t2 = -qJD(5) * t6 - t154 * t10 + t156 * t14;
t249 = qJD(5) + t77;
t251 = t249 * t6 + t2;
t211 = qJD(1) * t149;
t191 = t156 * t211;
t176 = t150 * t191;
t205 = qJD(5) * t154;
t20 = -qJD(5) * t176 - t154 * t186 - t156 * t76 + t205 * t80;
t41 = t154 * t80 - t176;
t250 = t249 * t41 - t20;
t232 = qJDD(1) * pkin(1);
t135 = qJDD(2) - t232;
t155 = sin(qJ(1));
t157 = cos(qJ(1));
t245 = -g(2) * t157 - g(3) * t155;
t248 = t245 - t135;
t224 = t150 * t157;
t218 = t157 * t152;
t221 = t155 * t149;
t98 = t153 * t218 + t221;
t57 = t148 * t224 + t98 * t151;
t219 = t157 * t149;
t220 = t155 * t152;
t97 = t153 * t219 - t220;
t247 = t57 * t154 - t97 * t156;
t246 = t97 * t154 + t57 * t156;
t244 = t249 - qJD(5);
t204 = qJ(2) * qJDD(1);
t243 = g(2) * t155;
t242 = g(3) * t157;
t192 = t154 * t211;
t44 = t150 * t192 + t156 * t80;
t240 = t44 * t41;
t82 = (t148 * t150 + t152 * t223) * qJD(1);
t227 = t149 * t156;
t96 = t151 * t227 - t154 * t152;
t239 = t96 * qJD(5) + t153 * t191 - t154 * t82;
t228 = t149 * t154;
t168 = t151 * t228 + t156 * t152;
t238 = t168 * qJD(5) + t153 * t192 + t156 * t82;
t233 = qJ(2) * t153;
t69 = t149 * t102 + t152 * t233;
t59 = -t153 * qJ(4) + t69;
t73 = (qJ(2) + t172) * t150;
t28 = t148 * t73 + t151 * t59;
t71 = qJDD(5) + t75;
t237 = t154 * t71;
t236 = t156 * t71;
t145 = t150 ^ 2;
t158 = qJD(1) ^ 2;
t231 = t145 * t158;
t230 = t148 * t149;
t229 = t149 * t150;
t225 = t150 * t155;
t222 = t153 * t158;
t216 = g(1) * t153 + g(3) * t224;
t215 = t157 * pkin(1) + t155 * qJ(2);
t144 = t149 ^ 2;
t146 = t152 ^ 2;
t213 = -t144 - t146;
t147 = t153 ^ 2;
t212 = t145 + t147;
t209 = qJD(1) * t152;
t207 = qJD(2) * t153;
t201 = qJDD(1) * t149;
t199 = t145 * qJDD(1);
t196 = t144 * t231;
t195 = t150 * t227;
t187 = t144 * t198;
t182 = t249 ^ 2;
t181 = t154 * t76 - t156 * t186;
t180 = t212 * t158;
t68 = t152 * t102 - t149 * t233;
t179 = 0.2e1 * t150 * t197;
t178 = qJ(3) * t224 + t157 * t241 + t215;
t177 = 0.2e1 * t212;
t95 = t153 * t220 - t219;
t54 = t95 * t148 - t151 * t225;
t56 = t98 * t148 - t151 * t224;
t175 = g(2) * t56 + g(3) * t54;
t94 = t153 * t221 + t218;
t174 = g(2) * t97 + g(3) * t94;
t61 = t153 * pkin(3) - t68;
t92 = -t149 * t206 + t152 * t207;
t171 = t156 * t10 + t154 * t14;
t11 = -t148 * t26 + t151 * t37;
t18 = -t148 * t40 + t151 * t62;
t27 = -t148 * t59 + t151 * t73;
t5 = t156 * t15 - t154 * t17;
t90 = -t153 * t148 + t151 * t226;
t22 = t89 * pkin(4) - t90 * pkin(6) + t61;
t24 = pkin(6) * t229 + t28;
t7 = -t154 * t24 + t156 * t22;
t8 = t154 * t22 + t156 * t24;
t170 = (-qJD(4) * t152 + qJD(2)) * t150;
t140 = t155 * pkin(1);
t169 = -t157 * qJ(2) + qJ(3) * t225 + t155 * t241 + t140;
t52 = t150 * t228 + t90 * t156;
t91 = t149 * t207 + t152 * t206;
t167 = qJD(1) * t91 - qJDD(1) * t68 - t32;
t166 = qJD(1) * t92 + qJDD(1) * t69 + t33;
t165 = t98 * pkin(3) + t97 * qJ(4) + t178;
t164 = -g(1) * t229 - g(2) * t94 + g(3) * t97;
t163 = t232 + t248;
t162 = t95 * pkin(3) + t94 * qJ(4) + t169;
t161 = t177 * t203 + t242;
t160 = t99 * t150 + (t203 + t204) * t145;
t136 = t147 * qJDD(1);
t118 = t144 * t199;
t79 = t152 * t190 - t189;
t70 = -t153 * qJD(4) + t92;
t64 = t96 * t210;
t63 = t168 * t210;
t55 = t148 * t225 + t95 * t151;
t51 = t90 * t154 - t195;
t47 = t52 * qJD(5);
t46 = -qJD(5) * t195 + t205 * t90;
t36 = t148 * t170 + t151 * t70;
t35 = t148 * t70 - t151 * t170;
t31 = t94 * t154 + t55 * t156;
t30 = -t55 * t154 + t94 * t156;
t23 = -pkin(4) * t229 - t27;
t21 = qJD(5) * t44 + t181;
t16 = -pkin(4) * t193 - t18;
t9 = -pkin(4) * t186 - t11;
t4 = -qJD(5) * t8 - t154 * t36 + t156 * t91;
t3 = qJD(5) * t7 + t154 * t91 + t156 * t36;
t1 = t5 * qJD(5) + t171;
t13 = [0, 0, 0, 0, 0, qJDD(1), t245, -t242 + t243, 0, 0, t199, t179, 0, t136, 0, 0, t163 * t153, -t163 * t150, t177 * t204 + t161 - t243, -t135 * pkin(1) - g(2) * t215 - g(3) * t140 + (t212 * t204 + t161) * qJ(2), t146 * t199, -0.2e1 * t152 * t149 * t199, -0.2e1 * t150 * t183, t118, t149 * t179, t136, -g(2) * t98 - g(3) * t95 + t149 * t160 + t153 * t167, t152 * t160 + t153 * t166 + t174, (-t149 * t166 + t152 * t167 + t245) * t150, t33 * t69 + t50 * t92 + t32 * t68 - t49 * t91 - g(2) * t178 - g(3) * t169 + (t99 * qJ(2) + t116 * qJD(2)) * t150, t76 * t90, -t90 * t75 - t76 * t89, (qJDD(1) * t90 + t76) * t229, t75 * t89, -0.2e1 * t75 * t229, t118, -g(2) * t57 - g(3) * t55 + t29 * t89 + t61 * t75 + t91 * t77 + (-qJD(1) * t35 + qJDD(1) * t27 + t11) * t229, t29 * t90 + t61 * t76 + t91 * t80 + (-qJD(1) * t36 - qJDD(1) * t28 - t12) * t229 + t175, -t11 * t90 - t12 * t89 - t27 * t76 - t28 * t75 + t35 * t80 - t36 * t77 - t174, -g(2) * t165 - g(3) * t162 + t11 * t27 + t12 * t28 - t18 * t35 + t19 * t36 + t29 * t61 + t39 * t91, -t20 * t52 - t44 * t46, t20 * t51 - t52 * t21 + t46 * t41 - t44 * t47, -t20 * t89 - t249 * t46 + t52 * t71, t21 * t51 + t41 * t47, -t21 * t89 - t249 * t47 - t51 * t71, t71 * t89, -g(2) * t246 - g(3) * t31 + t16 * t47 + t2 * t89 + t23 * t21 + t249 * t4 + t35 * t41 + t9 * t51 + t7 * t71, g(2) * t247 - g(3) * t30 - t1 * t89 - t16 * t46 - t23 * t20 - t249 * t3 + t35 * t44 + t9 * t52 - t8 * t71, -t1 * t51 - t2 * t52 + t7 * t20 - t8 * t21 - t3 * t41 - t4 * t44 + t5 * t46 - t6 * t47 - t175, t1 * t8 + t6 * t3 + t2 * t7 + t5 * t4 + t9 * t23 + t16 * t35 - g(2) * (t57 * pkin(4) + t56 * pkin(6) + t165) - g(3) * (t55 * pkin(4) + t54 * pkin(6) + t162); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t198, -t180, -qJ(2) * t180 - t248, 0, 0, 0, 0, 0, 0, -t149 * t180 - t183, -t152 * t180 + t184, t213 * t198, t33 * t149 + t32 * t152 + (-t116 * t150 + (t149 * t49 - t152 * t50) * t153) * qJD(1) - t245, 0, 0, 0, 0, 0, 0, -t148 * t187 - t152 * t75 + (t150 * t79 - t153 * t77) * t211, -t151 * t187 - t152 * t76 + (t150 * t82 - t153 * t80) * t211, t82 * t77 - t79 * t80 + (t148 * t76 - t151 * t75) * t149, -t29 * t152 + t18 * t79 - t19 * t82 + (-t11 * t148 + t12 * t151 - t208 * t39) * t149 - t245, 0, 0, 0, 0, 0, 0, -t168 * t71 + t21 * t230 - t239 * t249 - t79 * t41, -t20 * t230 + t238 * t249 - t79 * t44 - t96 * t71, -t168 * t20 - t96 * t21 + t238 * t41 + t239 * t44, t1 * t96 - t16 * t79 - t168 * t2 + t9 * t230 - t238 * t6 - t239 * t5 - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t152 * t222 + t201) * t150, (t149 * t222 + t200) * t150, t213 * t231, (-t243 + (t149 * t50 + t152 * t49) * qJD(1)) * t150 + t99 + t216, 0, 0, 0, 0, 0, 0, -t148 * t196 + (t151 * t201 - t209 * t77) * t150, -t151 * t196 + (-t148 * t201 - t209 * t80) * t150, -t148 * t75 - t151 * t76 + (t148 * t80 - t151 * t77) * t193, t11 * t151 + t12 * t148 + (-t243 + (-t152 * t39 + (-t148 * t18 + t151 * t19) * t149) * qJD(1)) * t150 + t216, 0, 0, 0, 0, 0, 0, -t151 * t21 - t63 * t249 + (-qJD(5) * t156 * t249 + t193 * t41 - t237) * t148, t151 * t20 - t64 * t249 + (t193 * t44 + t205 * t249 - t236) * t148, -t64 * t41 + t63 * t44 + (-t154 * t20 - t156 * t21 + (t154 * t41 + t156 * t44) * qJD(5)) * t148, -g(2) * t225 - t9 * t151 - t5 * t63 + t6 * t64 + (t16 * t193 + t1 * t156 - t154 * t2 + (-t154 * t6 - t156 * t5) * qJD(5)) * t148 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 * t80 + t75, -t117 + (-t211 * t77 + t185) * t150, -t77 ^ 2 - t80 ^ 2, t18 * t80 + t19 * t77 + t164 + t29, 0, 0, 0, 0, 0, 0, -t154 * t182 - t80 * t41 + t236, -t156 * t182 - t80 * t44 - t237, -t250 * t156 + (t249 * t44 - t21) * t154, -t16 * t80 + t251 * t156 + (-t249 * t5 + t1) * t154 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, -t41 ^ 2 + t44 ^ 2, t250, -t240, t244 * t44 - t181, t71, g(1) * t51 - g(2) * t30 - g(3) * t247 - t16 * t44 + t251, g(1) * t52 + g(2) * t31 - g(3) * t246 + t16 * t41 + t244 * t5 - t171, 0, 0;];
tau_reg = t13;
