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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:46
% EndTime: 2022-01-23 08:59:52
% DurationCPUTime: 3.04s
% Computational Cost: add. (2911->347), mult. (7255->530), div. (0->0), fcn. (5627->10), ass. (0->184)
t142 = sin(pkin(9));
t145 = cos(pkin(9));
t147 = cos(pkin(7));
t199 = qJD(1) * qJD(4);
t143 = sin(pkin(8));
t146 = cos(pkin(8));
t194 = t147 * qJDD(1);
t180 = t146 * t194;
t200 = qJD(1) * qJD(2);
t185 = t147 * t200;
t144 = sin(pkin(7));
t178 = t144 * qJ(3) + pkin(1);
t102 = pkin(2) * t147 + t178;
t203 = qJD(3) * t144;
t64 = -qJD(1) * t203 - qJDD(1) * t102 + qJDD(2);
t31 = qJ(2) * t180 + t143 * t64 + t146 * t185;
t26 = (-qJ(4) * qJDD(1) - t199) * t147 + t31;
t226 = qJ(4) * t146;
t166 = pkin(3) * t143 - t226;
t195 = t144 * qJDD(1);
t98 = qJ(2) * t195 + t144 * t200 + qJDD(3);
t36 = (qJDD(1) * t166 - t146 * t199) * t144 + t98;
t12 = t142 * t36 + t145 * t26;
t183 = t143 * t195;
t10 = pkin(6) * t183 + t12;
t182 = t143 * t194;
t30 = -qJ(2) * t182 - t143 * t185 + t146 * t64;
t29 = pkin(3) * t194 + qJDD(4) - t30;
t218 = t144 * t146;
t87 = t142 * t218 + t145 * t147;
t74 = qJDD(1) * t87;
t117 = t142 * t194;
t197 = qJDD(1) * t146;
t181 = t145 * t197;
t75 = t144 * t181 - t117;
t14 = pkin(4) * t74 - pkin(6) * t75 + t29;
t148 = sin(qJ(5));
t150 = cos(qJ(5));
t205 = qJD(1) * t147;
t191 = qJ(2) * t205;
t86 = -qJD(1) * t102 + qJD(2);
t48 = -t143 * t191 + t146 * t86;
t38 = pkin(3) * t205 + qJD(4) - t48;
t76 = t87 * qJD(1);
t207 = qJD(1) * t144;
t186 = t145 * t207;
t187 = t142 * t205;
t79 = t146 * t186 - t187;
t15 = pkin(4) * t76 - pkin(6) * t79 + t38;
t49 = t143 * t86 + t146 * t191;
t39 = -qJ(4) * t205 + t49;
t116 = qJ(2) * t207 + qJD(3);
t60 = t166 * t207 + t116;
t19 = t142 * t60 + t145 * t39;
t190 = t143 * t207;
t17 = pkin(6) * t190 + t19;
t6 = t148 * t15 + t150 * t17;
t2 = -qJD(5) * t6 - t148 * t10 + t150 * t14;
t240 = qJD(5) + t76;
t243 = t6 * t240 + t2;
t164 = t148 * t17 - t150 * t15;
t1 = -t164 * qJD(5) + t150 * t10 + t148 * t14;
t242 = t164 * t240 + t1;
t208 = qJD(1) * t143;
t188 = t150 * t208;
t171 = t144 * t188;
t202 = qJD(5) * t148;
t20 = -qJD(5) * t171 - t148 * t183 - t150 * t75 + t202 * t79;
t40 = t148 * t79 - t171;
t241 = t240 * t40 - t20;
t149 = sin(qJ(1));
t151 = cos(qJ(1));
t177 = g(1) * t149 - g(2) * t151;
t220 = t143 * t148;
t215 = t146 * t147;
t222 = t142 * t144;
t89 = t145 * t215 + t222;
t161 = t147 * t220 + t150 * t89;
t219 = t143 * t150;
t93 = t145 * t219 - t146 * t148;
t239 = t149 * t161 - t151 * t93;
t238 = -0.2e1 * t146;
t225 = qJDD(1) * pkin(1);
t237 = t225 + t177;
t52 = t147 * t219 - t148 * t89;
t92 = t145 * t220 + t146 * t150;
t236 = -t149 * t92 + t151 * t52;
t201 = qJ(2) * qJDD(1);
t189 = t148 * t208;
t43 = t144 * t189 + t150 * t79;
t234 = t43 * t40;
t233 = t143 * qJ(4) + pkin(2);
t81 = t89 * qJD(1);
t232 = t93 * qJD(5) + t147 * t188 - t148 * t81;
t231 = t92 * qJD(5) + t147 * t189 + t150 * t81;
t67 = qJ(2) * t215 - t143 * t102;
t57 = -qJ(4) * t147 + t67;
t162 = qJ(2) + t166;
t72 = t162 * t144;
t28 = t142 * t72 + t145 * t57;
t69 = qJDD(5) + t74;
t230 = t148 * t69;
t228 = t150 * t69;
t139 = t144 ^ 2;
t152 = qJD(1) ^ 2;
t224 = t139 * t152;
t223 = t142 * t143;
t221 = t143 * t144;
t217 = t145 * t149;
t216 = t145 * t151;
t214 = t147 * t149;
t213 = t147 * t151;
t212 = t147 * t152;
t138 = t143 ^ 2;
t140 = t146 ^ 2;
t210 = -t138 - t140;
t141 = t147 ^ 2;
t209 = t139 + t141;
t206 = qJD(1) * t146;
t204 = qJD(2) * t147;
t198 = qJDD(1) * t143;
t196 = t139 * qJDD(1);
t193 = t138 * t224;
t192 = t144 * t219;
t184 = t138 * t195;
t179 = t144 * t194;
t176 = t240 ^ 2;
t175 = t148 * t75 - t150 * t183;
t174 = t209 * t152;
t66 = -t143 * t147 * qJ(2) - t102 * t146;
t173 = 0.2e1 * t179;
t172 = 0.2e1 * t209;
t95 = t143 * t151 - t146 * t214;
t97 = t143 * t149 + t146 * t213;
t170 = g(1) * (t142 * t95 + t144 * t217) + g(2) * (t142 * t97 - t144 * t216);
t94 = t143 * t214 + t146 * t151;
t96 = t143 * t213 - t146 * t149;
t169 = -g(1) * t94 + g(2) * t96;
t168 = g(1) * t151 + g(2) * t149;
t59 = t147 * pkin(3) - t66;
t91 = -t143 * t203 + t146 * t204;
t11 = -t142 * t26 + t145 * t36;
t18 = -t142 * t39 + t145 * t60;
t27 = -t142 * t57 + t145 * t72;
t88 = -t147 * t142 + t145 * t218;
t22 = pkin(4) * t87 - pkin(6) * t88 + t59;
t24 = pkin(6) * t221 + t28;
t7 = -t148 * t24 + t150 * t22;
t8 = t148 * t22 + t150 * t24;
t163 = (-qJD(4) * t146 + qJD(2)) * t144;
t51 = t144 * t220 + t150 * t88;
t105 = pkin(4) * t145 + pkin(6) * t142 + pkin(3);
t160 = -t105 * t143 - qJ(2) + t226;
t90 = t143 * t204 + t146 * t203;
t159 = qJD(1) * t90 - qJDD(1) * t66 - t30;
t158 = qJD(1) * t91 + qJDD(1) * t67 + t31;
t157 = t172 * t200;
t156 = -g(1) * t96 - g(2) * t94 - g(3) * t221;
t130 = qJDD(2) - t225;
t155 = -t130 + t237;
t154 = t144 * t98 + (t200 + t201) * t139;
t136 = g(3) * t147;
t135 = t151 * qJ(2);
t133 = t149 * qJ(2);
t131 = t141 * qJDD(1);
t118 = t138 * t196;
t78 = t146 * t187 - t186;
t71 = (pkin(3) * t146 + t233) * t147 + t178;
t68 = -qJD(4) * t147 + t91;
t62 = t93 * t207;
t61 = t92 * t207;
t50 = t148 * t88 - t192;
t46 = t51 * qJD(5);
t45 = -qJD(5) * t192 + t202 * t88;
t35 = t142 * t163 + t145 * t68;
t34 = t142 * t68 - t145 * t163;
t33 = (t105 * t146 + t233) * t147 + pkin(1) + (pkin(4) * t142 - pkin(6) * t145 + qJ(3)) * t144;
t23 = -pkin(4) * t221 - t27;
t21 = qJD(5) * t43 + t175;
t16 = -pkin(4) * t190 - t18;
t9 = -pkin(4) * t183 - t11;
t4 = -qJD(5) * t8 - t148 * t35 + t150 * t90;
t3 = qJD(5) * t7 + t148 * t90 + t150 * t35;
t5 = [0, 0, 0, 0, 0, qJDD(1), t177, t168, 0, 0, t196, t173, 0, t131, 0, 0, t155 * t147, -t155 * t144, t172 * t201 + t157 - t168, -t130 * pkin(1) - g(1) * (-pkin(1) * t149 + t135) - g(2) * (pkin(1) * t151 + t133) + (t209 * t201 + t157) * qJ(2), t140 * t196, t143 * t196 * t238, t179 * t238, t118, t143 * t173, t131, -g(1) * t95 - g(2) * t97 + t143 * t154 + t147 * t159, t146 * t154 + t147 * t158 + t169, (-t143 * t158 + t146 * t159 + t177) * t144, t31 * t67 + t49 * t91 + t30 * t66 - t48 * t90 - g(1) * (-t102 * t149 + t135) - g(2) * (t102 * t151 + t133) + (t98 * qJ(2) + t116 * qJD(2)) * t144, t75 * t88, -t74 * t88 - t75 * t87, (qJDD(1) * t88 + t75) * t221, t74 * t87, -0.2e1 * t74 * t221, t118, t90 * t76 + t59 * t74 + t29 * t87 - g(1) * (t145 * t95 - t149 * t222) - g(2) * (t145 * t97 + t151 * t222) + (-qJD(1) * t34 + qJDD(1) * t27 + t11) * t221, t29 * t88 + t59 * t75 + t79 * t90 + (-qJD(1) * t35 - qJDD(1) * t28 - t12) * t221 + t170, -t11 * t88 - t12 * t87 - t27 * t75 - t28 * t74 + t34 * t79 - t35 * t76 - t169, t12 * t28 + t19 * t35 + t11 * t27 - t18 * t34 + t29 * t59 + t38 * t90 - g(1) * (-t149 * t71 + t151 * t162) - g(2) * (t149 * t162 + t151 * t71), -t20 * t51 - t43 * t45, t20 * t50 - t21 * t51 + t40 * t45 - t43 * t46, -t20 * t87 - t240 * t45 + t51 * t69, t21 * t50 + t40 * t46, -t21 * t87 - t240 * t46 - t50 * t69, t69 * t87, t4 * t240 + t7 * t69 + t2 * t87 + t34 * t40 + t23 * t21 + t9 * t50 + t16 * t46 + g(1) * t239 - g(2) * ((t143 * t217 + t151 * t89) * t150 + t96 * t148), -t3 * t240 - t8 * t69 - t1 * t87 + t34 * t43 - t23 * t20 + t9 * t51 - t16 * t45 - g(1) * (-t149 * t52 - t151 * t92) - g(2) * t236, -t1 * t50 - t164 * t45 - t2 * t51 + t20 * t7 - t21 * t8 - t3 * t40 - t4 * t43 - t46 * t6 - t170, t1 * t8 + t6 * t3 + t2 * t7 - t164 * t4 + t9 * t23 + t16 * t34 - g(1) * (-t33 * t149 - t151 * t160) - g(2) * (-t149 * t160 + t33 * t151); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, t195, -t174, -qJ(2) * t174 + qJDD(2) - t237, 0, 0, 0, 0, 0, 0, -t143 * t174 - t180, -t146 * t174 + t182, t210 * t195, t143 * t31 + t146 * t30 + (-t116 * t144 + (t143 * t48 - t146 * t49) * t147) * qJD(1) - t177, 0, 0, 0, 0, 0, 0, -t142 * t184 - t146 * t74 + (t144 * t78 - t147 * t76) * t208, -t145 * t184 - t146 * t75 + (t144 * t81 - t147 * t79) * t208, t76 * t81 - t78 * t79 + (t142 * t75 - t145 * t74) * t143, -t146 * t29 + t18 * t78 - t19 * t81 + (-t11 * t142 + t12 * t145 - t205 * t38) * t143 - t177, 0, 0, 0, 0, 0, 0, t21 * t223 - t232 * t240 - t40 * t78 - t69 * t92, -t20 * t223 + t231 * t240 - t43 * t78 - t69 * t93, -t20 * t92 - t21 * t93 + t231 * t40 + t232 * t43, t1 * t93 - t16 * t78 + t164 * t232 - t2 * t92 + t9 * t223 - t231 * t6 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t146 * t212 + t198) * t144, (t143 * t212 + t197) * t144, t210 * t224, t136 + ((t143 * t49 + t146 * t48) * qJD(1) - t168) * t144 + t98, 0, 0, 0, 0, 0, 0, -t142 * t193 + (t145 * t198 - t206 * t76) * t144, -t145 * t193 + (-t142 * t198 - t206 * t79) * t144, -t142 * t74 - t145 * t75 + (t142 * t79 - t145 * t76) * t190, t11 * t145 + t12 * t142 + t136 + ((-t146 * t38 + (-t142 * t18 + t145 * t19) * t143) * qJD(1) - t168) * t144, 0, 0, 0, 0, 0, 0, -t145 * t21 - t61 * t240 + (-qJD(5) * t150 * t240 + t190 * t40 - t230) * t142, t145 * t20 - t62 * t240 + (t190 * t43 + t202 * t240 - t228) * t142, -t40 * t62 + t43 * t61 + (-t148 * t20 - t150 * t21 + (t148 * t40 + t150 * t43) * qJD(5)) * t142, -t145 * t9 + t164 * t61 + t6 * t62 + t136 - t168 * t144 + (t16 * t190 + t1 * t150 - t148 * t2 + (-t148 * t6 + t150 * t164) * qJD(5)) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190 * t79 + t74, -t117 + (-t208 * t76 + t181) * t144, -t76 ^ 2 - t79 ^ 2, t18 * t79 + t19 * t76 + t156 + t29, 0, 0, 0, 0, 0, 0, -t148 * t176 - t40 * t79 + t228, -t150 * t176 - t43 * t79 - t230, -t241 * t150 + (t240 * t43 - t21) * t148, t242 * t148 + t150 * t243 - t16 * t79 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, -t40 ^ 2 + t43 ^ 2, t241, -t234, -t175 + (-qJD(5) + t240) * t43, t69, -t16 * t43 - g(1) * t236 - g(2) * (-(-t143 * t216 + t89 * t149) * t148 + t94 * t150) + g(3) * t50 + t243, t16 * t40 - g(1) * (-t149 * t93 - t151 * t161) + g(2) * t239 + g(3) * t51 - t242, 0, 0;];
tau_reg = t5;
