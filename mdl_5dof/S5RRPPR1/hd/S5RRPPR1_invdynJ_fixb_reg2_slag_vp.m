% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:04
% EndTime: 2020-01-03 11:56:06
% DurationCPUTime: 1.49s
% Computational Cost: add. (2420->237), mult. (3882->294), div. (0->0), fcn. (2553->16), ass. (0->161)
t145 = qJ(1) + qJ(2);
t133 = pkin(8) + t145;
t118 = sin(t133);
t119 = cos(t133);
t223 = g(2) * t119 + g(3) * t118;
t140 = qJDD(1) + qJDD(2);
t147 = sin(pkin(8));
t149 = cos(pkin(8));
t155 = cos(qJ(2));
t210 = pkin(1) * qJD(2);
t191 = qJD(1) * t210;
t152 = sin(qJ(2));
t198 = qJDD(1) * t152;
t225 = pkin(1) * t198 + t155 * t191;
t216 = t155 * pkin(1);
t126 = qJDD(1) * t216;
t66 = t140 * pkin(2) - t152 * t191 + t126;
t194 = t147 * t225 - t149 * t66;
t178 = qJDD(4) + t194;
t33 = -t140 * pkin(3) + t178;
t226 = t33 + t223;
t211 = pkin(1) * qJD(1);
t193 = t152 * t211;
t102 = t147 * t193;
t192 = t155 * t211;
t71 = t149 * t192 - t102;
t201 = qJD(4) - t71;
t148 = cos(pkin(9));
t154 = cos(qJ(5));
t202 = t154 * t148;
t146 = sin(pkin(9));
t151 = sin(qJ(5));
t205 = t151 * t146;
t80 = -t202 + t205;
t81 = t154 * t146 + t151 * t148;
t144 = qJD(1) + qJD(2);
t63 = t81 * t144;
t141 = t146 ^ 2;
t142 = t148 ^ 2;
t199 = t141 + t142;
t224 = t144 * t199;
t222 = g(2) * t118 - g(3) * t119;
t221 = t63 ^ 2;
t220 = pkin(1) * t152;
t219 = t147 * pkin(2);
t218 = t148 * pkin(4);
t136 = t148 * pkin(7);
t217 = t149 * pkin(2);
t190 = t144 * t205;
t61 = -t144 * t202 + t190;
t215 = t63 * t61;
t116 = qJ(4) + t219;
t77 = (-pkin(7) - t116) * t146;
t78 = t148 * t116 + t136;
t45 = -t151 * t78 + t154 * t77;
t214 = t45 * qJD(5) - t201 * t80;
t46 = t151 * t77 + t154 * t78;
t213 = -t46 * qJD(5) - t201 * t81;
t177 = t80 * t140;
t76 = t81 * qJD(5);
t39 = t144 * t76 + t177;
t186 = qJD(5) * t202;
t187 = qJD(5) * t205;
t75 = -t186 + t187;
t212 = -t81 * t39 + t75 * t61;
t206 = t149 * t152;
t159 = pkin(1) * (t147 * t155 + t206);
t69 = qJD(1) * t159;
t209 = t69 * t144;
t70 = qJD(2) * t159;
t208 = t70 * t144;
t125 = pkin(2) + t216;
t74 = pkin(1) * t206 + t147 * t125;
t41 = t147 * t66 + t149 * t225;
t31 = t140 * qJ(4) + t144 * qJD(4) + t41;
t21 = t146 * qJDD(3) + t148 * t31;
t86 = t144 * pkin(2) + t192;
t56 = t147 * t86 + t149 * t193;
t52 = t144 * qJ(4) + t56;
t43 = t146 * qJD(3) + t148 * t52;
t207 = t148 * t140;
t197 = t226 * t146;
t196 = t140 * t81 + t144 * t186;
t114 = t147 * t220;
t120 = pkin(3) + t218;
t134 = sin(t145);
t123 = pkin(2) * t134;
t150 = -pkin(7) - qJ(4);
t195 = t118 * t120 + t119 * t150 + t123;
t135 = cos(t145);
t124 = pkin(2) * t135;
t188 = t119 * pkin(3) + t118 * qJ(4) + t124;
t185 = g(2) * t134 - g(3) * t135;
t128 = t148 * qJDD(3);
t15 = t128 + (-pkin(7) * t140 - t31) * t146;
t16 = pkin(7) * t207 + t21;
t184 = t154 * t15 - t151 * t16;
t143 = pkin(9) + qJ(5);
t131 = sin(t143);
t24 = -t120 * t140 + t178;
t55 = t149 * t86 - t102;
t169 = qJD(4) - t55;
t44 = -t120 * t144 + t169;
t183 = t223 * t131 + t24 * t81 - t44 * t75;
t182 = t199 * t140;
t73 = t149 * t125 - t114;
t181 = qJD(1) * (-qJD(2) + t144);
t180 = qJD(2) * (-qJD(1) - t144);
t68 = -pkin(3) - t73;
t176 = -t118 * t150 + t119 * t120 + t124;
t173 = -g(2) * t135 - g(3) * t134;
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t172 = -g(2) * t156 - g(3) * t153;
t38 = t144 * t187 - t196;
t171 = -t80 * t38 + t63 * t76;
t170 = t118 * pkin(3) - t119 * qJ(4) + t123;
t168 = -t41 + t222;
t167 = -t194 - t223;
t130 = t148 * qJD(3);
t36 = t130 + (-pkin(7) * t144 - t52) * t146;
t37 = t144 * t136 + t43;
t10 = -t151 * t37 + t154 * t36;
t11 = t151 * t36 + t154 * t37;
t164 = t151 * t15 + t154 * t16;
t3 = t10 * qJD(5) + t164;
t4 = -t11 * qJD(5) + t184;
t166 = t10 * t75 - t11 * t76 - t3 * t80 - t4 * t81 - t222;
t165 = t140 * t68 + t208;
t67 = qJ(4) + t74;
t53 = (-pkin(7) - t67) * t146;
t54 = t148 * t67 + t136;
t22 = -t151 * t54 + t154 * t53;
t23 = t151 * t53 + t154 * t54;
t121 = -pkin(3) - t217;
t163 = t121 * t140 - t209;
t20 = -t146 * t31 + t128;
t162 = -t20 * t146 + t21 * t148 - t222;
t72 = t149 * t155 * t210 - qJD(2) * t114;
t160 = t126 + t173;
t132 = cos(t143);
t158 = -t132 * t223 + t24 * t80 + t44 * t76;
t138 = t156 * pkin(1);
t137 = t153 * pkin(1);
t113 = t142 * t140;
t112 = t141 * t140;
t89 = -t120 - t217;
t87 = 0.2e1 * t146 * t207;
t65 = qJD(4) + t72;
t60 = t61 ^ 2;
t57 = t68 - t218;
t51 = -t144 * pkin(3) + t169;
t49 = -t76 * qJD(5) - t80 * qJDD(5);
t48 = -t75 * qJD(5) + t81 * qJDD(5);
t42 = -t146 * t52 + t130;
t13 = t39 * t80 + t61 * t76;
t12 = -t38 * t81 - t63 * t75;
t9 = -t23 * qJD(5) - t81 * t65;
t8 = t22 * qJD(5) - t80 * t65;
t5 = -t171 + t212;
t1 = [0, 0, 0, 0, 0, qJDD(1), t172, g(2) * t153 - g(3) * t156, 0, 0, 0, 0, 0, 0, 0, t140, (t140 * t155 + t152 * t180) * pkin(1) + t160, ((-qJDD(1) - t140) * t152 + t155 * t180) * pkin(1) + t185, 0, (t172 + (t152 ^ 2 + t155 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t140, t73 * t140 + t167 - t208, -t74 * t140 - t72 * t144 + t168, 0, t41 * t74 + t56 * t72 - t194 * t73 - t55 * t70 - g(2) * (t124 + t138) - g(3) * (t123 + t137), t112, t87, 0, t113, 0, 0, (-t226 - t165) * t148, t146 * t165 + t197, t67 * t182 + t65 * t224 + t162, t33 * t68 + t51 * t70 - g(2) * (t138 + t188) - g(3) * (t137 + t170) + (t21 * t67 + t43 * t65) * t148 + (-t20 * t67 - t42 * t65) * t146, t12, t5, t48, t13, t49, 0, t9 * qJD(5) + t22 * qJDD(5) + t57 * t39 + t70 * t61 + t158, -t8 * qJD(5) - t23 * qJDD(5) - t57 * t38 + t70 * t63 + t183, t22 * t38 - t23 * t39 - t8 * t61 - t9 * t63 + t166, t3 * t23 + t11 * t8 + t4 * t22 + t10 * t9 + t24 * t57 + t44 * t70 - g(2) * (t138 + t176) - g(3) * (t137 + t195); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t181 * t220 + t160, (t155 * t181 - t198) * pkin(1) + t185, 0, 0, 0, 0, 0, 0, 0, t140, t140 * t217 + t167 + t209, -t140 * t219 + t71 * t144 + t168, 0, t55 * t69 - t56 * t71 + (t147 * t41 - t149 * t194 + t173) * pkin(2), t112, t87, 0, t113, 0, 0, (-t226 - t163) * t148, t146 * t163 + t197, t116 * t182 + t201 * t224 + t162, t33 * t121 - t51 * t69 - g(2) * t188 - g(3) * t170 + (t21 * t116 + t201 * t43) * t148 + (-t20 * t116 - t201 * t42) * t146, t12, t5, t48, t13, t49, 0, t213 * qJD(5) + t45 * qJDD(5) + t89 * t39 - t69 * t61 + t158, -t214 * qJD(5) - t46 * qJDD(5) - t89 * t38 - t69 * t63 + t183, -t213 * t63 - t214 * t61 + t45 * t38 - t46 * t39 + t166, -g(2) * t176 - g(3) * t195 + t213 * t10 + t214 * t11 + t24 * t89 + t3 * t46 + t4 * t45 - t44 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t146 + t20 * t148 - g(1), 0, 0, 0, 0, 0, 0, t49, -t48, t171 + t212, -t10 * t76 - t11 * t75 + t3 * t81 - t4 * t80 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, t146 * t140, -t199 * t144 ^ 2, (t42 * t146 - t43 * t148) * t144 + t226, 0, 0, 0, 0, 0, 0, 0.2e1 * t63 * qJD(5) + t177, (-t61 - t190) * qJD(5) + t196, -t60 - t221, t10 * t63 + t11 * t61 + t223 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t60 + t221, (t61 - t190) * qJD(5) + t196, -t215, -t177, qJDD(5), -g(1) * t132 + t131 * t222 - t44 * t63 + t184, g(1) * t131 + t132 * t222 + t44 * t61 - t164, 0, 0;];
tau_reg = t1;
