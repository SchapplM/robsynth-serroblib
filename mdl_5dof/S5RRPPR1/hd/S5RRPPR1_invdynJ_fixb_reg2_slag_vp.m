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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:18:35
% EndTime: 2019-12-05 18:18:37
% DurationCPUTime: 1.34s
% Computational Cost: add. (2420->241), mult. (3882->295), div. (0->0), fcn. (2553->16), ass. (0->161)
t141 = sin(pkin(8));
t146 = sin(qJ(2));
t202 = pkin(1) * qJD(1);
t186 = t146 * t202;
t100 = t141 * t186;
t143 = cos(pkin(8));
t149 = cos(qJ(2));
t185 = t149 * t202;
t71 = t143 * t185 - t100;
t192 = qJD(4) - t71;
t139 = qJ(1) + qJ(2);
t129 = pkin(8) + t139;
t114 = sin(t129);
t115 = cos(t129);
t221 = g(2) * t115 + g(3) * t114;
t201 = pkin(1) * qJD(2);
t184 = qJD(1) * t201;
t189 = qJDD(1) * t146;
t220 = pkin(1) * t189 + t149 * t184;
t142 = cos(pkin(9));
t148 = cos(qJ(5));
t193 = t148 * t142;
t140 = sin(pkin(9));
t145 = sin(qJ(5));
t198 = t140 * t145;
t80 = -t193 + t198;
t81 = t140 * t148 + t142 * t145;
t138 = qJD(1) + qJD(2);
t63 = t81 * t138;
t135 = t140 ^ 2;
t136 = t142 ^ 2;
t190 = t135 + t136;
t219 = t138 * t190;
t191 = g(2) * t114 - g(3) * t115;
t130 = sin(t139);
t131 = cos(t139);
t218 = g(2) * t131 + g(3) * t130;
t217 = t63 ^ 2;
t216 = pkin(1) * t146;
t147 = sin(qJ(1));
t215 = pkin(1) * t147;
t214 = pkin(1) * t149;
t150 = cos(qJ(1));
t213 = pkin(1) * t150;
t212 = pkin(2) * t130;
t211 = pkin(2) * t131;
t210 = pkin(2) * t141;
t209 = pkin(2) * t143;
t208 = pkin(4) * t142;
t132 = t142 * pkin(7);
t183 = t138 * t198;
t61 = -t138 * t193 + t183;
t207 = t63 * t61;
t112 = qJ(4) + t210;
t77 = (-pkin(7) - t112) * t140;
t78 = t112 * t142 + t132;
t45 = -t145 * t78 + t148 * t77;
t206 = t45 * qJD(5) - t192 * t80;
t46 = t145 * t77 + t148 * t78;
t205 = -t46 * qJD(5) - t192 * t81;
t134 = qJDD(1) + qJDD(2);
t170 = t80 * t134;
t76 = t81 * qJD(5);
t39 = t138 * t76 + t170;
t179 = qJD(5) * t193;
t180 = qJD(5) * t198;
t75 = -t179 + t180;
t204 = -t81 * t39 + t75 * t61;
t203 = t221 * t142;
t194 = t143 * t146;
t155 = pkin(1) * (t141 * t149 + t194);
t69 = qJD(1) * t155;
t200 = t138 * t69;
t70 = qJD(2) * t155;
t199 = t138 * t70;
t121 = pkin(2) + t214;
t74 = pkin(1) * t194 + t141 * t121;
t122 = qJDD(1) * t214;
t66 = pkin(2) * t134 - t146 * t184 + t122;
t41 = t141 * t66 + t220 * t143;
t31 = qJ(4) * t134 + qJD(4) * t138 + t41;
t21 = t140 * qJDD(3) + t142 * t31;
t84 = pkin(2) * t138 + t185;
t56 = t141 * t84 + t143 * t186;
t52 = qJ(4) * t138 + t56;
t43 = t140 * qJD(3) + t142 * t52;
t196 = t142 * t134;
t188 = t81 * t134 + t138 * t179;
t187 = t220 * t141 - t143 * t66;
t110 = t141 * t216;
t181 = t122 + t218;
t116 = pkin(3) + t208;
t178 = -g(2) * t130 + g(3) * t131;
t124 = t142 * qJDD(3);
t15 = t124 + (-pkin(7) * t134 - t31) * t140;
t16 = pkin(7) * t196 + t21;
t177 = -t145 * t16 + t148 * t15;
t137 = pkin(9) + qJ(5);
t128 = cos(t137);
t172 = qJDD(4) + t187;
t24 = -t116 * t134 + t172;
t55 = t143 * t84 - t100;
t163 = qJD(4) - t55;
t44 = -t116 * t138 + t163;
t176 = t221 * t128 + t24 * t80 + t44 * t76;
t175 = t190 * t134;
t73 = t121 * t143 - t110;
t174 = qJD(1) * (-qJD(2) + t138);
t173 = qJD(2) * (-qJD(1) - t138);
t68 = -pkin(3) - t73;
t169 = -t41 - t191;
t168 = t187 - t221;
t165 = g(2) * t150 + g(3) * t147;
t38 = t138 * t180 - t188;
t164 = -t38 * t80 + t63 * t76;
t126 = t142 * qJD(3);
t36 = t126 + (-pkin(7) * t138 - t52) * t140;
t37 = t138 * t132 + t43;
t10 = -t145 * t37 + t148 * t36;
t11 = t145 * t36 + t148 * t37;
t160 = t145 * t15 + t148 * t16;
t3 = t10 * qJD(5) + t160;
t4 = -t11 * qJD(5) + t177;
t162 = t10 * t75 - t11 * t76 - t3 * t80 - t4 * t81 + t191;
t161 = -t134 * t68 - t199;
t67 = qJ(4) + t74;
t53 = (-pkin(7) - t67) * t140;
t54 = t142 * t67 + t132;
t22 = -t145 * t54 + t148 * t53;
t23 = t145 * t53 + t148 * t54;
t117 = -pkin(3) - t209;
t159 = -t117 * t134 + t200;
t20 = -t140 * t31 + t124;
t158 = -t20 * t140 + t21 * t142 + t191;
t72 = t143 * t149 * t201 - qJD(2) * t110;
t157 = -pkin(3) * t114 + t115 * qJ(4) - t212;
t144 = -pkin(7) - qJ(4);
t156 = t114 * t144 - t115 * t116 - t211;
t33 = -t134 * pkin(3) + t172;
t154 = -pkin(3) * t115 - qJ(4) * t114 - t211;
t153 = -t114 * t116 - t115 * t144 - t212;
t127 = sin(t137);
t152 = -t127 * t221 + t24 * t81 - t44 * t75;
t109 = t136 * t134;
t108 = t135 * t134;
t87 = -t116 - t209;
t85 = 0.2e1 * t140 * t196;
t65 = qJD(4) + t72;
t60 = t61 ^ 2;
t57 = t68 - t208;
t51 = -pkin(3) * t138 + t163;
t49 = -qJD(5) * t76 - qJDD(5) * t80;
t48 = -qJD(5) * t75 + qJDD(5) * t81;
t42 = -t140 * t52 + t126;
t32 = t33 * t140;
t13 = t39 * t80 + t61 * t76;
t12 = -t38 * t81 - t63 * t75;
t9 = -t23 * qJD(5) - t81 * t65;
t8 = t22 * qJD(5) - t80 * t65;
t5 = -t164 + t204;
t1 = [0, 0, 0, 0, 0, qJDD(1), t165, -g(2) * t147 + g(3) * t150, 0, 0, 0, 0, 0, 0, 0, t134, (t134 * t149 + t146 * t173) * pkin(1) + t181, ((-qJDD(1) - t134) * t146 + t149 * t173) * pkin(1) + t178, 0, (t165 + (t146 ^ 2 + t149 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t134, t134 * t73 - t168 - t199, -t134 * t74 - t138 * t72 + t169, 0, t41 * t74 + t56 * t72 - t187 * t73 - t55 * t70 - g(2) * (-t211 - t213) - g(3) * (-t212 - t215), t108, t85, 0, t109, 0, 0, (t161 - t33) * t142 + t203, t32 + (-t161 - t221) * t140, t67 * t175 + t65 * t219 + t158, t33 * t68 + t51 * t70 - g(2) * (t154 - t213) - g(3) * (t157 - t215) + (t21 * t67 + t43 * t65) * t142 + (-t20 * t67 - t42 * t65) * t140, t12, t5, t48, t13, t49, 0, qJD(5) * t9 + qJDD(5) * t22 + t39 * t57 + t61 * t70 + t176, -t8 * qJD(5) - t23 * qJDD(5) - t57 * t38 + t70 * t63 + t152, t22 * t38 - t23 * t39 - t61 * t8 - t63 * t9 + t162, t3 * t23 + t11 * t8 + t4 * t22 + t10 * t9 + t24 * t57 + t44 * t70 - g(2) * (t156 - t213) - g(3) * (t153 - t215); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t174 * t216 + t181, (t149 * t174 - t189) * pkin(1) + t178, 0, 0, 0, 0, 0, 0, 0, t134, t134 * t209 - t168 + t200, -t134 * t210 + t138 * t71 + t169, 0, t55 * t69 - t56 * t71 + (t141 * t41 - t143 * t187 + t218) * pkin(2), t108, t85, 0, t109, 0, 0, (t159 - t33) * t142 + t203, t32 + (-t159 - t221) * t140, t112 * t175 + t192 * t219 + t158, t33 * t117 - t51 * t69 - g(2) * t154 - g(3) * t157 + (t21 * t112 + t192 * t43) * t142 + (-t20 * t112 - t192 * t42) * t140, t12, t5, t48, t13, t49, 0, t205 * qJD(5) + t45 * qJDD(5) + t87 * t39 - t69 * t61 + t176, -t206 * qJD(5) - t46 * qJDD(5) - t87 * t38 - t69 * t63 + t152, -t205 * t63 - t206 * t61 + t38 * t45 - t39 * t46 + t162, -g(2) * t156 - g(3) * t153 + t205 * t10 + t206 * t11 + t24 * t87 + t3 * t46 + t4 * t45 - t44 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t21 + t142 * t20 - g(1), 0, 0, 0, 0, 0, 0, t49, -t48, t164 + t204, -t10 * t76 - t11 * t75 + t3 * t81 - t4 * t80 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, t140 * t134, -t190 * t138 ^ 2, (t140 * t42 - t142 * t43) * t138 + t33 - t221, 0, 0, 0, 0, 0, 0, 0.2e1 * t63 * qJD(5) + t170, (-t61 - t183) * qJD(5) + t188, -t60 - t217, t10 * t63 + t11 * t61 - t221 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, -t60 + t217, (t61 - t183) * qJD(5) + t188, -t207, -t170, qJDD(5), -g(1) * t128 - t127 * t191 - t44 * t63 + t177, g(1) * t127 - t128 * t191 + t44 * t61 - t160, 0, 0;];
tau_reg = t1;
