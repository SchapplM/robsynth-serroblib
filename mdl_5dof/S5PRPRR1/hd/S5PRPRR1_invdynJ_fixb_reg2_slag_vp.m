% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:11
% EndTime: 2019-12-05 15:43:14
% DurationCPUTime: 1.49s
% Computational Cost: add. (2666->246), mult. (6079->311), div. (0->0), fcn. (4721->12), ass. (0->138)
t136 = sin(qJ(5));
t192 = cos(qJ(5));
t138 = cos(qJ(4));
t134 = cos(pkin(9));
t174 = qJD(2) * t134;
t161 = t138 * t174;
t137 = sin(qJ(4));
t133 = sin(pkin(9));
t175 = qJD(2) * t133;
t163 = t137 * t175;
t76 = -t161 + t163;
t87 = t133 * t138 + t134 * t137;
t78 = t87 * qJD(2);
t148 = t136 * t76 - t192 * t78;
t173 = qJD(4) * t137;
t162 = t133 * t173;
t167 = t134 * qJDD(2);
t168 = t133 * qJDD(2);
t165 = qJD(4) * t161 + t137 * t167 + t138 * t168;
t45 = qJD(2) * t162 - t165;
t152 = t137 * t168 - t138 * t167;
t81 = t87 * qJD(4);
t46 = qJD(2) * t81 + t152;
t142 = t148 * qJD(5) + t136 * t45 - t192 * t46;
t131 = qJD(4) + qJD(5);
t179 = t148 * t131;
t203 = t142 - t179;
t160 = qJD(5) * t192;
t171 = qJD(5) * t136;
t147 = -t136 * t46 - t76 * t160 - t78 * t171 - t192 * t45;
t40 = -t136 * t78 - t192 * t76;
t182 = t131 * t40;
t202 = t147 - t182;
t187 = t40 ^ 2;
t188 = t148 ^ 2;
t201 = -t187 + t188;
t186 = t40 * t148;
t178 = qJDD(2) * pkin(2);
t130 = pkin(8) + qJ(2);
t121 = sin(t130);
t123 = cos(t130);
t198 = g(1) * t121 - g(2) * t123;
t149 = -qJDD(3) + t178 + t198;
t156 = g(1) * t123 + g(2) * t121;
t94 = qJ(3) * t174 + t133 * qJD(1);
t72 = pkin(6) * t174 + t94;
t180 = t137 * t72;
t119 = t134 * qJD(1);
t184 = pkin(6) + qJ(3);
t97 = t184 * t133;
t71 = -qJD(2) * t97 + t119;
t33 = t138 * t71 - t180;
t27 = -pkin(7) * t78 + t33;
t24 = qJD(4) * pkin(4) + t27;
t34 = t137 * t71 + t138 * t72;
t28 = -pkin(7) * t76 + t34;
t117 = t134 * qJDD(1);
t169 = qJD(2) * qJD(3);
t58 = t117 + (-t184 * qJDD(2) - t169) * t133;
t70 = qJ(3) * t167 + t133 * qJDD(1) + t134 * t169;
t59 = pkin(6) * t167 + t70;
t157 = -t137 * t59 + t138 * t58;
t18 = -qJD(4) * t34 + t157;
t7 = qJDD(4) * pkin(4) + pkin(7) * t45 + t18;
t172 = qJD(4) * t138;
t166 = -t137 * t58 - t138 * t59 - t71 * t172;
t17 = -t72 * t173 - t166;
t8 = -pkin(7) * t46 + t17;
t1 = t136 * t7 + t24 * t160 - t28 * t171 + t192 * t8;
t129 = pkin(9) + qJ(4);
t124 = qJ(5) + t129;
t112 = sin(t124);
t113 = cos(t124);
t114 = pkin(3) * t134 + pkin(2);
t96 = -t114 * qJD(2) + qJD(3);
t52 = pkin(4) * t76 + t96;
t200 = g(3) * t112 + t156 * t113 - t52 * t40 - t1;
t164 = t192 * t28;
t10 = t136 * t24 + t164;
t2 = -t10 * qJD(5) - t136 * t8 + t192 * t7;
t199 = -g(3) * t113 + t156 * t112 + t52 * t148 + t2;
t197 = qJ(3) * qJDD(2) + t169;
t196 = t78 ^ 2;
t195 = pkin(4) * t46;
t194 = pkin(4) * t81;
t80 = -t134 * t172 + t162;
t177 = t138 * t134;
t86 = t133 * t137 - t177;
t21 = t136 * t81 + t86 * t160 + t87 * t171 + t192 * t80;
t50 = -t136 * t86 + t192 * t87;
t193 = t142 * t50 - t21 * t40;
t185 = t78 * t76;
t183 = -t87 * t46 + t80 * t76;
t98 = t184 * t134;
t57 = -t137 * t97 + t138 * t98;
t181 = t136 * t28;
t127 = t133 ^ 2;
t128 = t134 ^ 2;
t176 = t127 + t128;
t56 = -t137 * t98 - t138 * t97;
t22 = t50 * qJD(5) - t136 * t80 + t192 * t81;
t49 = t136 * t87 + t192 * t86;
t154 = t147 * t49 - t148 * t22;
t153 = -t45 * t86 + t78 * t81;
t125 = qJDD(4) + qJDD(5);
t151 = t125 * t50 - t131 * t21;
t150 = t133 * (-qJ(3) * t175 + t119) - t134 * t94;
t35 = -pkin(7) * t87 + t56;
t36 = -pkin(7) * t86 + t57;
t19 = -t136 * t36 + t192 * t35;
t20 = t136 * t35 + t192 * t36;
t95 = -t114 * qJDD(2) + qJDD(3);
t146 = t149 + t178;
t120 = sin(t129);
t122 = cos(t129);
t145 = -g(3) * t122 + t156 * t120;
t31 = qJD(3) * t177 - t97 * t172 + (-qJD(3) * t133 - qJD(4) * t98) * t137;
t69 = -t197 * t133 + t117;
t144 = -t133 * t69 + t134 * t70 - t156;
t143 = -t198 + t95;
t32 = -t87 * qJD(3) - t57 * qJD(4);
t132 = qJDD(1) - g(3);
t126 = -pkin(7) - t184;
t88 = pkin(4) * t122 + t114;
t74 = t76 ^ 2;
t62 = pkin(4) * t86 - t114;
t48 = -qJD(4) * t81 - qJDD(4) * t86;
t47 = -qJD(4) * t80 + qJDD(4) * t87;
t30 = t95 + t195;
t26 = pkin(7) * t80 + t32;
t25 = -pkin(7) * t81 + t31;
t13 = -t125 * t49 - t131 * t22;
t12 = t192 * t27 - t181;
t11 = -t136 * t27 - t164;
t9 = t192 * t24 - t181;
t4 = -t20 * qJD(5) - t136 * t25 + t192 * t26;
t3 = t19 * qJD(5) + t136 * t26 + t192 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t70 + t134 * t69 - g(3), 0, 0, 0, 0, 0, 0, t48, -t47, t153 + t183, t17 * t87 - t18 * t86 - t33 * t81 - t34 * t80 - g(3), 0, 0, 0, 0, 0, 0, t13, -t151, t154 + t193, t1 * t50 - t10 * t21 - t2 * t49 - t22 * t9 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t198, t156, 0, 0, t127 * qJDD(2), 0.2e1 * t133 * t167, 0, t128 * qJDD(2), 0, 0, t146 * t134, -t146 * t133, t197 * t176 + t144, t149 * pkin(2) + t144 * qJ(3) - t150 * qJD(3), -t45 * t87 - t78 * t80, -t153 + t183, t47, t46 * t86 + t76 * t81, t48, 0, qJD(4) * t32 + qJDD(4) * t56 - t114 * t46 + t122 * t198 + t81 * t96 + t86 * t95, -qJD(4) * t31 - qJDD(4) * t57 + t114 * t45 - t120 * t198 - t80 * t96 + t87 * t95, -t17 * t86 - t18 * t87 - t31 * t76 - t32 * t78 + t33 * t80 - t34 * t81 + t45 * t56 - t46 * t57 - t156, t17 * t57 + t34 * t31 + t18 * t56 + t33 * t32 - t95 * t114 - g(1) * (-t114 * t121 + t123 * t184) - g(2) * (t114 * t123 + t121 * t184), t147 * t50 + t148 * t21, -t154 + t193, t151, -t142 * t49 - t22 * t40, t13, 0, t113 * t198 + t125 * t19 + t131 * t4 - t142 * t62 - t194 * t40 + t22 * t52 + t30 * t49, -t112 * t198 - t125 * t20 - t131 * t3 + t147 * t62 - t148 * t194 - t21 * t52 + t30 * t50, -t1 * t49 - t10 * t22 + t142 * t20 - t147 * t19 + t148 * t4 - t2 * t50 + t21 * t9 + t3 * t40 - t156, t1 * t20 + t10 * t3 + t2 * t19 + t9 * t4 + t30 * t62 + t52 * t194 - g(1) * (-t121 * t88 - t123 * t126) - g(2) * (-t121 * t126 + t123 * t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t168, -t176 * qJD(2) ^ 2, t150 * qJD(2) - t149, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t78 + t152, (-t76 - t163) * qJD(4) + t165, -t74 - t196, t33 * t78 + t34 * t76 + t143, 0, 0, 0, 0, 0, 0, -t142 - t179, t147 + t182, -t187 - t188, -t10 * t40 - t148 * t9 + t143 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, -t74 + t196, (t76 - t163) * qJD(4) + t165, -t185, -t152, qJDD(4), -t78 * t96 + t145 + t157, g(3) * t120 + t76 * t96 + t156 * t122 + (t33 + t180) * qJD(4) + t166, 0, 0, t186, t201, t202, -t186, t203, t125, -t11 * t131 + (t192 * t125 - t131 * t171 + t40 * t78) * pkin(4) + t199, t12 * t131 + (-t125 * t136 - t131 * t160 + t148 * t78) * pkin(4) + t200, -t10 * t148 - t11 * t148 - t12 * t40 + t9 * t40 + (-t192 * t147 + t136 * t142 + (-t136 * t148 + t192 * t40) * qJD(5)) * pkin(4), -t10 * t12 - t9 * t11 + (t192 * t2 + t1 * t136 - t52 * t78 + (t192 * t10 - t136 * t9) * qJD(5) + t145) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t201, t202, -t186, t203, t125, t10 * t131 + t199, t9 * t131 + t200, 0, 0;];
tau_reg = t5;
