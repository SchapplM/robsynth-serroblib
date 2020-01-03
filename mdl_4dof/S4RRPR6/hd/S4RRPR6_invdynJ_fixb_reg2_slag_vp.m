% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:57
% EndTime: 2019-12-31 17:05:00
% DurationCPUTime: 1.54s
% Computational Cost: add. (2194->255), mult. (5387->345), div. (0->0), fcn. (3875->12), ass. (0->139)
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t128 = sin(pkin(7));
t129 = cos(pkin(7));
t135 = cos(qJ(2));
t173 = t129 * t135;
t162 = qJD(1) * t173;
t132 = sin(qJ(2));
t170 = qJD(1) * t132;
t80 = t128 * t170 - t162;
t92 = t128 * t135 + t129 * t132;
t83 = t92 * qJD(1);
t149 = t131 * t80 - t134 * t83;
t165 = t135 * qJDD(1);
t166 = t132 * qJDD(1);
t150 = t128 * t166 - t129 * t165;
t82 = t92 * qJD(2);
t46 = qJD(1) * t82 + t150;
t167 = qJD(1) * qJD(2);
t160 = t132 * t167;
t146 = qJDD(1) * t92 - t128 * t160;
t159 = t135 * t167;
t47 = t129 * t159 + t146;
t143 = qJD(4) * t149 - t131 * t47 - t134 * t46;
t124 = qJD(2) + qJD(4);
t176 = t149 * t124;
t202 = t143 - t176;
t168 = qJD(4) * t134;
t169 = qJD(4) * t131;
t147 = -t131 * t46 + t134 * t47 - t80 * t168 - t83 * t169;
t37 = -t131 * t83 - t134 * t80;
t175 = t37 * t124;
t201 = t147 - t175;
t184 = t149 ^ 2;
t185 = t37 ^ 2;
t200 = t184 - t185;
t183 = t37 * t149;
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t152 = g(1) * t136 + g(2) * t133;
t181 = qJ(3) + pkin(5);
t105 = t181 * t132;
t155 = qJD(2) * t181;
t77 = -t132 * qJD(3) - t135 * t155;
t43 = qJDD(2) * pkin(2) + t77 * qJD(1) - qJDD(1) * t105;
t106 = t181 * t135;
t76 = t135 * qJD(3) - t132 * t155;
t50 = t76 * qJD(1) + qJDD(1) * t106;
t17 = t128 * t43 + t129 * t50;
t11 = -t46 * pkin(6) + t17;
t193 = t83 * pkin(6);
t98 = qJD(1) * t106;
t86 = t128 * t98;
t178 = qJD(2) * pkin(2);
t97 = qJD(1) * t105;
t90 = -t97 + t178;
t44 = t129 * t90 - t86;
t23 = qJD(2) * pkin(3) - t193 + t44;
t194 = t80 * pkin(6);
t177 = t129 * t98;
t45 = t128 * t90 + t177;
t25 = t45 - t194;
t16 = -t128 * t50 + t129 * t43;
t8 = qJDD(2) * pkin(3) - t47 * pkin(6) + t16;
t1 = (qJD(4) * t23 + t11) * t134 + t131 * t8 - t25 * t169;
t125 = qJ(2) + pkin(7);
t120 = qJ(4) + t125;
t114 = sin(t120);
t115 = cos(t120);
t186 = t135 * pkin(2);
t117 = pkin(1) + t186;
t100 = -t117 * qJD(1) + qJD(3);
t53 = t80 * pkin(3) + t100;
t199 = g(3) * t114 + t152 * t115 - t53 * t37 - t1;
t10 = t131 * t23 + t134 * t25;
t2 = -qJD(4) * t10 - t131 * t11 + t134 * t8;
t198 = -g(3) * t115 + t152 * t114 + t149 * t53 + t2;
t197 = g(1) * t133 - g(2) * t136;
t196 = t83 ^ 2;
t195 = t46 * pkin(3);
t192 = pkin(2) * t128;
t188 = g(3) * t135;
t187 = t132 * pkin(2);
t182 = t83 * t80;
t29 = t128 * t77 + t129 * t76;
t51 = t128 * t97 - t177;
t26 = t51 + t194;
t52 = -t129 * t97 - t86;
t27 = t52 - t193;
t116 = t129 * pkin(2) + pkin(3);
t74 = t134 * t116 - t131 * t192;
t180 = t74 * qJD(4) - t131 * t26 - t134 * t27;
t75 = t131 * t116 + t134 * t192;
t179 = -t75 * qJD(4) + t131 * t27 - t134 * t26;
t55 = -t128 * t105 + t129 * t106;
t174 = pkin(5) * qJDD(1);
t126 = t132 ^ 2;
t127 = t135 ^ 2;
t172 = t126 - t127;
t171 = t126 + t127;
t164 = t132 * t178;
t138 = qJD(1) ^ 2;
t163 = t132 * t138 * t135;
t119 = cos(t125);
t158 = pkin(3) * t119 + t186;
t28 = -t128 * t76 + t129 * t77;
t54 = -t129 * t105 - t128 * t106;
t153 = t132 * t159;
t30 = -t92 * pkin(6) + t54;
t91 = t128 * t132 - t173;
t31 = -t91 * pkin(6) + t55;
t14 = -t131 * t31 + t134 * t30;
t15 = t131 * t30 + t134 * t31;
t49 = -t131 * t91 + t134 * t92;
t148 = -0.2e1 * pkin(1) * t167 - pkin(5) * qJDD(2);
t72 = pkin(2) * t160 - t117 * qJDD(1) + qJDD(3);
t137 = qJD(2) ^ 2;
t145 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t137 + t197;
t144 = pkin(1) * t138 + t152 - t174;
t142 = -t197 + t72;
t123 = -pkin(6) - t181;
t122 = qJDD(2) + qJDD(4);
t118 = sin(t125);
t96 = pkin(1) + t158;
t85 = t91 * qJD(2);
t78 = t80 ^ 2;
t59 = t91 * pkin(3) - t117;
t57 = t82 * pkin(3) + t164;
t56 = pkin(2) * t170 + t83 * pkin(3);
t48 = t131 * t92 + t134 * t91;
t24 = t72 + t195;
t21 = -t82 * pkin(6) + t29;
t20 = t85 * pkin(6) + t28;
t19 = qJD(4) * t49 - t131 * t85 + t134 * t82;
t18 = t131 * t82 + t134 * t85 + t91 * t168 + t92 * t169;
t9 = -t131 * t25 + t134 * t23;
t4 = -qJD(4) * t15 - t131 * t21 + t134 * t20;
t3 = qJD(4) * t14 + t131 * t20 + t134 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), t197, t152, 0, 0, t126 * qJDD(1) + 0.2e1 * t153, 0.2e1 * t132 * t165 - 0.2e1 * t172 * t167, qJDD(2) * t132 + t137 * t135, t127 * qJDD(1) - 0.2e1 * t153, qJDD(2) * t135 - t137 * t132, 0, t132 * t148 + t135 * t145, -t132 * t145 + t135 * t148, 0.2e1 * t171 * t174 - t152, -g(1) * (-t133 * pkin(1) + t136 * pkin(5)) - g(2) * (t136 * pkin(1) + t133 * pkin(5)) + (t171 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t47 * t92 - t83 * t85, -t92 * t46 - t47 * t91 + t85 * t80 - t83 * t82, -t85 * qJD(2) + t92 * qJDD(2), t46 * t91 + t80 * t82, -t82 * qJD(2) - t91 * qJDD(2), 0, t54 * qJDD(2) + t100 * t82 - t117 * t46 + t72 * t91 + t197 * t119 + (t80 * t187 + t28) * qJD(2), -t55 * qJDD(2) - t100 * t85 - t117 * t47 + t72 * t92 - t197 * t118 + (t83 * t187 - t29) * qJD(2), -t16 * t92 - t17 * t91 - t28 * t83 - t29 * t80 + t44 * t85 - t45 * t82 - t55 * t46 - t54 * t47 - t152, t17 * t55 + t45 * t29 + t16 * t54 + t44 * t28 - t72 * t117 + t100 * t164 - g(1) * (-t133 * t117 + t136 * t181) - g(2) * (t136 * t117 + t133 * t181), t147 * t49 + t149 * t18, t143 * t49 - t147 * t48 + t149 * t19 - t18 * t37, t49 * t122 - t18 * t124, -t143 * t48 - t19 * t37, -t48 * t122 - t19 * t124, 0, t115 * t197 + t14 * t122 + t4 * t124 - t143 * t59 + t53 * t19 + t24 * t48 - t37 * t57, -t114 * t197 - t15 * t122 - t3 * t124 + t147 * t59 - t149 * t57 - t53 * t18 + t24 * t49, -t1 * t48 - t10 * t19 - t14 * t147 + t143 * t15 + t149 * t4 + t9 * t18 - t2 * t49 + t3 * t37 - t152, t1 * t15 + t10 * t3 + t2 * t14 + t9 * t4 + t24 * t59 + t53 * t57 - g(1) * (-t136 * t123 - t133 * t96) - g(2) * (-t133 * t123 + t136 * t96); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t172 * t138, t166, t163, t165, qJDD(2), t132 * t144 - t188, g(3) * t132 + t135 * t144, 0, 0, t182, -t78 + t196, (t80 + t162) * qJD(2) + t146, -t182, -t150, qJDD(2), -g(3) * t119 - t51 * qJD(2) - t100 * t83 + t152 * t118 + (qJDD(2) * t129 - t80 * t170) * pkin(2) + t16, g(3) * t118 + t52 * qJD(2) + t100 * t80 + t152 * t119 + (-qJDD(2) * t128 - t83 * t170) * pkin(2) - t17, (t45 + t51) * t83 + (-t44 + t52) * t80 + (-t128 * t46 - t129 * t47) * pkin(2), -t44 * t51 - t45 * t52 + (-t188 + t128 * t17 + t129 * t16 + (-qJD(1) * t100 + t152) * t132) * pkin(2), t183, t200, t201, -t183, t202, t122, t74 * t122 + t179 * t124 + t37 * t56 + t198, -t75 * t122 - t124 * t180 + t149 * t56 + t199, -t147 * t74 + t75 * t143 + (t180 + t9) * t37 + (-t10 + t179) * t149, t1 * t75 + t2 * t74 - t53 * t56 - g(3) * t158 - t152 * (-pkin(3) * t118 - t187) + t179 * t9 + t180 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t83 * qJD(2) + t150, (-t80 + t162) * qJD(2) + t146, -t78 - t196, t44 * t83 + t45 * t80 + t142, 0, 0, 0, 0, 0, 0, -t143 - t176, t147 + t175, -t184 - t185, -t10 * t37 - t149 * t9 + t142 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t200, t201, -t183, t202, t122, t10 * t124 + t198, t9 * t124 + t199, 0, 0;];
tau_reg = t5;
