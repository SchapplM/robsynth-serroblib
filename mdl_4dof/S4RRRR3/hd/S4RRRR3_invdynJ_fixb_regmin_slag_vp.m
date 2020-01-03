% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR3
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:38
% EndTime: 2019-12-31 17:24:42
% DurationCPUTime: 1.33s
% Computational Cost: add. (1456->200), mult. (3443->288), div. (0->0), fcn. (2521->12), ass. (0->129)
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t117 = cos(qJ(3));
t118 = cos(qJ(2));
t166 = qJD(1) * t118;
t154 = t117 * t166;
t113 = sin(qJ(3));
t114 = sin(qJ(2));
t167 = qJD(1) * t114;
t155 = t113 * t167;
t55 = -t154 + t155;
t57 = -t113 * t166 - t117 * t167;
t137 = t112 * t55 + t116 * t57;
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t143 = g(1) * t119 + g(2) * t115;
t107 = qJDD(2) + qJDD(3);
t183 = pkin(5) + pkin(6);
t81 = t183 * t118;
t73 = qJD(1) * t81;
t62 = t117 * t73;
t173 = qJD(2) * pkin(2);
t80 = t183 * t114;
t71 = qJD(1) * t80;
t64 = -t71 + t173;
t136 = -t113 * t64 - t62;
t161 = qJD(1) * qJD(2);
t152 = t118 * t161;
t160 = t114 * qJDD(1);
t185 = t152 + t160;
t42 = qJDD(2) * pkin(2) - t183 * t185;
t153 = t114 * t161;
t159 = t118 * qJDD(1);
t43 = t183 * (-t153 + t159);
t127 = t136 * qJD(3) - t113 * t43 + t117 * t42;
t108 = qJD(2) + qJD(3);
t23 = qJD(3) * t154 - t108 * t155 + t113 * t159 + t185 * t117;
t4 = t107 * pkin(3) - t23 * pkin(7) + t127;
t101 = -t118 * pkin(2) - pkin(1);
t79 = t101 * qJD(1);
t44 = t55 * pkin(3) + t79;
t165 = qJD(3) * t113;
t184 = (qJD(3) * t64 + t43) * t117 + t113 * t42 - t73 * t165;
t141 = t113 * t160 - t117 * t159;
t67 = t113 * t118 + t117 * t114;
t41 = t108 * t67;
t24 = qJD(1) * t41 + t141;
t6 = -t24 * pkin(7) + t184;
t111 = qJ(2) + qJ(3);
t106 = qJ(4) + t111;
t98 = sin(t106);
t99 = cos(t106);
t189 = -g(3) * t99 - t112 * t6 + t116 * t4 + t44 * t137 + t143 * t98;
t163 = qJD(4) * t112;
t182 = t55 * pkin(7);
t22 = -t136 - t182;
t33 = t112 * t57 - t116 * t55;
t188 = g(3) * t98 + t143 * t99 + t22 * t163 - t44 * t33;
t58 = t113 * t73;
t148 = t117 * t64 - t58;
t51 = t57 * pkin(7);
t21 = t148 + t51;
t16 = t108 * pkin(3) + t21;
t172 = t116 * t22;
t140 = -t112 * t16 - t172;
t187 = t140 * qJD(4) + t189;
t158 = -qJD(3) - qJD(4);
t103 = qJD(2) - t158;
t186 = (-t22 * t103 - t4) * t112 + t188;
t179 = t137 * t33;
t5 = t137 ^ 2 - t33 ^ 2;
t162 = qJD(4) * t116;
t7 = -t112 * t24 + t116 * t23 - t55 * t162 + t57 * t163;
t1 = -t33 * t103 + t7;
t125 = t137 * qJD(4) - t112 * t23 - t116 * t24;
t2 = -t103 * t137 + t125;
t176 = t57 * t55;
t175 = -t117 * t71 - t58;
t174 = -t113 * t80 + t117 * t81;
t102 = qJDD(4) + t107;
t171 = t113 * t102;
t170 = t113 * t116;
t169 = t116 * t102;
t109 = t114 ^ 2;
t168 = -t118 ^ 2 + t109;
t164 = qJD(3) * t117;
t157 = t114 * t173;
t156 = qJD(2) * t183;
t151 = -qJD(4) * t16 - t6;
t147 = t113 * t71 - t62;
t146 = -t113 * t81 - t117 * t80;
t142 = g(1) * t115 - g(2) * t119;
t28 = -t67 * pkin(7) + t146;
t66 = t113 * t114 - t117 * t118;
t29 = -t66 * pkin(7) + t174;
t139 = -t112 * t29 + t116 * t28;
t138 = t112 * t28 + t116 * t29;
t38 = t112 * t67 + t116 * t66;
t39 = -t112 * t66 + t116 * t67;
t135 = -0.2e1 * pkin(1) * t161 - pkin(5) * qJDD(2);
t72 = t114 * t156;
t74 = t118 * t156;
t134 = -t113 * t74 - t117 * t72 - t80 * t164 - t81 * t165;
t52 = pkin(2) * t153 + t101 * qJDD(1);
t120 = qJD(2) ^ 2;
t129 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t120 + t142;
t121 = qJD(1) ^ 2;
t128 = pkin(1) * t121 - pkin(5) * qJDD(1) + t143;
t126 = -t174 * qJD(3) + t113 * t72 - t117 * t74;
t104 = sin(t111);
t105 = cos(t111);
t124 = g(3) * t104 + t143 * t105 + t79 * t55 - t184;
t122 = -g(3) * t105 + t143 * t104 + t79 * t57 + t127;
t100 = t117 * pkin(2) + pkin(3);
t47 = t66 * pkin(3) + t101;
t45 = pkin(2) * t167 - t57 * pkin(3);
t40 = t108 * t66;
t35 = t41 * pkin(3) + t157;
t27 = t51 + t175;
t26 = t147 + t182;
t25 = -t55 ^ 2 + t57 ^ 2;
t15 = t24 * pkin(3) + t52;
t14 = -t141 + (-qJD(1) * t67 - t57) * t108;
t13 = t55 * t108 + t23;
t12 = t40 * pkin(7) + t126;
t11 = -t41 * pkin(7) + t134;
t10 = t39 * qJD(4) - t112 * t40 + t116 * t41;
t9 = -t38 * qJD(4) - t112 * t41 - t116 * t40;
t3 = [qJDD(1), t142, t143, t109 * qJDD(1) + 0.2e1 * t114 * t152, 0.2e1 * t114 * t159 - 0.2e1 * t168 * t161, qJDD(2) * t114 + t120 * t118, qJDD(2) * t118 - t120 * t114, 0, t135 * t114 + t129 * t118, -t129 * t114 + t135 * t118, t23 * t67 + t57 * t40, -t23 * t66 - t67 * t24 + t40 * t55 + t57 * t41, t67 * t107 - t40 * t108, -t66 * t107 - t41 * t108, 0, t101 * t24 + t142 * t105 + t146 * t107 + t126 * t108 + t55 * t157 + t79 * t41 + t52 * t66, t101 * t23 - t142 * t104 - t174 * t107 - t134 * t108 - t57 * t157 - t79 * t40 + t52 * t67, -t137 * t9 + t7 * t39, t10 * t137 + t125 * t39 + t33 * t9 - t7 * t38, t39 * t102 + t9 * t103, -t10 * t103 - t38 * t102, 0, -t35 * t33 - t47 * t125 + t15 * t38 + t44 * t10 + (-t138 * qJD(4) - t112 * t11 + t116 * t12) * t103 + t139 * t102 + t142 * t99, -t35 * t137 + t47 * t7 + t15 * t39 + t44 * t9 - (t139 * qJD(4) + t116 * t11 + t112 * t12) * t103 - t138 * t102 - t142 * t98; 0, 0, 0, -t114 * t121 * t118, t168 * t121, t160, t159, qJDD(2), -g(3) * t118 + t114 * t128, g(3) * t114 + t118 * t128, -t176, t25, t13, t14, t107, -t147 * t108 + (t117 * t107 - t108 * t165 - t55 * t167) * pkin(2) + t122, t175 * t108 + (-t113 * t107 - t108 * t164 + t57 * t167) * pkin(2) + t124, t179, t5, t1, t2, t102, t100 * t169 + t45 * t33 - (-t112 * t27 + t116 * t26) * t103 + (-t112 * t171 + (-t112 * t117 - t170) * t103 * qJD(3)) * pkin(2) + ((-pkin(2) * t170 - t112 * t100) * t103 + t140) * qJD(4) + t189, t45 * t137 + (-t100 * t102 - t4 + (-pkin(2) * t113 * t158 + t26) * t103) * t112 + (-pkin(2) * t171 + (-pkin(2) * t164 - qJD(4) * t100 + t27) * t103 + t151) * t116 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t25, t13, t14, t107, -t136 * t108 + t122, t148 * t108 + t124, t179, t5, t1, t2, t102, -(-t112 * t21 - t172) * t103 + (-t103 * t163 - t33 * t57 + t169) * pkin(3) + t187, (t21 * t103 + t151) * t116 + (-t112 * t102 - t103 * t162 - t137 * t57) * pkin(3) + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t5, t1, t2, t102, -t140 * t103 + t187, (-t6 + (-qJD(4) + t103) * t16) * t116 + t186;];
tau_reg = t3;
