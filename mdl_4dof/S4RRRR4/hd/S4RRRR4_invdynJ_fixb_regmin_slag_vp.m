% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRR4
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:17
% EndTime: 2019-12-31 17:26:20
% DurationCPUTime: 1.57s
% Computational Cost: add. (1612->220), mult. (3718->316), div. (0->0), fcn. (2676->10), ass. (0->136)
t103 = cos(qJ(2));
t174 = cos(qJ(3));
t139 = qJD(1) * t174;
t100 = sin(qJ(2));
t152 = qJD(1) * t100;
t99 = sin(qJ(3));
t183 = t103 * t139 - t99 * t152;
t94 = qJD(2) + qJD(3);
t185 = t183 * t94;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t127 = g(1) * t104 + g(2) * t101;
t97 = qJ(2) + qJ(3);
t91 = sin(t97);
t184 = t127 * t91;
t92 = cos(t97);
t175 = g(3) * t92;
t176 = pkin(5) + pkin(6);
t67 = t176 * t103;
t61 = qJD(1) * t67;
t144 = t174 * t61;
t160 = qJD(2) * pkin(2);
t66 = t176 * t100;
t59 = qJD(1) * t66;
t56 = -t59 + t160;
t33 = t99 * t56 + t144;
t150 = qJD(1) * qJD(2);
t136 = t103 * t150;
t149 = t100 * qJDD(1);
t39 = qJDD(2) * pkin(2) + t176 * (-t136 - t149);
t137 = t100 * t150;
t148 = t103 * qJDD(1);
t40 = t176 * (-t137 + t148);
t178 = t33 * qJD(3) - t174 * t39 + t99 * t40;
t93 = qJDD(2) + qJDD(3);
t9 = -t93 * pkin(3) + t178;
t181 = t9 + t175;
t157 = t99 * t103;
t55 = -qJD(1) * t157 - t100 * t139;
t30 = -t55 * pkin(3) - pkin(7) * t183;
t50 = qJD(4) - t183;
t88 = t99 * pkin(2) + pkin(7);
t180 = (pkin(2) * t152 + qJD(4) * t88 + t30) * t50;
t179 = (pkin(7) * qJD(4) + t30) * t50;
t102 = cos(qJ(4));
t98 = sin(qJ(4));
t124 = t102 * t55 - t98 * t94;
t132 = qJDD(1) * t174;
t21 = t100 * t132 + t99 * t148 + t185;
t11 = -t124 * qJD(4) - t102 * t93 + t98 * t21;
t118 = -t99 * t100 + t174 * t103;
t138 = t174 * qJD(3);
t156 = qJD(3) * t99;
t110 = t56 * t138 - t61 * t156 + t174 * t40 + t99 * t39;
t90 = -t103 * pkin(2) - pkin(1);
t65 = t90 * qJD(1);
t26 = -pkin(3) * t183 + t55 * pkin(7) + t65;
t134 = t93 * pkin(7) + qJD(4) * t26 + t110;
t120 = -t174 * t66 - t99 * t67;
t142 = qJD(2) * t176;
t60 = t100 * t142;
t62 = t103 * t142;
t17 = t120 * qJD(3) - t174 * t60 - t99 * t62;
t125 = -t103 * t132 + t99 * t149;
t58 = t174 * t100 + t157;
t38 = t94 * t58;
t22 = qJD(1) * t38 + t125;
t20 = qJDD(4) + t22;
t162 = t99 * t61;
t32 = t174 * t56 - t162;
t28 = -t94 * pkin(3) - t32;
t31 = -pkin(3) * t118 - t58 * pkin(7) + t90;
t37 = t94 * t118;
t45 = t174 * t67 - t99 * t66;
t177 = -(qJD(4) * t31 + t17) * t50 + t134 * t118 - t45 * t20 + t28 * t37 + t9 * t58;
t86 = g(3) * t91;
t151 = qJD(4) * t102;
t155 = qJD(4) * t98;
t10 = t102 * t21 + t94 * t151 + t55 * t155 + t98 * t93;
t171 = t10 * t98;
t170 = t28 * t183;
t169 = t28 * t58;
t168 = t31 * t20;
t41 = -t102 * t94 - t98 * t55;
t167 = t41 * t50;
t166 = t124 * t50;
t165 = t50 * t55;
t164 = t55 * t183;
t163 = t98 * t20;
t95 = t100 ^ 2;
t161 = -t103 ^ 2 + t95;
t159 = t101 * t98;
t158 = t104 * t98;
t154 = t101 * t102;
t153 = t104 * t102;
t147 = t100 * t160;
t146 = t58 * t155;
t29 = t94 * pkin(7) + t33;
t51 = pkin(2) * t137 + t90 * qJDD(1);
t6 = t22 * pkin(3) - t21 * pkin(7) + t51;
t133 = qJD(4) * t29 - t6;
t129 = t102 * t50;
t35 = -t99 * t59 + t144;
t128 = pkin(2) * t156 - t35;
t126 = g(1) * t101 - g(2) * t104;
t13 = t102 * t29 + t98 * t26;
t123 = -t13 * t55 + t28 * t151 + t181 * t98;
t12 = t102 * t26 - t98 * t29;
t122 = t102 * t184 + t12 * t55 + t28 * t155;
t121 = -t134 + t86;
t117 = -0.2e1 * pkin(1) * t150 - pkin(5) * qJDD(2);
t114 = -pkin(7) * t20 + t32 * t50 - t170;
t105 = qJD(2) ^ 2;
t112 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t105 + t126;
t106 = qJD(1) ^ 2;
t111 = pkin(1) * t106 - pkin(5) * qJDD(1) + t127;
t36 = -t174 * t59 - t162;
t109 = -t88 * t20 - t170 + (-pkin(2) * t138 + t36) * t50;
t108 = t127 * t92 - t183 * t65 - t110 + t86;
t107 = t65 * t55 - t175 - t178 + t184;
t89 = -t174 * pkin(2) - pkin(3);
t49 = t92 * t153 + t159;
t48 = -t92 * t158 + t154;
t47 = -t92 * t154 + t158;
t46 = t92 * t159 + t153;
t23 = -t183 ^ 2 + t55 ^ 2;
t18 = t45 * qJD(3) + t174 * t62 - t99 * t60;
t16 = t38 * pkin(3) - t37 * pkin(7) + t147;
t15 = -t125 + (-qJD(1) * t58 - t55) * t94;
t14 = t21 - t185;
t5 = t102 * t6;
t4 = -t124 * t55 + t50 * t129 + t163;
t3 = -t50 ^ 2 * t98 + t102 * t20 - t41 * t55;
t2 = -t124 * t129 + t171;
t1 = (-t11 + t166) * t98 + (t10 - t167) * t102;
t7 = [qJDD(1), t126, t127, t95 * qJDD(1) + 0.2e1 * t100 * t136, 0.2e1 * t100 * t148 - 0.2e1 * t161 * t150, qJDD(2) * t100 + t105 * t103, qJDD(2) * t103 - t105 * t100, 0, t117 * t100 + t112 * t103, -t112 * t100 + t117 * t103, t21 * t58 - t55 * t37, t118 * t21 + t183 * t37 - t58 * t22 + t55 * t38, t37 * t94 + t58 * t93, t118 * t93 - t38 * t94, 0, -t118 * t51 + t120 * t93 + t126 * t92 - t147 * t183 - t18 * t94 + t90 * t22 + t65 * t38, -t126 * t91 - t55 * t147 - t17 * t94 + t90 * t21 + t65 * t37 - t45 * t93 + t51 * t58, t124 * t146 + (t10 * t58 - t124 * t37) * t102, (-t102 * t41 + t124 * t98) * t37 + (-t171 - t102 * t11 + (t102 * t124 + t41 * t98) * qJD(4)) * t58, -t50 * t146 - t10 * t118 - t124 * t38 + (t20 * t58 + t37 * t50) * t102, -t58 * t163 + t11 * t118 - t41 * t38 + (-t58 * t151 - t98 * t37) * t50, -t118 * t20 + t50 * t38, -g(1) * t47 - g(2) * t49 - t120 * t11 + t12 * t38 + t18 * t41 - t5 * t118 + (t16 * t50 + t168 + (t118 * t29 - t45 * t50 + t169) * qJD(4)) * t102 + t177 * t98, -g(1) * t46 - g(2) * t48 - t120 * t10 - t13 * t38 - t18 * t124 + (-(-qJD(4) * t45 + t16) * t50 - t168 - t133 * t118 - qJD(4) * t169) * t98 + t177 * t102; 0, 0, 0, -t100 * t106 * t103, t161 * t106, t149, t148, qJDD(2), -g(3) * t103 + t111 * t100, g(3) * t100 + t111 * t103, t164, t23, t14, t15, t93, t35 * t94 + (t152 * t183 - t94 * t156 + t174 * t93) * pkin(2) + t107, t36 * t94 + (-t94 * t138 + t55 * t152 - t93 * t99) * pkin(2) + t108, t2, t1, t4, t3, t165, t89 * t11 + t128 * t41 + t109 * t98 + (-t181 - t180) * t102 + t122, t89 * t10 - t128 * t124 + (-t184 + t180) * t98 + t109 * t102 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t23, t14, t15, t93, t33 * t94 + t107, t32 * t94 + t108, t2, t1, t4, t3, t165, -pkin(3) * t11 - t33 * t41 + t114 * t98 + (-t181 - t179) * t102 + t122, -pkin(3) * t10 + t33 * t124 + t114 * t102 + (-t184 + t179) * t98 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 * t41, t124 ^ 2 - t41 ^ 2, t10 + t167, -t11 - t166, t20, -g(1) * t48 + g(2) * t46 + t121 * t98 + t124 * t28 + t13 * t50 - t29 * t151 + t5, g(1) * t49 - g(2) * t47 + t121 * t102 + t12 * t50 + t133 * t98 + t28 * t41;];
tau_reg = t7;
