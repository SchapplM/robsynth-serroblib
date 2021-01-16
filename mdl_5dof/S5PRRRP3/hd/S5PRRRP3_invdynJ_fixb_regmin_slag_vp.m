% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:32
% DurationCPUTime: 1.09s
% Computational Cost: add. (1532->217), mult. (3295->270), div. (0->0), fcn. (2250->8), ass. (0->137)
t170 = cos(qJ(4));
t101 = sin(qJ(4));
t102 = sin(qJ(3));
t147 = t102 * qJD(1);
t103 = cos(qJ(3));
t176 = pkin(6) + pkin(7);
t69 = t176 * t103;
t39 = qJD(2) * t69 + t147;
t33 = t101 * t39;
t158 = qJD(3) * pkin(3);
t140 = qJD(2) * t176;
t38 = t103 * qJD(1) - t102 * t140;
t36 = t38 + t158;
t133 = t170 * t36 - t33;
t51 = t101 * t103 + t170 * t102;
t44 = t51 * qJD(2);
t153 = t44 * qJ(5);
t10 = t133 - t153;
t97 = qJD(3) + qJD(4);
t136 = t170 * qJD(4);
t175 = pkin(3) * t97;
t94 = qJDD(3) + qJDD(4);
t181 = -t101 * pkin(3) * t94 - t136 * t175;
t152 = t101 * t102;
t120 = t97 * t152;
t138 = t170 * t103;
t125 = qJD(2) * t138;
t130 = qJDD(2) * t170;
t144 = t103 * qJDD(2);
t129 = -t101 * t144 - t102 * t130 - t97 * t125;
t16 = qJD(2) * t120 + t129;
t156 = t16 * qJ(5);
t86 = t94 * pkin(4);
t180 = t156 + t86;
t68 = t176 * t102;
t89 = t103 * qJDD(1);
t21 = qJDD(3) * pkin(3) + t89 - qJDD(2) * t68 + (-t103 * t140 - t147) * qJD(3);
t25 = t38 * qJD(3) + t102 * qJDD(1) + qJDD(2) * t69;
t179 = -t101 * t25 + t170 * t21;
t96 = pkin(8) + qJ(2);
t88 = cos(t96);
t100 = qJ(3) + qJ(4);
t91 = sin(t100);
t164 = t88 * t91;
t87 = sin(t96);
t166 = t87 * t91;
t92 = cos(t100);
t174 = g(3) * t92;
t178 = g(1) * t164 + g(2) * t166 - t174;
t177 = t44 ^ 2;
t93 = t103 * pkin(3);
t173 = pkin(2) + t93;
t8 = t97 * pkin(4) + t10;
t172 = t10 - t8;
t145 = t102 * qJDD(2);
t121 = t101 * t145 - t103 * t130;
t27 = t97 * t51;
t17 = t27 * qJD(2) + t121;
t26 = -qJD(3) * t138 - t103 * t136 + t120;
t149 = qJD(2) * t102;
t137 = t101 * t149;
t42 = -t125 + t137;
t171 = -t51 * t17 + t26 * t42;
t169 = t42 * t97;
t168 = t44 * t42;
t165 = t87 * t92;
t163 = t88 * t92;
t162 = t170 * t38 - t33;
t161 = -t101 * t68 + t170 * t69;
t160 = pkin(4) * t92 + t93;
t98 = t102 ^ 2;
t159 = -t103 ^ 2 + t98;
t155 = t17 * qJ(5);
t154 = t42 * qJ(5);
t134 = t42 * pkin(4) + qJD(5);
t67 = t173 * qJD(2);
t28 = t134 - t67;
t151 = qJD(5) + t28;
t150 = qJDD(1) - g(3);
t148 = qJD(4) * t101;
t146 = qJD(2) * qJD(3);
t143 = t102 * t158;
t142 = pkin(3) * t149;
t35 = t170 * t39;
t139 = qJD(3) * t176;
t135 = t102 * t146;
t132 = -t101 * t38 - t35;
t131 = -t101 * t69 - t170 * t68;
t128 = -g(1) * t166 + g(2) * t164;
t127 = g(1) * t165 - g(2) * t163;
t124 = g(1) * t88 + g(2) * t87;
t123 = g(1) * t87 - g(2) * t88;
t50 = -t138 + t152;
t122 = -t50 * t16 + t44 * t27;
t14 = t26 * t97 - t51 * t94;
t119 = -t101 * t36 - t35;
t118 = -0.2e1 * pkin(2) * t146 - pkin(6) * qJDD(3);
t56 = t102 * t139;
t57 = t103 * t139;
t117 = -t101 * t57 - t68 * t136 - t69 * t148 - t170 * t56;
t40 = pkin(3) * t135 - qJDD(2) * t173;
t104 = qJD(3) ^ 2;
t116 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t104 + t123;
t105 = qJD(2) ^ 2;
t115 = t105 * pkin(2) - qJDD(2) * pkin(6) + t124;
t5 = t17 * pkin(4) + qJDD(5) + t40;
t114 = -t97 * t137 - t129;
t113 = t119 * qJD(4) + t179;
t112 = -t161 * qJD(4) + t101 * t56 - t170 * t57;
t111 = t101 * t21 + t36 * t136 - t39 * t148 + t170 * t25;
t110 = g(1) * t163 + g(2) * t165 + g(3) * t91 - t111;
t109 = t113 + t178;
t108 = -t67 * t42 + t110;
t107 = t67 * t44 + t109;
t106 = t151 * t42 + t110 + t155;
t95 = -qJ(5) - t176;
t84 = t170 * pkin(3) + pkin(4);
t66 = qJDD(3) * t103 - t104 * t102;
t65 = qJDD(3) * t102 + t104 * t103;
t55 = pkin(2) + t160;
t41 = t42 ^ 2;
t32 = t50 * pkin(4) - t173;
t29 = t44 * pkin(4) + t142;
t24 = t27 * pkin(4) + t143;
t23 = -t50 * qJ(5) + t161;
t22 = -t51 * qJ(5) + t131;
t18 = -t41 + t177;
t15 = -t27 * t97 - t50 * t94;
t13 = -t153 + t162;
t12 = t132 + t154;
t11 = -t119 - t154;
t6 = t114 + t169;
t4 = t26 * qJ(5) - t51 * qJD(5) + t112;
t3 = -t27 * qJ(5) - t50 * qJD(5) + t117;
t2 = -t42 * qJD(5) + t111 - t155;
t1 = -t44 * qJD(5) + t113 + t180;
t7 = [t150, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t65, 0, 0, 0, 0, 0, t15, t14, t15, t14, t122 + t171, -t1 * t50 - t11 * t26 + t2 * t51 - t8 * t27 - g(3); 0, qJDD(2), t123, t124, t98 * qJDD(2) + 0.2e1 * t103 * t135, 0.2e1 * t102 * t144 - 0.2e1 * t159 * t146, t65, t66, 0, t118 * t102 + t116 * t103, -t116 * t102 + t118 * t103, -t16 * t51 - t44 * t26, -t122 + t171, -t14, t15, 0, t112 * t97 + t131 * t94 + t42 * t143 - t17 * t173 - t67 * t27 + t40 * t50 + t127, -t117 * t97 + t44 * t143 + t16 * t173 - t161 * t94 + t67 * t26 + t40 * t51 + t128, t32 * t17 + t22 * t94 + t24 * t42 + t28 * t27 + t4 * t97 + t5 * t50 + t127, -t32 * t16 - t23 * t94 + t24 * t44 - t28 * t26 - t3 * t97 + t5 * t51 + t128, -t1 * t51 - t11 * t27 + t22 * t16 - t23 * t17 - t2 * t50 + t8 * t26 - t3 * t42 - t4 * t44 - t124, t2 * t23 + t11 * t3 + t1 * t22 + t8 * t4 + t5 * t32 + t28 * t24 - g(1) * (-t87 * t55 - t88 * t95) - g(2) * (t88 * t55 - t87 * t95); 0, 0, 0, 0, -t102 * t105 * t103, t159 * t105, t145, t144, qJDD(3), -g(3) * t103 + t115 * t102 + t89, -t150 * t102 + t115 * t103, t168, t18, t6, -t121, t94, -t132 * t97 + (-t97 * t148 - t42 * t149 + t170 * t94) * pkin(3) + t107, -t44 * t142 + t162 * t97 + t108 + t181, -t12 * t97 - t29 * t42 + t84 * t94 - t151 * t44 + (-t35 + (-t36 - t175) * t101) * qJD(4) + t178 + t179 + t180, t13 * t97 - t29 * t44 + t106 + t181, t84 * t16 + (t11 + t12) * t44 + (t13 - t8) * t42 + (-t101 * t17 + (t101 * t44 - t170 * t42) * qJD(4)) * pkin(3), t1 * t84 - t11 * t13 - t8 * t12 - t28 * t29 - g(3) * t160 - t124 * (-t102 * pkin(3) - pkin(4) * t91) + (t2 * t101 + (-t101 * t8 + t170 * t11) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t18, t6, -t121, t94, -t119 * t97 + t107, t133 * t97 + t108, t156 + t11 * t97 + 0.2e1 * t86 + (-t134 - t28) * t44 + t109, -t177 * pkin(4) + t10 * t97 + t106, t16 * pkin(4) + t172 * t42, -t172 * t11 + (t124 * t91 - t28 * t44 + t1 - t174) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t97 + t17, t114 - t169, -t41 - t177, t11 * t42 + t8 * t44 - t123 + t5;];
tau_reg = t7;
