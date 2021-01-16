% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:41
% EndTime: 2021-01-15 15:52:48
% DurationCPUTime: 1.14s
% Computational Cost: add. (1118->185), mult. (2979->285), div. (0->0), fcn. (2206->8), ass. (0->115)
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t106 = cos(qJ(3));
t142 = t101 * t106;
t127 = qJD(2) * t142;
t103 = sin(qJ(3));
t138 = qJD(2) * t103;
t73 = t100 * t138 - t127;
t62 = t105 * t73;
t82 = t100 * t106 + t101 * t103;
t76 = t82 * qJD(2);
t33 = -t102 * t76 - t62;
t97 = qJD(3) + qJD(5);
t149 = t33 * t97;
t134 = qJD(5) * t102;
t75 = t82 * qJD(3);
t65 = qJD(2) * t75;
t131 = qJD(2) * qJD(3);
t125 = t106 * t131;
t126 = t103 * t131;
t88 = t100 * t126;
t66 = t101 * t125 - t88;
t4 = -qJD(5) * t62 - t102 * t65 + t105 * t66 - t76 * t134;
t164 = t4 - t149;
t114 = t102 * t73 - t105 * t76;
t110 = t114 * qJD(5) - t102 * t66 - t105 * t65;
t150 = t114 * t97;
t163 = t110 - t150;
t162 = t114 * t33;
t161 = t114 ^ 2 - t33 ^ 2;
t154 = t73 * pkin(7);
t104 = sin(qJ(2));
t133 = t104 * qJD(1);
t89 = qJD(2) * pkin(6) + t133;
t122 = qJ(4) * qJD(2) + t89;
t70 = t122 * t106;
t143 = t101 * t70;
t69 = t122 * t103;
t54 = qJD(3) * pkin(3) - t69;
t18 = t100 * t54 + t143;
t12 = t18 - t154;
t107 = cos(qJ(2));
t132 = t107 * qJD(1);
t96 = -t106 * pkin(3) - pkin(2);
t79 = t96 * qJD(2) + qJD(4) - t132;
t40 = t73 * pkin(4) + t79;
t160 = t12 * t134 - t40 * t33;
t118 = qJD(4) + t132;
t136 = qJD(3) * t103;
t35 = -t89 * t136 + (-qJ(4) * t136 + t118 * t106) * qJD(2);
t135 = qJD(3) * t106;
t36 = -t89 * t135 + (-qJ(4) * t135 - t118 * t103) * qJD(2);
t6 = -t100 * t35 + t101 * t36;
t2 = -t66 * pkin(7) + t6;
t7 = t100 * t36 + t101 * t35;
t3 = -t65 * pkin(7) + t7;
t159 = -t102 * t3 + t105 * t2 + t40 * t114;
t158 = -0.2e1 * t131;
t148 = qJ(4) + pkin(6);
t123 = qJD(3) * t148;
t71 = t106 * qJD(4) - t103 * t123;
t72 = -t103 * qJD(4) - t106 * t123;
t147 = -t100 * t71 + t101 * t72 + t82 * t132;
t81 = t100 * t103 - t142;
t112 = t81 * t107;
t146 = qJD(1) * t112 + t100 * t72 + t101 * t71;
t108 = qJD(3) ^ 2;
t109 = qJD(2) ^ 2;
t157 = (t108 + t109) * t104;
t78 = t81 * qJD(3);
t156 = qJD(5) - t97;
t155 = pkin(3) * t136 - t133;
t153 = t76 * pkin(7);
t152 = pkin(3) * t100;
t151 = pkin(3) * t103;
t50 = t100 * t70;
t20 = -t101 * t69 - t50;
t86 = t148 * t103;
t87 = t148 * t106;
t42 = -t100 * t86 + t101 * t87;
t80 = pkin(3) * t126 + qJD(2) * t133;
t145 = t103 ^ 2 - t106 ^ 2;
t144 = qJD(2) * pkin(2);
t141 = t108 * t103;
t140 = t108 * t106;
t137 = qJD(2) * t104;
t129 = pkin(3) * t138;
t17 = t101 * t54 - t50;
t19 = t100 * t69 - t143;
t41 = -t100 * t87 - t101 * t86;
t121 = t75 * pkin(4) + t155;
t119 = t107 * t158;
t117 = -t78 * pkin(7) + qJD(5) * (-t81 * pkin(7) + t42) - t147;
t116 = t75 * pkin(7) - qJD(5) * (-t82 * pkin(7) + t41) - t146;
t11 = qJD(3) * pkin(4) - t153 + t17;
t115 = -t102 * t11 - t105 * t12;
t38 = t102 * t82 + t105 * t81;
t39 = -t102 * t81 + t105 * t82;
t113 = qJD(2) * t144;
t67 = t82 * t104;
t111 = -0.2e1 * qJD(3) * t144;
t94 = t101 * pkin(3) + pkin(4);
t68 = t81 * t104;
t58 = t81 * pkin(4) + t96;
t45 = t76 * pkin(4) + t129;
t37 = t65 * pkin(4) + t80;
t22 = -qJD(2) * t112 - qJD(3) * t67;
t21 = t104 * t78 - t107 * t76;
t14 = t20 - t153;
t13 = t19 + t154;
t9 = t39 * qJD(5) - t102 * t78 + t105 * t75;
t8 = -t38 * qJD(5) - t102 * t75 - t105 * t78;
t1 = [0, 0, -t109 * t104, -t109 * t107, 0, 0, 0, 0, 0, t103 * t119 - t106 * t157, t103 * t157 + t106 * t119, t21 * qJD(3) - t107 * t65 + t73 * t137, -t22 * qJD(3) - t107 * t66 + t76 * t137, -t21 * t76 - t22 * t73 + t68 * t65 + t67 * t66, -t80 * t107 + t79 * t137 + t17 * t21 + t18 * t22 - t6 * t67 - t7 * t68, 0, 0, 0, 0, 0, (-t102 * t22 + t105 * t21 + (t102 * t67 + t105 * t68) * qJD(5)) * t97 - t33 * t137 + t107 * t110, -(t102 * t21 + t105 * t22 + (t102 * t68 - t105 * t67) * qJD(5)) * t97 - t114 * t137 - t107 * t4; 0, 0, 0, 0, 0.2e1 * t103 * t125, t145 * t158, t140, -t141, 0, -pkin(6) * t140 + t103 * t111, pkin(6) * t141 + t106 * t111, -t73 * t133 + t96 * t65 + t79 * t75 + t80 * t81 + (t73 * t151 + t147) * qJD(3), -t76 * t133 + t96 * t66 - t79 * t78 + t80 * t82 + (t76 * t151 - t146) * qJD(3), -t146 * t73 - t147 * t76 + t17 * t78 - t18 * t75 - t41 * t66 - t42 * t65 - t6 * t82 - t7 * t81, t146 * t18 + t147 * t17 + t155 * t79 + t6 * t41 + t7 * t42 + t80 * t96, -t114 * t8 + t4 * t39, t110 * t39 + t114 * t9 + t33 * t8 - t4 * t38, t8 * t97, -t9 * t97, 0, t37 * t38 + t40 * t9 - t58 * t110 + (t116 * t102 - t117 * t105) * t97 - t121 * t33, t37 * t39 + t58 * t4 + t40 * t8 + (t117 * t102 + t116 * t105) * t97 - t121 * t114; 0, 0, 0, 0, -t103 * t109 * t106, t145 * t109, 0, 0, 0, t103 * t113, t106 * t113, -t19 * qJD(3) - t73 * t129 - t79 * t76 + t6, t20 * qJD(3) - t76 * t129 + t79 * t73 - t7, (t18 + t19) * t76 + (-t17 + t20) * t73 + (-t100 * t65 - t101 * t66) * pkin(3), -t17 * t19 - t18 * t20 + (t100 * t7 + t101 * t6 - t79 * t138) * pkin(3), t162, t161, t164, t163, 0, -(-t102 * t14 + t105 * t13) * t97 + t45 * t33 + ((-t102 * t94 - t105 * t152) * t97 + t115) * qJD(5) + t159, -t105 * t3 - t102 * t2 + (t102 * t13 + t105 * t14) * t97 + t45 * t114 + (-(-t102 * t152 + t105 * t94) * t97 - t105 * t11) * qJD(5) + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t76 * qJD(3), -t88 + (-t73 + t127) * qJD(3), -t73 ^ 2 - t76 ^ 2, t17 * t76 + t18 * t73 + t80, 0, 0, 0, 0, 0, -t110 - t150, t4 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t161, t164, t163, 0, t156 * t115 + t159, (-t12 * t97 - t2) * t102 + (-t156 * t11 - t3) * t105 + t160;];
tauc_reg = t1;
