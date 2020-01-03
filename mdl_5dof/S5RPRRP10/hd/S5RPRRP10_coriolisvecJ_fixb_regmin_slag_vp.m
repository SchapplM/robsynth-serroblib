% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:08
% DurationCPUTime: 1.22s
% Computational Cost: add. (1880->213), mult. (4969->291), div. (0->0), fcn. (3624->6), ass. (0->110)
t130 = pkin(6) + qJ(2);
t91 = sin(pkin(8));
t77 = t130 * t91;
t74 = qJD(1) * t77;
t92 = cos(pkin(8));
t78 = t130 * t92;
t75 = qJD(1) * t78;
t94 = sin(qJ(3));
t96 = cos(qJ(3));
t46 = -t94 * t74 + t96 * t75;
t149 = t46 * qJD(3);
t73 = t96 * t91 + t94 * t92;
t136 = t94 * t91;
t115 = qJD(1) * t136;
t131 = t96 * t92;
t84 = qJD(1) * t131;
t67 = t84 - t115;
t61 = qJD(4) - t67;
t93 = sin(qJ(4));
t112 = t61 * t93;
t68 = t73 * qJD(1);
t95 = cos(qJ(4));
t54 = t93 * qJD(3) + t95 * t68;
t148 = t54 * t112;
t147 = -t96 * t74 - t94 * t75;
t49 = t96 * t77 + t94 * t78;
t83 = qJD(3) * t84;
t105 = qJD(3) * t115 - t83;
t119 = t95 * qJD(3);
t121 = qJD(4) * t93;
t21 = -qJD(4) * t119 + t95 * t105 + t68 * t121;
t70 = t73 * qJD(3);
t58 = qJD(1) * t70;
t86 = -t92 * pkin(2) - pkin(1);
t76 = t86 * qJD(1) + qJD(2);
t25 = -t67 * pkin(3) - t68 * pkin(7) + t76;
t41 = qJD(3) * pkin(7) + t46;
t13 = t93 * t25 + t95 * t41;
t72 = -t131 + t136;
t100 = t72 * qJD(2);
t17 = -qJD(1) * t100 + qJD(3) * t147;
t30 = t58 * pkin(3) + t105 * pkin(7);
t29 = t95 * t30;
t98 = -t13 * qJD(4) - t93 * t17 + t29;
t1 = t58 * pkin(4) + t21 * qJ(5) - t54 * qJD(5) + t98;
t52 = t93 * t68 - t119;
t7 = -t52 * qJ(5) + t13;
t146 = t61 * t7 + t1;
t145 = t54 ^ 2;
t12 = t95 * t25 - t93 * t41;
t6 = -t54 * qJ(5) + t12;
t5 = t61 * pkin(4) + t6;
t144 = t5 - t6;
t129 = -qJ(5) - pkin(7);
t113 = qJD(4) * t129;
t122 = qJ(5) * t95;
t42 = t68 * pkin(3) - t67 * pkin(7);
t34 = t95 * t42;
t143 = -t68 * pkin(4) + t95 * t113 + t67 * t122 - t34 + (-qJD(5) + t147) * t93;
t142 = t21 * t93;
t141 = t52 * t67;
t140 = t54 * t68;
t139 = t68 * t52;
t138 = t73 * t95;
t137 = t93 * t58;
t50 = -t94 * t77 + t96 * t78;
t47 = t95 * t50;
t57 = t95 * t58;
t120 = qJD(4) * t95;
t99 = t93 * t105;
t22 = t54 * qJD(4) - t99;
t128 = -t52 * t120 - t93 * t22;
t127 = t147 * t95 + t93 * t42;
t44 = t72 * pkin(3) - t73 * pkin(7) + t86;
t126 = t93 * t44 + t47;
t123 = qJ(5) * t93;
t125 = t95 * qJD(5) + t93 * t113 + t67 * t123 - t127;
t124 = t91 ^ 2 + t92 ^ 2;
t118 = qJD(1) * qJD(2);
t26 = -t49 * qJD(3) - t100;
t69 = t72 * qJD(3);
t43 = t70 * pkin(3) + t69 * pkin(7);
t117 = t44 * t120 + t95 * t26 + t93 * t43;
t116 = t73 * t120;
t114 = t124 * qJD(1) ^ 2;
t111 = t61 * t95;
t18 = t73 * t118 + t149;
t102 = t25 * t120 - t41 * t121 + t95 * t17 + t93 * t30;
t2 = -t22 * qJ(5) - t52 * qJD(5) + t102;
t110 = -t61 * t5 + t2;
t108 = qJ(5) * t69 - qJD(5) * t73;
t107 = 0.2e1 * t124 * t118;
t106 = t67 * t112 - t61 * t121 + t57;
t40 = -qJD(3) * pkin(3) - t147;
t9 = t22 * pkin(4) + t18;
t104 = -t69 * t93 + t116;
t103 = -t73 * t121 - t95 * t69;
t101 = -pkin(7) * t58 + t61 * t40;
t27 = t73 * qJD(2) + t50 * qJD(3);
t80 = t129 * t95;
t79 = t129 * t93;
t51 = t52 ^ 2;
t37 = t95 * t44;
t35 = t95 * t43;
t15 = t52 * pkin(4) + qJD(5) + t40;
t14 = -t73 * t123 + t126;
t10 = t72 * pkin(4) - t73 * t122 - t93 * t50 + t37;
t4 = -qJ(5) * t116 + (-qJD(4) * t50 + t108) * t93 + t117;
t3 = t70 * pkin(4) - t93 * t26 + t35 + t108 * t95 + (-t47 + (qJ(5) * t73 - t44) * t93) * qJD(4);
t8 = [0, 0, 0, 0, 0, t107, qJ(2) * t107, -t105 * t73 - t68 * t69, t105 * t72 - t73 * t58 - t69 * t67 - t68 * t70, -t69 * qJD(3), -t70 * qJD(3), 0, -t27 * qJD(3) + t86 * t58 + t76 * t70, -t26 * qJD(3) - t86 * t105 - t76 * t69, t103 * t54 - t21 * t138, -(-t52 * t95 - t54 * t93) * t69 + (t142 - t22 * t95 + (t52 * t93 - t54 * t95) * qJD(4)) * t73, t103 * t61 - t21 * t72 + t54 * t70 + t57 * t73, -t104 * t61 - t73 * t137 - t22 * t72 - t52 * t70, t58 * t72 + t61 * t70, (-t120 * t50 + t35) * t61 + t37 * t58 + (-t120 * t41 + t29) * t72 + t12 * t70 + t27 * t52 + t49 * t22 + t40 * t116 + ((-qJD(4) * t44 - t26) * t61 - t50 * t58 + (-qJD(4) * t25 - t17) * t72 + t18 * t73 - t40 * t69) * t93, -(-t121 * t50 + t117) * t61 - t126 * t58 - t102 * t72 - t13 * t70 + t27 * t54 - t49 * t21 + t18 * t138 + t103 * t40, t10 * t21 - t14 * t22 - t3 * t54 - t4 * t52 - (-t5 * t95 - t7 * t93) * t69 + (-t1 * t95 - t2 * t93 + (t5 * t93 - t7 * t95) * qJD(4)) * t73, t2 * t14 + t7 * t4 + t1 * t10 + t5 * t3 + t9 * (t93 * t73 * pkin(4) + t49) + t15 * (pkin(4) * t104 + t27); 0, 0, 0, 0, 0, -t114, -qJ(2) * t114, 0, 0, 0, 0, 0, 0.2e1 * t68 * qJD(3), t83 + (t67 - t115) * qJD(3), 0, 0, 0, 0, 0, t106 - t139, -t61 ^ 2 * t95 - t137 - t140, (t21 + t141) * t95 + t148 + t128, t110 * t93 + t146 * t95 - t15 * t68; 0, 0, 0, 0, 0, 0, 0, -t68 * t67, -t67 ^ 2 + t68 ^ 2, t83 + (-t67 - t115) * qJD(3), 0, 0, -t76 * t68 + t149 - t18, t72 * t118 - t76 * t67, t54 * t111 - t142, (-t21 + t141) * t95 - t148 + t128, t111 * t61 + t137 - t140, t106 + t139, -t61 * t68, -pkin(3) * t22 - t12 * t68 - t18 * t95 - t46 * t52 + (-pkin(7) * t120 - t34) * t61 + (t147 * t61 + t101) * t93, pkin(3) * t21 + t13 * t68 + t18 * t93 - t46 * t54 + (pkin(7) * t121 + t127) * t61 + t101 * t95, t110 * t95 - t125 * t52 - t143 * t54 - t146 * t93 + t79 * t21 + t80 * t22, -t2 * t80 + t1 * t79 + t9 * (-t95 * pkin(4) - pkin(3)) + t125 * t7 + t143 * t5 + (pkin(4) * t112 - t46) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t51 + t145, t52 * t61 - t21, t99 + (-qJD(4) + t61) * t54, t58, t13 * t61 - t40 * t54 + t98, t12 * t61 + t40 * t52 - t102, pkin(4) * t21 - t144 * t52, t144 * t7 + (-t15 * t54 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 - t145, t5 * t54 + t7 * t52 + t9;];
tauc_reg = t8;
