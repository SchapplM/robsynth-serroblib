% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:20
% DurationCPUTime: 0.73s
% Computational Cost: add. (723->138), mult. (1953->212), div. (0->0), fcn. (1484->6), ass. (0->93)
t78 = cos(qJ(5));
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t48 = t73 * t77 + t74 * t79;
t84 = qJD(1) * t48;
t117 = t78 * t84;
t114 = qJD(1) * t73;
t101 = t79 * t114;
t113 = qJD(1) * t74;
t103 = t77 * t113;
t44 = t101 - t103;
t76 = sin(qJ(5));
t18 = t76 * t44 + t117;
t72 = qJD(4) + qJD(5);
t129 = t18 * t72;
t108 = qJD(5) * t76;
t40 = t48 * qJD(4);
t35 = qJD(1) * t40;
t55 = qJD(4) * t103;
t36 = qJD(4) * t101 - t55;
t85 = qJD(5) * t117 + t44 * t108 + t78 * t35 + t76 * t36;
t135 = -t85 + t129;
t93 = t78 * t44 - t76 * t84;
t134 = t93 * t18;
t111 = qJD(2) * t79;
t112 = qJD(2) * t77;
t133 = t73 * t111 - t74 * t112;
t132 = -t18 ^ 2 + t93 ^ 2;
t56 = qJ(2) * t114 + qJD(3);
t47 = -pkin(6) * t114 + t56;
t116 = -pkin(6) + qJ(2);
t53 = t116 * t74;
t50 = qJD(1) * t53;
t92 = -t77 * t47 - t79 * t50;
t13 = -pkin(7) * t84 - t92;
t38 = -qJD(1) * pkin(1) - pkin(2) * t113 - qJ(3) * t114 + qJD(2);
t27 = pkin(3) * t113 - t38;
t16 = pkin(4) * t84 + t27;
t88 = t133 * qJD(1);
t5 = t35 * pkin(7) + qJD(4) * t92 + t88;
t131 = t16 * t18 + t13 * t108 + (-t13 * t72 - t5) * t76;
t128 = t72 * t93;
t127 = 0.2e1 * t84;
t126 = t79 * t47 - t77 * t50;
t125 = qJD(5) - t72;
t2 = qJD(5) * t93 - t76 * t35 + t78 * t36;
t124 = -t2 + t128;
t4 = -t36 * pkin(7) + qJD(2) * t84 + t126 * qJD(4);
t123 = -t16 * t93 - t76 * t4 + t78 * t5;
t122 = t72 ^ 2;
t121 = pkin(4) * t44;
t118 = t78 * t13;
t70 = t73 ^ 2;
t71 = t74 ^ 2;
t115 = t70 + t71;
t110 = qJD(4) * t77;
t109 = qJD(4) * t79;
t107 = t73 * qJD(3);
t106 = qJD(1) * qJD(2);
t105 = qJD(1) * qJD(3);
t104 = -t74 * pkin(2) - t73 * qJ(3) - pkin(1);
t12 = -t44 * pkin(7) + t126;
t11 = qJD(4) * pkin(4) + t12;
t99 = -pkin(4) * t72 - t11;
t81 = qJD(1) ^ 2;
t51 = t115 * t81;
t98 = t73 * t105;
t96 = qJ(2) * t106;
t45 = t74 * pkin(3) - t104;
t95 = t71 * t96;
t49 = t73 * t79 - t74 * t77;
t22 = t78 * t48 + t76 * t49;
t23 = -t76 * t48 + t78 * t49;
t52 = t116 * t73;
t91 = -t77 * t52 - t79 * t53;
t86 = t52 * t109 - t53 * t110 + t74 * t111 + t73 * t112;
t82 = qJD(4) * t91 + t133;
t80 = qJD(4) ^ 2;
t57 = t70 * t96;
t41 = t73 * t109 - t74 * t110;
t39 = 0.2e1 * t115 * t106;
t26 = t41 * pkin(4) + t107;
t25 = t36 * pkin(4) + t98;
t24 = t48 * pkin(4) + t45;
t15 = -t48 * pkin(7) - t91;
t14 = -t49 * pkin(7) + t79 * t52 - t77 * t53;
t9 = t40 * pkin(7) + t82;
t8 = -t41 * pkin(7) + t86;
t7 = qJD(5) * t23 - t76 * t40 + t78 * t41;
t6 = -qJD(5) * t22 - t78 * t40 - t76 * t41;
t1 = [0, 0, 0, 0, 0, t39, 0.2e1 * t57 + 0.2e1 * t95, 0.2e1 * t74 * t98, t39, 0.2e1 * t70 * t105, 0.2e1 * t95 + t57 + (t56 * qJD(2) + (-qJD(1) * t104 - t38) * qJD(3)) * t73, -t35 * t49 - t44 * t40, t35 * t48 - t49 * t36 + t40 * t84 - t44 * t41, -t40 * qJD(4), -t41 * qJD(4), 0, t82 * qJD(4) + t127 * t107 + t27 * t41 + t45 * t36, -t45 * t35 - t27 * t40 - t86 * qJD(4) + (qJD(1) * t49 + t44) * t107, -t23 * t85 + t6 * t93, -t6 * t18 - t23 * t2 + t22 * t85 - t7 * t93, t6 * t72, -t7 * t72, 0, t26 * t18 + t24 * t2 + t25 * t22 + t16 * t7 + (-t76 * t8 + t78 * t9 + (-t14 * t76 - t15 * t78) * qJD(5)) * t72, t26 * t93 - t24 * t85 + t25 * t23 + t16 * t6 - (t76 * t9 + t78 * t8 + (t14 * t78 - t15 * t76) * qJD(5)) * t72; 0, 0, 0, 0, 0, -t51, -qJ(2) * t51, 0, -t51, 0, -t71 * t81 * qJ(2) + (-qJD(3) - t56) * t114, 0, 0, 0, 0, 0, t55 + (-t44 - t101) * qJD(4), t127 * qJD(4), 0, 0, 0, 0, 0, -t2 - t128, t85 + t129; 0, 0, 0, 0, 0, 0, 0, -t73 * t81 * t74, 0, -t70 * t81, (qJD(2) + t38) * t114, 0, 0, 0, 0, 0, -t114 * t84 - t80 * t77, -t44 * t114 - t80 * t79, 0, 0, 0, 0, 0, -t18 * t114 + (-t76 * t79 - t77 * t78) * t122, -t93 * t114 - (-t76 * t77 + t78 * t79) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t84, t44 ^ 2 - t84 ^ 2, 0, t55 + (t44 - t101) * qJD(4), 0, -t27 * t44 + t88, -(qJD(2) - t27) * t84, t134, t132, t135, t124, 0, -t18 * t121 - (-t76 * t12 - t118) * t72 + (t76 * t99 - t118) * qJD(5) + t123, -t93 * t121 + (qJD(5) * t99 + t12 * t72 - t4) * t78 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t132, t135, t124, 0, t125 * (-t76 * t11 - t118) + t123, (-t125 * t11 - t4) * t78 + t131;];
tauc_reg = t1;
