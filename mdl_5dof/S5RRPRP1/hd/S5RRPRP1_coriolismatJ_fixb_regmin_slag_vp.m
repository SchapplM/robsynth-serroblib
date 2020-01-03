% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:59:14
% EndTime: 2020-01-03 11:59:15
% DurationCPUTime: 0.61s
% Computational Cost: add. (1022->104), mult. (1939->153), div. (0->0), fcn. (1546->6), ass. (0->101)
t81 = sin(pkin(8));
t84 = sin(qJ(2));
t73 = t81 * t84 * pkin(1);
t86 = cos(qJ(2));
t130 = pkin(1) * t86;
t76 = pkin(2) + t130;
t82 = cos(pkin(8));
t98 = t76 * t82 - t73;
t53 = -pkin(3) - t98;
t75 = -pkin(2) * t82 - pkin(3);
t137 = t53 / 0.2e1 + t75 / 0.2e1;
t106 = qJD(1) + qJD(2);
t83 = sin(qJ(4));
t79 = t83 ^ 2;
t85 = cos(qJ(4));
t80 = t85 ^ 2;
t70 = t79 + t80;
t135 = t106 * t70;
t71 = t80 - t79;
t134 = t106 * t71;
t58 = t130 * t82 - t73;
t133 = -t58 / 0.2e1;
t129 = pkin(4) * t83;
t128 = pkin(4) * t85;
t122 = t82 * t84;
t95 = pkin(1) * t122 + t76 * t81;
t54 = pkin(7) + t95;
t116 = qJ(5) + t54;
t35 = t116 * t83;
t127 = t35 * t83;
t36 = t116 * t85;
t126 = t36 * t85;
t74 = pkin(2) * t81 + pkin(7);
t115 = qJ(5) + t74;
t59 = t115 * t83;
t125 = t59 * t83;
t60 = t115 * t85;
t124 = t60 * t85;
t123 = t81 * t86;
t26 = t70 * t58;
t66 = t70 * qJD(5);
t121 = qJD(2) * t26 + t66;
t120 = pkin(1) * qJD(1);
t119 = pkin(1) * qJD(2);
t118 = qJD(4) * pkin(4);
t14 = t126 + t127;
t39 = t53 - t128;
t57 = (t122 + t123) * pkin(1);
t4 = t14 * t58 + t39 * t57;
t117 = t4 * qJD(1);
t114 = qJD(1) * t14;
t113 = qJD(1) * t26;
t112 = qJD(1) * t83;
t111 = qJD(1) * t85;
t110 = qJD(2) * t83;
t109 = qJD(2) * t85;
t108 = qJD(4) * t83;
t78 = qJD(4) * t85;
t10 = -t57 * t98 + t58 * t95;
t107 = t10 * qJD(1);
t105 = pkin(4) * t78;
t104 = qJD(5) * t129;
t102 = t53 * t112;
t101 = t53 * t111;
t100 = t57 * t112;
t99 = t83 * t133;
t97 = pkin(1) * t106;
t96 = t106 * t83;
t94 = t133 - t137;
t28 = t124 + t125;
t88 = (t123 / 0.2e1 + t122 / 0.2e1) * pkin(1);
t5 = (-t60 / 0.2e1 - t36 / 0.2e1) * t85 + (-t59 / 0.2e1 - t35 / 0.2e1) * t83 + t88;
t93 = -qJD(1) * t5 + qJD(2) * t28;
t31 = t126 / 0.2e1;
t12 = t31 - t126 / 0.2e1;
t3 = t39 * t129;
t92 = qJD(1) * t3 + qJD(3) * t12;
t50 = t124 / 0.2e1;
t25 = t50 - t124 / 0.2e1;
t91 = -qJD(1) * t12 - qJD(2) * t25;
t15 = t94 * t83;
t90 = qJD(1) * t15 - t110 * t75;
t16 = t94 * t85;
t89 = qJD(1) * t16 - t109 * t75;
t64 = t75 - t128;
t1 = (t133 - t64 / 0.2e1 - t39 / 0.2e1) * t129;
t9 = t64 * t129;
t87 = -qJD(1) * t1 + qJD(2) * t9 + qJD(3) * t25;
t77 = pkin(4) * t108;
t72 = t83 * t78;
t67 = t71 * qJD(4);
t61 = pkin(4) * t96;
t56 = t85 * t96;
t43 = t57 * t110;
t23 = t25 * qJD(4);
t18 = (t133 + t137) * t85;
t17 = t137 * t83 + t99;
t11 = t12 * qJD(4);
t6 = t50 + t31 + t125 / 0.2e1 + t127 / 0.2e1 + t88;
t2 = pkin(4) * t99 + (t39 + t64) * t129 / 0.2e1;
t7 = [0, 0, 0, 0, -t84 * t119, -t86 * t119, t10 * qJD(2), t72, t67, 0, 0, 0, t108 * t53 - t109 * t57, t53 * t78 + t43, t121, qJD(2) * t4 + qJD(4) * t3 + qJD(5) * t14; 0, 0, 0, 0, -t84 * t97, -t86 * t97, t107 + (-t57 * t82 + t58 * t81) * qJD(2) * pkin(2), t72, t67, 0, 0, 0, -t106 * t57 * t85 + qJD(4) * t17, qJD(4) * t18 + t100 + t43, t113 + t121, t117 + (t28 * t58 + t57 * t64) * qJD(2) + t2 * qJD(4) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, t56, t134, t78, -t108, 0, qJD(2) * t17 - t54 * t78 + t102, qJD(2) * t18 + t108 * t54 + t101, -t105, qJD(2) * t2 - t118 * t36 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, qJD(2) * t6 + t114; 0, 0, 0, 0, t84 * t120, t86 * t120, -t107, t72, t67, 0, 0, 0, -qJD(4) * t15 + t111 * t57, -qJD(4) * t16 - t100, t66 - t113, -qJD(4) * t1 - qJD(5) * t5 - t117; 0, 0, 0, 0, 0, 0, 0, t72, t67, 0, 0, 0, t75 * t108, t75 * t78, t66, qJD(4) * t9 + qJD(5) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, t56, t134, t78, -t108, 0, -t74 * t78 - t90, t108 * t74 - t89, -t105, -t118 * t60 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t78, 0, -t77 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t56, -t134, 0, 0, 0, qJD(2) * t15 - t102, qJD(2) * t16 - t101, 0, qJD(2) * t1 - t104 - t92; 0, 0, 0, 0, 0, 0, 0, -t56, -t134, 0, 0, 0, t90, t89, 0, -t87 - t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, qJD(2) * t5 - t114 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t77 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
