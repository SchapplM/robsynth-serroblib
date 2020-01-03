% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:55
% EndTime: 2019-12-31 17:04:57
% DurationCPUTime: 0.62s
% Computational Cost: add. (906->88), mult. (1854->139), div. (0->0), fcn. (2062->6), ass. (0->87)
t125 = cos(qJ(4));
t118 = -qJ(3) - pkin(5);
t72 = sin(qJ(2));
t63 = t118 * t72;
t73 = cos(qJ(2));
t64 = t118 * t73;
t69 = sin(pkin(7));
t70 = cos(pkin(7));
t130 = t63 * t70 + t64 * t69;
t56 = t69 * t73 + t70 * t72;
t132 = -t56 * pkin(6) + t130;
t136 = t125 * t132;
t79 = -t136 / 0.2e1;
t38 = -t63 * t69 + t70 * t64;
t55 = -t69 * t72 + t70 * t73;
t25 = -pkin(6) * t55 + t38;
t138 = t125 * t25;
t80 = t138 / 0.2e1;
t71 = sin(qJ(4));
t139 = t71 * t25;
t141 = -t136 - t139;
t137 = t71 * t132;
t140 = t138 - t137;
t119 = t71 * t55;
t54 = t125 * t56;
t129 = t54 + t119;
t33 = -t125 * t55 + t56 * t71;
t131 = -t129 ^ 2 + t33 ^ 2;
t135 = t131 * qJD(1);
t99 = t33 * qJD(4);
t11 = -qJD(2) * t33 - t99;
t134 = qJD(1) * t33;
t133 = qJD(3) * t33;
t100 = t129 * qJD(1);
t78 = t54 / 0.2e1;
t128 = pkin(2) * t69;
t126 = t72 * pkin(2);
t117 = qJD(2) * pkin(2);
t86 = -pkin(2) * t73 - pkin(1);
t7 = t86 * t126;
t114 = t7 * qJD(1);
t42 = -pkin(3) * t55 + t86;
t43 = pkin(3) * t56 + t126;
t8 = t129 * t42 + t33 * t43;
t113 = t8 * qJD(1);
t9 = t129 * t43 - t33 * t42;
t112 = t9 * qJD(1);
t111 = qJD(1) * t42;
t110 = qJD(1) * t73;
t109 = qJD(2) * t72;
t108 = qJD(2) * t73;
t12 = -t130 * t56 - t38 * t55;
t106 = t12 * qJD(1);
t14 = 0.2e1 * t78 + t119;
t104 = t14 * qJD(1);
t74 = t69 * t55 / 0.2e1 - t70 * t56 / 0.2e1;
t24 = (-t72 / 0.2e1 + t74) * pkin(2);
t103 = t24 * qJD(1);
t31 = t78 - t54 / 0.2e1;
t102 = t31 * qJD(1);
t101 = t31 * qJD(4);
t96 = t129 * qJD(4);
t37 = t55 ^ 2 + t56 ^ 2;
t95 = t37 * qJD(1);
t65 = -t72 ^ 2 + t73 ^ 2;
t94 = t65 * qJD(1);
t93 = pkin(1) * t72 * qJD(1);
t92 = pkin(1) * t110;
t91 = t33 * t100;
t90 = t129 * t134;
t89 = t33 * t111;
t88 = t129 * t111;
t87 = t72 * t110;
t1 = t80 - t138 / 0.2e1;
t67 = pkin(2) * t70 + pkin(3);
t52 = t125 * t128 + t67 * t71;
t77 = qJD(1) * t1 + qJD(2) * t52;
t2 = t79 + t136 / 0.2e1;
t51 = -t125 * t67 + t128 * t71;
t76 = -qJD(1) * t2 - qJD(2) * t51;
t75 = qJD(2) * t129 + qJD(4) * t14;
t45 = t52 * qJD(4);
t44 = t51 * qJD(4);
t23 = t126 / 0.2e1 + t74 * pkin(2);
t4 = 0.2e1 * t80 - t137;
t3 = -t139 + 0.2e1 * t79;
t5 = [0, 0, 0, t72 * t108, t65 * qJD(2), 0, 0, 0, -pkin(1) * t109, -pkin(1) * t108, t37 * qJD(3), qJD(2) * t7 + qJD(3) * t12, t11 * t129, (qJD(2) + qJD(4)) * t131, 0, 0, 0, qJD(2) * t8 + t42 * t96, qJD(2) * t9 - t42 * t99; 0, 0, 0, t87, t94, t108, -t109, 0, -pkin(5) * t108 - t93, pkin(5) * t109 - t92, (-t55 * t70 - t56 * t69) * t117, t114 + t23 * qJD(3) + (t130 * t69 + t38 * t70) * t117, -t90, t135, t11, -t75, 0, qJD(2) * t140 + t4 * qJD(4) + t113, qJD(2) * t141 + t3 * qJD(4) + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, qJD(2) * t23 + t106, 0, 0, 0, 0, 0, t101, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t135, t11, -qJD(2) * t14 - t96, 0, t4 * qJD(2) + t31 * qJD(3) + qJD(4) * t140 + t88, t3 * qJD(2) + qJD(4) * t141 - t89; 0, 0, 0, -t87, -t94, 0, 0, 0, t93, t92, 0, qJD(3) * t24 - t114, t90, -t135, 0, -t101, 0, -qJD(3) * t129 - qJD(4) * t1 - t113, qJD(4) * t2 - t112 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, -t100, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, 0, -t45 - t77, t44 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -qJD(2) * t24 - t106, 0, 0, 0, 0, 0, t75, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, 0, 0, 0, 0, 0, t100, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t135, 0, t31 * qJD(2), 0, qJD(2) * t1 - qJD(3) * t14 - t88, -qJD(2) * t2 + t133 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, t77, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
