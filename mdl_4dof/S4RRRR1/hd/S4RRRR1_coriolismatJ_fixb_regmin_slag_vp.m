% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:15
% EndTime: 2019-12-31 17:22:17
% DurationCPUTime: 0.51s
% Computational Cost: add. (532->132), mult. (1249->181), div. (0->0), fcn. (902->6), ass. (0->112)
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t59 = -t70 ^ 2 + t73 ^ 2;
t99 = -qJD(2) - qJD(3);
t90 = qJD(1) - t99;
t131 = t90 * t59;
t130 = pkin(3) / 0.2e1;
t75 = cos(qJ(2));
t121 = t75 * pkin(1);
t64 = pkin(2) + t121;
t74 = cos(qJ(3));
t51 = t74 * t64;
t71 = sin(qJ(3));
t72 = sin(qJ(2));
t61 = t71 * t72 * pkin(1);
t40 = t61 - t51;
t38 = -pkin(3) + t40;
t129 = -t38 / 0.2e1;
t122 = t74 * pkin(2);
t63 = -pkin(3) - t122;
t128 = -t63 / 0.2e1;
t127 = -t70 / 0.2e1;
t126 = t70 / 0.2e1;
t125 = -t73 / 0.2e1;
t124 = t73 / 0.2e1;
t123 = t71 * pkin(2);
t117 = t74 * t72;
t118 = t71 * t75;
t44 = (t117 + t118) * pkin(1);
t120 = t44 * t73;
t119 = t71 * t64;
t41 = pkin(1) * t117 + t119;
t37 = t41 * qJD(3);
t42 = t44 * qJD(2);
t116 = -t42 - t37;
t115 = pkin(1) * qJD(1);
t114 = pkin(1) * qJD(2);
t113 = pkin(2) * qJD(2);
t112 = pkin(2) * qJD(3);
t111 = pkin(3) * qJD(3);
t110 = qJD(1) * t38;
t109 = qJD(2) * t63;
t108 = qJD(4) * t70;
t69 = qJD(4) * t73;
t86 = t123 / 0.2e1 + t41 / 0.2e1;
t81 = -t44 / 0.2e1 + t86;
t11 = t81 * t73;
t107 = t11 * qJD(1);
t95 = -t121 / 0.2e1;
t87 = t95 + pkin(2) / 0.2e1;
t18 = (t64 / 0.2e1 + t87) * t71;
t106 = t18 * qJD(1);
t20 = t51 / 0.2e1 + t87 * t74;
t105 = t20 * qJD(1);
t104 = t40 * qJD(1);
t103 = t41 * qJD(1);
t102 = t44 * qJD(1);
t45 = t74 * t121 - t61;
t101 = t45 * qJD(1);
t100 = -qJD(1) - qJD(2);
t98 = t71 * t112;
t97 = t71 * t113;
t96 = -t122 / 0.2e1;
t94 = t70 * t110;
t93 = t73 * t110;
t92 = t70 * t103;
t91 = t70 * t102;
t89 = pkin(1) * t100;
t88 = pkin(2) * t99;
t85 = t40 / 0.2e1 + t130 + t129;
t84 = -t45 / 0.2e1 + t128 + t129;
t83 = t71 * t88;
t82 = t96 + t130 + t128;
t1 = t84 * t70;
t80 = t1 * qJD(1) - t70 * t109;
t2 = t84 * t73;
t79 = t2 * qJD(1) - t73 * t109;
t10 = t81 * t70;
t78 = -t10 * qJD(1) - t70 * t97;
t17 = t82 * t70;
t5 = t85 * t70;
t77 = t5 * qJD(1) + t17 * qJD(2) + t70 * t111;
t19 = t82 * t73;
t6 = t85 * t73;
t76 = t6 * qJD(1) + t19 * qJD(2) + t73 * t111;
t68 = pkin(3) * t125;
t67 = pkin(3) * t127;
t62 = pkin(7) + t123;
t60 = t70 * t69;
t58 = t70 * t98;
t50 = t59 * qJD(4);
t49 = t63 * t124;
t48 = t63 * t126;
t43 = t45 * qJD(2);
t39 = pkin(7) + t41;
t36 = t40 * qJD(3);
t35 = t70 * t42;
t30 = t70 * t37;
t25 = t38 * t124;
t24 = t38 * t126;
t23 = t90 * t73 * t70;
t22 = t73 * t96 + t49 + t68;
t21 = t70 * t96 + t48 + t67;
t14 = t96 + t61 - t51 / 0.2e1 + t74 * t95;
t13 = -t123 / 0.2e1 - t119 / 0.2e1 + (-t117 - t118 / 0.2e1) * pkin(1);
t12 = -t120 / 0.2e1 - t86 * t73;
t9 = t44 * t126 + t86 * t70;
t8 = t40 * t124 + t25 + t68;
t7 = t40 * t126 + t24 + t67;
t4 = t45 * t125 + t25 + t49;
t3 = t45 * t127 + t24 + t48;
t15 = [0, 0, 0, 0, -t72 * t114, -t75 * t114, 0, t116, -t43 + t36, t60, t50, 0, 0, 0, t38 * t108 + t116 * t73, t38 * t69 + t30 + t35; 0, 0, 0, 0, t72 * t89, t75 * t89, 0, t13 * qJD(3) - t102 - t42, t14 * qJD(3) - t101 - t43, t60, t50, 0, 0, 0, t12 * qJD(3) + t3 * qJD(4) + t100 * t120, t9 * qJD(3) + t4 * qJD(4) + t35 + t91; 0, 0, 0, 0, 0, 0, 0, t13 * qJD(2) - t103 - t37, t14 * qJD(2) + t104 + t36, t60, t50, 0, 0, 0, t12 * qJD(2) + t7 * qJD(4) + (-qJD(1) - qJD(3)) * t73 * t41, t9 * qJD(2) + t8 * qJD(4) + t30 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t131, t69, -t108, 0, t3 * qJD(2) + t7 * qJD(3) - t39 * t69 + t94, t4 * qJD(2) + t8 * qJD(3) + t39 * t108 + t93; 0, 0, 0, 0, t72 * t115, t75 * t115, 0, -t18 * qJD(3) + t102, -t20 * qJD(3) + t101, t60, t50, 0, 0, 0, -t11 * qJD(3) - t1 * qJD(4) + t102 * t73, t10 * qJD(3) - t2 * qJD(4) - t91; 0, 0, 0, 0, 0, 0, 0, -t98, -t74 * t112, t60, t50, 0, 0, 0, t63 * t108 - t73 * t98, t63 * t69 + t58; 0, 0, 0, 0, 0, 0, 0, t83 - t106, t74 * t88 - t105, t60, t50, 0, 0, 0, t21 * qJD(4) + t73 * t83 - t107, t22 * qJD(4) + t58 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t131, t69, -t108, 0, t21 * qJD(3) - t62 * t69 - t80, t22 * qJD(3) + t62 * t108 - t79; 0, 0, 0, 0, 0, 0, 0, t18 * qJD(2) + t103, t20 * qJD(2) - t104, t60, t50, 0, 0, 0, t11 * qJD(2) - t5 * qJD(4) + t103 * t73, -t10 * qJD(2) - t6 * qJD(4) - t92; 0, 0, 0, 0, 0, 0, 0, t97 + t106, t74 * t113 + t105, t60, t50, 0, 0, 0, -t17 * qJD(4) + t73 * t97 + t107, -t19 * qJD(4) + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t50, 0, 0, 0, -pkin(3) * t108, -pkin(3) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t131, t69, -t108, 0, -pkin(7) * t69 - t77, pkin(7) * t108 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t131, 0, 0, 0, t1 * qJD(2) + t5 * qJD(3) - t94, t2 * qJD(2) + t6 * qJD(3) - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t131, 0, 0, 0, t17 * qJD(3) + t80, t19 * qJD(3) + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t131, 0, 0, 0, t77, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
