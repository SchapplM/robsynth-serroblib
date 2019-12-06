% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:21
% EndTime: 2019-12-05 17:06:23
% DurationCPUTime: 0.60s
% Computational Cost: add. (534->134), mult. (1251->182), div. (0->0), fcn. (904->6), ass. (0->113)
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t59 = -t70 ^ 2 + t73 ^ 2;
t99 = -qJD(3) - qJD(4);
t90 = qJD(2) - t99;
t132 = t90 * t59;
t131 = pkin(4) / 0.2e1;
t75 = cos(qJ(3));
t122 = t75 * pkin(2);
t64 = pkin(3) + t122;
t74 = cos(qJ(4));
t51 = t74 * t64;
t71 = sin(qJ(4));
t72 = sin(qJ(3));
t61 = t71 * t72 * pkin(2);
t40 = t61 - t51;
t38 = -pkin(4) + t40;
t130 = -t38 / 0.2e1;
t123 = t74 * pkin(3);
t63 = -pkin(4) - t123;
t129 = -t63 / 0.2e1;
t128 = -t70 / 0.2e1;
t127 = t70 / 0.2e1;
t126 = -t73 / 0.2e1;
t125 = t73 / 0.2e1;
t124 = t71 * pkin(3);
t118 = t74 * t72;
t119 = t71 * t75;
t44 = (t118 + t119) * pkin(2);
t121 = t44 * t73;
t120 = t71 * t64;
t41 = pkin(2) * t118 + t120;
t37 = t41 * qJD(4);
t42 = t44 * qJD(3);
t117 = -t42 - t37;
t116 = pkin(2) * qJD(2);
t115 = pkin(2) * qJD(3);
t114 = pkin(3) * qJD(3);
t113 = pkin(3) * qJD(4);
t112 = pkin(4) * qJD(4);
t86 = t124 / 0.2e1 + t41 / 0.2e1;
t81 = -t44 / 0.2e1 + t86;
t11 = t81 * t73;
t111 = qJD(2) * t11;
t95 = -t122 / 0.2e1;
t87 = t95 + pkin(3) / 0.2e1;
t18 = (t64 / 0.2e1 + t87) * t71;
t110 = qJD(2) * t18;
t20 = t51 / 0.2e1 + t87 * t74;
t109 = qJD(2) * t20;
t108 = qJD(2) * t40;
t107 = qJD(2) * t41;
t106 = qJD(2) * t44;
t45 = t74 * t122 - t61;
t105 = qJD(2) * t45;
t104 = qJD(2) * t70;
t103 = qJD(2) * t73;
t102 = qJD(3) * t63;
t101 = qJD(5) * t70;
t69 = qJD(5) * t73;
t100 = -qJD(2) - qJD(3);
t98 = t71 * t114;
t97 = t71 * t113;
t96 = -t123 / 0.2e1;
t94 = t38 * t104;
t93 = t38 * t103;
t92 = t41 * t104;
t91 = t44 * t104;
t89 = pkin(2) * t100;
t88 = pkin(3) * t99;
t85 = t40 / 0.2e1 + t131 + t130;
t84 = -t45 / 0.2e1 + t129 + t130;
t83 = t71 * t88;
t82 = t96 + t131 + t129;
t1 = t84 * t70;
t80 = qJD(2) * t1 - t70 * t102;
t2 = t84 * t73;
t79 = qJD(2) * t2 - t73 * t102;
t10 = t81 * t70;
t78 = -qJD(2) * t10 - t70 * t98;
t17 = t82 * t70;
t5 = t85 * t70;
t77 = qJD(2) * t5 + qJD(3) * t17 + t70 * t112;
t19 = t82 * t73;
t6 = t85 * t73;
t76 = qJD(2) * t6 + qJD(3) * t19 + t73 * t112;
t68 = pkin(4) * t126;
t67 = pkin(4) * t128;
t62 = pkin(8) + t124;
t60 = t70 * t69;
t58 = t70 * t97;
t50 = t59 * qJD(5);
t49 = t63 * t125;
t48 = t63 * t127;
t43 = t45 * qJD(3);
t39 = pkin(8) + t41;
t36 = t40 * qJD(4);
t35 = t70 * t42;
t30 = t70 * t37;
t25 = t38 * t125;
t24 = t38 * t127;
t23 = t90 * t73 * t70;
t22 = t73 * t96 + t49 + t68;
t21 = t70 * t96 + t48 + t67;
t14 = t96 + t61 - t51 / 0.2e1 + t74 * t95;
t13 = -t124 / 0.2e1 - t120 / 0.2e1 + (-t118 - t119 / 0.2e1) * pkin(2);
t12 = -t121 / 0.2e1 - t86 * t73;
t9 = t44 * t127 + t86 * t70;
t8 = t40 * t125 + t25 + t68;
t7 = t40 * t127 + t24 + t67;
t4 = t45 * t126 + t25 + t49;
t3 = t45 * t128 + t24 + t48;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t72 * t115, -t75 * t115, 0, t117, -t43 + t36, t60, t50, 0, 0, 0, t38 * t101 + t117 * t73, t38 * t69 + t30 + t35; 0, 0, 0, 0, 0, t72 * t89, t75 * t89, 0, qJD(4) * t13 - t106 - t42, qJD(4) * t14 - t105 - t43, t60, t50, 0, 0, 0, qJD(4) * t12 + qJD(5) * t3 + t100 * t121, qJD(4) * t9 + qJD(5) * t4 + t35 + t91; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t13 - t107 - t37, qJD(3) * t14 + t108 + t36, t60, t50, 0, 0, 0, qJD(3) * t12 + qJD(5) * t7 + (-qJD(2) - qJD(4)) * t73 * t41, qJD(3) * t9 + qJD(5) * t8 + t30 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t132, t69, -t101, 0, qJD(3) * t3 + qJD(4) * t7 - t39 * t69 + t94, qJD(3) * t4 + qJD(4) * t8 + t39 * t101 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t72 * t116, t75 * t116, 0, -qJD(4) * t18 + t106, -qJD(4) * t20 + t105, t60, t50, 0, 0, 0, -qJD(4) * t11 - qJD(5) * t1 + t44 * t103, qJD(4) * t10 - qJD(5) * t2 - t91; 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t74 * t113, t60, t50, 0, 0, 0, t63 * t101 - t73 * t97, t63 * t69 + t58; 0, 0, 0, 0, 0, 0, 0, 0, t83 - t110, t74 * t88 - t109, t60, t50, 0, 0, 0, qJD(5) * t21 + t73 * t83 - t111, qJD(5) * t22 + t58 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t132, t69, -t101, 0, qJD(4) * t21 - t62 * t69 - t80, qJD(4) * t22 + t62 * t101 - t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t18 + t107, qJD(3) * t20 - t108, t60, t50, 0, 0, 0, qJD(3) * t11 - qJD(5) * t5 + t41 * t103, -qJD(3) * t10 - qJD(5) * t6 - t92; 0, 0, 0, 0, 0, 0, 0, 0, t98 + t110, t74 * t114 + t109, t60, t50, 0, 0, 0, -qJD(5) * t17 + t73 * t98 + t111, -qJD(5) * t19 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t50, 0, 0, 0, -pkin(4) * t101, -pkin(4) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t132, t69, -t101, 0, -pkin(8) * t69 - t77, pkin(8) * t101 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t132, 0, 0, 0, qJD(3) * t1 + qJD(4) * t5 - t94, qJD(3) * t2 + qJD(4) * t6 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t132, 0, 0, 0, qJD(4) * t17 + t80, qJD(4) * t19 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t132, 0, 0, 0, t77, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
