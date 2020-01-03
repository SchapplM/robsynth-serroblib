% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:23
% EndTime: 2020-01-03 11:52:25
% DurationCPUTime: 0.67s
% Computational Cost: add. (1112->139), mult. (2155->178), div. (0->0), fcn. (1784->8), ass. (0->115)
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t70 = -t80 ^ 2 + t83 ^ 2;
t107 = -qJD(3) - qJD(4);
t98 = qJD(1) - t107;
t140 = t98 * t70;
t139 = pkin(4) / 0.2e1;
t129 = cos(qJ(3));
t132 = pkin(1) * sin(pkin(9));
t82 = sin(qJ(3));
t92 = cos(pkin(9)) * pkin(1) + pkin(2);
t60 = t129 * t132 + t82 * t92;
t81 = sin(qJ(4));
t125 = t81 * t60;
t59 = -t129 * t92 + t82 * t132;
t58 = pkin(3) - t59;
t84 = cos(qJ(4));
t54 = t84 * t58;
t35 = -t54 + t125;
t33 = -pkin(4) + t35;
t138 = -t33 / 0.2e1;
t130 = t84 * pkin(3);
t73 = -pkin(4) - t130;
t137 = -t73 / 0.2e1;
t136 = -t80 / 0.2e1;
t135 = t80 / 0.2e1;
t134 = -t83 / 0.2e1;
t133 = t83 / 0.2e1;
t131 = t81 * pkin(3);
t126 = t81 * t59;
t55 = t84 * t60;
t37 = t55 - t126;
t128 = t37 * t83;
t127 = t81 * t58;
t124 = t84 * t59;
t36 = t55 + t127;
t30 = t36 * qJD(4);
t31 = t37 * qJD(3);
t123 = -t31 - t30;
t122 = pkin(3) * qJD(3);
t121 = pkin(3) * qJD(4);
t120 = pkin(4) * qJD(4);
t96 = t131 / 0.2e1 + t36 / 0.2e1;
t90 = -t37 / 0.2e1 + t96;
t7 = t90 * t83;
t119 = t7 * qJD(1);
t118 = qJD(1) * t33;
t117 = qJD(3) * t73;
t116 = qJD(5) * t80;
t78 = qJD(5) * t83;
t103 = t59 / 0.2e1 + pkin(3) / 0.2e1;
t13 = (t58 / 0.2e1 + t103) * t81;
t115 = t13 * qJD(1);
t15 = t54 / 0.2e1 + t103 * t84;
t114 = t15 * qJD(1);
t113 = t35 * qJD(1);
t112 = t36 * qJD(1);
t111 = t37 * qJD(1);
t38 = -t124 - t125;
t110 = t38 * qJD(1);
t109 = t59 * qJD(1);
t108 = t60 * qJD(1);
t106 = t81 * t121;
t105 = t81 * t122;
t104 = -t130 / 0.2e1;
t102 = t80 * t118;
t101 = t83 * t118;
t100 = t80 * t112;
t99 = t80 * t111;
t97 = pkin(3) * t107;
t95 = t35 / 0.2e1 + t139 + t138;
t94 = -t38 / 0.2e1 + t137 + t138;
t93 = t81 * t97;
t91 = t104 + t139 + t137;
t1 = t94 * t80;
t89 = t1 * qJD(1) - t80 * t117;
t2 = t94 * t83;
t88 = t2 * qJD(1) - t83 * t117;
t6 = t90 * t80;
t87 = -t6 * qJD(1) - t80 * t105;
t41 = t91 * t80;
t9 = t95 * t80;
t86 = t9 * qJD(1) + t41 * qJD(3) + t80 * t120;
t10 = t95 * t83;
t42 = t91 * t83;
t85 = t10 * qJD(1) + t42 * qJD(3) + t83 * t120;
t77 = pkin(4) * t134;
t76 = pkin(4) * t136;
t72 = pkin(8) + t131;
t71 = t80 * t78;
t69 = t80 * t106;
t64 = t70 * qJD(5);
t62 = t73 * t133;
t61 = t73 * t135;
t57 = t60 * qJD(3);
t56 = t59 * qJD(3);
t47 = t98 * t83 * t80;
t44 = t83 * t104 + t62 + t77;
t43 = t80 * t104 + t61 + t76;
t34 = pkin(8) + t36;
t32 = t38 * qJD(3);
t29 = t35 * qJD(4);
t28 = t80 * t31;
t23 = t80 * t30;
t18 = t33 * t133;
t17 = t33 * t135;
t16 = t104 + t125 - t54 / 0.2e1 + t124 / 0.2e1;
t14 = -t131 / 0.2e1 - t55 - t127 / 0.2e1 + t126 / 0.2e1;
t12 = t35 * t133 + t18 + t77;
t11 = t35 * t135 + t17 + t76;
t8 = -t128 / 0.2e1 - t96 * t83;
t5 = t37 * t135 + t96 * t80;
t4 = t38 * t134 + t18 + t62;
t3 = t38 * t136 + t17 + t61;
t19 = [0, 0, 0, 0, 0, -t57, t56, 0, t123, -t32 + t29, t71, t64, 0, 0, 0, t33 * t116 + t123 * t83, t33 * t78 + t23 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t57 - t108, t56 + t109, 0, t14 * qJD(4) - t111 - t31, t16 * qJD(4) - t110 - t32, t71, t64, 0, 0, 0, t8 * qJD(4) + t3 * qJD(5) + (-qJD(1) - qJD(3)) * t128, t5 * qJD(4) + t4 * qJD(5) + t28 + t99; 0, 0, 0, 0, 0, 0, 0, 0, t14 * qJD(3) - t112 - t30, t16 * qJD(3) + t113 + t29, t71, t64, 0, 0, 0, t8 * qJD(3) + t11 * qJD(5) + (-qJD(1) - qJD(4)) * t83 * t36, t5 * qJD(3) + t12 * qJD(5) + t100 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t140, t78, -t116, 0, t3 * qJD(3) + t11 * qJD(4) - t34 * t78 + t102, t4 * qJD(3) + t12 * qJD(4) + t34 * t116 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t78; 0, 0, 0, 0, 0, t108, -t109, 0, -t13 * qJD(4) + t111, -t15 * qJD(4) + t110, t71, t64, 0, 0, 0, -t7 * qJD(4) - t1 * qJD(5) + t83 * t111, t6 * qJD(4) - t2 * qJD(5) - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t84 * t121, t71, t64, 0, 0, 0, -t83 * t106 + t73 * t116, t73 * t78 + t69; 0, 0, 0, 0, 0, 0, 0, 0, t93 - t115, t84 * t97 - t114, t71, t64, 0, 0, 0, t43 * qJD(5) + t83 * t93 - t119, t44 * qJD(5) + t69 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t140, t78, -t116, 0, t43 * qJD(4) - t72 * t78 - t89, t44 * qJD(4) + t72 * t116 - t88; 0, 0, 0, 0, 0, 0, 0, 0, t13 * qJD(3) + t112, t15 * qJD(3) - t113, t71, t64, 0, 0, 0, t7 * qJD(3) - t9 * qJD(5) + t83 * t112, -t6 * qJD(3) - t10 * qJD(5) - t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t105 + t115, t84 * t122 + t114, t71, t64, 0, 0, 0, -t41 * qJD(5) + t83 * t105 + t119, -t42 * qJD(5) + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t64, 0, 0, 0, -pkin(4) * t116, -pkin(4) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t140, t78, -t116, 0, -pkin(8) * t78 - t86, pkin(8) * t116 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t140, 0, 0, 0, t1 * qJD(3) + t9 * qJD(4) - t102, t2 * qJD(3) + t10 * qJD(4) - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t140, 0, 0, 0, t41 * qJD(4) + t89, t42 * qJD(4) + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t140, 0, 0, 0, t86, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
