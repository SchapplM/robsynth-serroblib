% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:39
% EndTime: 2019-12-31 18:24:43
% DurationCPUTime: 0.90s
% Computational Cost: add. (456->119), mult. (970->182), div. (0->0), fcn. (793->6), ass. (0->116)
t59 = sin(qJ(5));
t130 = 0.2e1 * t59;
t62 = cos(qJ(3));
t63 = -pkin(3) - pkin(7);
t129 = t62 * t63;
t54 = t59 ^ 2;
t61 = cos(qJ(5));
t56 = t61 ^ 2;
t37 = t54 - t56;
t98 = t62 * qJD(1);
t93 = t61 * t98;
t66 = t37 * qJD(3) + t93 * t130;
t60 = sin(qJ(3));
t115 = t60 * qJ(4);
t97 = t62 * qJD(4);
t128 = (t115 - t129) * qJD(3) - t97;
t127 = t62 * pkin(3);
t48 = sin(pkin(8)) * pkin(1) + pkin(6);
t126 = pkin(4) + t48;
t22 = t126 * t62;
t125 = t22 * t59;
t21 = t126 * t60;
t124 = t59 * t21;
t114 = t62 * qJ(4);
t34 = -t60 * pkin(3) + t114;
t24 = t60 * pkin(7) - t34;
t123 = t59 * t24;
t122 = t61 * t21;
t121 = t61 * t24;
t120 = t61 * t62;
t55 = t60 ^ 2;
t57 = t62 ^ 2;
t38 = t57 - t55;
t49 = -cos(pkin(8)) * pkin(1) - pkin(2);
t73 = t49 - t115;
t14 = t73 + t129;
t7 = t59 * t14 - t122;
t1 = (-t7 - t122) * t62 - t123 * t60;
t119 = t1 * qJD(1);
t8 = t61 * t14 + t124;
t2 = t121 * t60 + (t8 - t124) * t62;
t118 = t2 * qJD(1);
t3 = t22 * t120 - t7 * t60;
t117 = t3 * qJD(1);
t4 = -t62 * t125 - t8 * t60;
t116 = t4 * qJD(1);
t113 = qJD(1) * t60;
t112 = qJD(1) * t61;
t111 = qJD(5) * t60;
t110 = qJD(5) * t61;
t109 = qJD(5) * t63;
t20 = t73 - t127;
t11 = t20 * t62 - t34 * t60;
t108 = t11 * qJD(1);
t12 = -t20 * t60 - t34 * t62;
t107 = t12 * qJD(1);
t28 = t38 * t59;
t106 = t28 * qJD(1);
t30 = t38 * t61;
t105 = t30 * qJD(1);
t104 = t38 * qJD(1);
t103 = t55 * qJD(1);
t102 = t55 * qJD(4);
t101 = t59 * qJD(3);
t51 = t60 * qJD(3);
t100 = t60 * qJD(4);
t99 = t61 * qJD(3);
t52 = t62 * qJD(3);
t96 = t62 * qJD(5);
t95 = qJ(4) * qJD(5);
t94 = qJD(3) * qJ(4);
t92 = t59 * t111;
t47 = t59 * t96;
t91 = t60 * t110;
t90 = t61 * t96;
t89 = t20 * t34 * qJD(1);
t88 = t20 * t113;
t87 = t49 * t113;
t86 = t49 * t98;
t85 = t59 * t52;
t43 = t60 * t52;
t84 = t48 * t51;
t42 = t60 * t98;
t83 = t59 * t110;
t82 = t59 * t99;
t81 = t120 * t130;
t80 = t59 * t42;
t78 = t62 * t82;
t77 = qJD(5) + t113;
t75 = -t103 - t111;
t74 = t34 * qJD(3) + t100;
t72 = -t63 * t60 / 0.2e1 - t114 / 0.2e1;
t65 = t24 / 0.2e1 + t72;
t9 = t65 * t59;
t71 = -t9 * qJD(1) - t61 * t94;
t23 = (t56 / 0.2e1 - t54 / 0.2e1) * t62;
t70 = t23 * qJD(1) + t82;
t10 = t65 * t61;
t69 = -t10 * qJD(1) + t59 * t94;
t68 = t59 * t57 * t112 - t23 * qJD(3);
t29 = t37 * t57;
t67 = -t29 * qJD(1) + 0.2e1 * t78;
t64 = (-t115 - t127) * qJD(3) + t97;
t50 = t52 / 0.2e1;
t46 = t61 * t52;
t45 = t60 * t99;
t44 = t60 * t112;
t41 = t59 * t113;
t31 = t48 * t52;
t27 = -t44 - t110;
t26 = -qJD(5) * t59 - t41;
t25 = t42 + t96 / 0.2e1;
t19 = t23 * qJD(5);
t6 = -t125 - t121 / 0.2e1 + t72 * t61;
t5 = t22 * t61 - t123 / 0.2e1 + t72 * t59;
t13 = [0, 0, 0, 0, t43, t38 * qJD(3), 0, 0, 0, t49 * t51, t49 * t52, 0, t12 * qJD(3) - t60 * t97, -t11 * qJD(3) + t102, -t74 * t20, -t54 * t43 + t57 * t83, -t29 * qJD(5) - 0.2e1 * t60 * t78, -t28 * qJD(3) - t60 * t90, -t30 * qJD(3) + t60 * t47, t43, t1 * qJD(3) + t4 * qJD(5) + t59 * t102, -t2 * qJD(3) - t3 * qJD(5) + t61 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t42, t104, t52, -t51, 0, -t31 + t87, t84 + t86, t64, t31 + t107, -t84 - t108, t64 * t48 - t89, -t19 + (-t54 * t98 + t82) * t60, qJD(5) * t81 - t66 * t60, t46 - t106, -t85 - t105, t25, t5 * qJD(5) - t21 * t101 - t128 * t61 + t119, t6 * qJD(5) + t128 * t59 - t21 * t99 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t42, t103, t31 - t88, 0, 0, 0, 0, 0, t59 * t103 + t46, t61 * t103 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, -t77 * t120, t47 + t80, t50, t5 * qJD(3) - t8 * qJD(5) + t116, t6 * qJD(3) + t7 * qJD(5) - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t52, 0, t51, t52, t74, 0, 0, 0, 0, 0, t85 + t91, t46 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 + t45, -t59 * t51 + t90; 0, 0, 0, 0, -t42, -t104, 0, 0, 0, -t87, -t86, 0, -t107, t108, t89, t54 * t42 - t19, t77 * t81, -t92 + t106, -t91 + t105, -t25, t9 * qJD(5) - t119, t10 * qJD(5) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t83, t37 * qJD(5), 0, 0, 0, qJD(4) * t59 + t61 * t95, qJD(4) * t61 - t59 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t94, 0, 0, 0, 0, 0, t101, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t66, t26, t27, -t98 / 0.2e1, -t59 * t109 - t71, -t61 * t109 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t103, t88, 0, 0, 0, 0, 0, t75 * t59, t75 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t94, 0, 0, 0, 0, 0, -t101, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, (t93 + t101) * t60, t45 - t80, t50, -t9 * qJD(3) + t59 * t100 - t116, -t10 * qJD(3) + t61 * t100 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t66, t41, t44, t98 / 0.2e1, t71, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
