% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:17
% DurationCPUTime: 0.79s
% Computational Cost: add. (1142->168), mult. (2846->217), div. (0->0), fcn. (1830->6), ass. (0->109)
t114 = cos(pkin(8));
t81 = sin(qJ(3));
t101 = t114 * t81;
t78 = sin(pkin(8));
t82 = cos(qJ(3));
t61 = t78 * t82 + t101;
t55 = t61 * qJD(3);
t42 = qJD(1) * t55;
t112 = qJD(1) * t81;
t100 = t114 * t82;
t96 = qJD(1) * t100;
t53 = t78 * t112 - t96;
t111 = qJD(3) * t81;
t58 = qJD(3) * t100 - t78 * t111;
t117 = -t61 * t42 - t58 * t53;
t106 = qJD(1) * qJD(3);
t102 = t81 * t106;
t65 = t78 * t102;
t43 = qJD(3) * t96 - t65;
t122 = t78 * t81;
t60 = -t100 + t122;
t89 = qJD(1) * t61;
t118 = t43 * t60 + t55 * t89;
t135 = t117 + t118;
t134 = t118 - t117;
t133 = 0.2e1 * t89;
t130 = t89 ^ 2;
t49 = t53 ^ 2;
t132 = -t49 - t130;
t131 = -t49 + t130;
t110 = t81 * qJD(2);
t71 = sin(pkin(7)) * pkin(1) + pkin(6);
t63 = t71 * qJD(1);
t98 = qJ(4) * qJD(1) + t63;
t38 = t98 * t82 + t110;
t129 = pkin(3) * t81;
t108 = t82 * qJD(4);
t75 = t82 * qJD(2);
t74 = qJD(3) * t75;
t40 = -t63 * t111 + t74;
t27 = (-qJ(4) * t111 + t108) * qJD(1) + t40;
t109 = t81 * qJD(4);
t85 = -qJD(1) * t109 - t38 * qJD(3);
t2 = -t114 * t85 + t78 * t27;
t115 = qJ(4) + t71;
t59 = t115 * t82;
t23 = t115 * t101 + t78 * t59;
t128 = t2 * t23;
t127 = t2 * t60;
t73 = -cos(pkin(7)) * pkin(1) - pkin(2);
t62 = -t82 * pkin(3) + t73;
t113 = qJD(1) * t62;
t52 = qJD(4) + t113;
t17 = t53 * pkin(4) - qJ(5) * t89 + t52;
t126 = t17 * t89;
t125 = t89 * t53;
t123 = t78 * t38;
t121 = t81 * t63;
t83 = qJD(3) ^ 2;
t120 = t83 * t81;
t119 = t83 * t82;
t3 = t114 * t27 + t78 * t85;
t33 = t114 * t38;
t37 = -t98 * t81 + t75;
t35 = qJD(3) * pkin(3) + t37;
t10 = t78 * t35 + t33;
t116 = t81 ^ 2 - t82 ^ 2;
t64 = qJD(1) * t73;
t47 = qJD(3) * t55;
t14 = t114 * t37 - t123;
t107 = qJD(5) - t14;
t84 = qJD(1) ^ 2;
t105 = t81 * t84 * t82;
t104 = pkin(3) * t111;
t103 = pkin(3) * t112;
t99 = qJD(3) * t115;
t97 = t82 * t102;
t95 = t42 * t60 + t53 * t55;
t45 = t82 * t63 + t110;
t93 = 0.2e1 * qJD(3) * t64;
t68 = pkin(3) * t102;
t92 = t42 * pkin(4) - t43 * qJ(5) + t68;
t13 = t78 * t37 + t33;
t91 = t13 * qJD(3) - t2;
t9 = t114 * t35 - t123;
t90 = -t82 * t99 - t109;
t39 = -t81 * t99 + t108;
t15 = -t114 * t90 + t78 * t39;
t16 = t114 * t39 + t78 * t90;
t24 = t114 * t59 - t115 * t122;
t88 = t15 * t89 - t16 * t53 + t2 * t61 + t23 * t43 - t24 * t42;
t87 = t133 * qJD(3);
t41 = t45 * qJD(3);
t44 = t75 - t121;
t86 = t40 * t82 + t41 * t81 + (-t44 * t82 - t45 * t81) * qJD(3);
t72 = -t114 * pkin(3) - pkin(4);
t69 = t78 * pkin(3) + qJ(5);
t48 = t58 * qJD(3);
t26 = -t65 + (t96 + t53) * qJD(3);
t25 = -t65 + (t96 - t53) * qJD(3);
t19 = t60 * pkin(4) - t61 * qJ(5) + t62;
t18 = pkin(4) * t89 + t53 * qJ(5) + t103;
t11 = t55 * pkin(4) - t58 * qJ(5) - t61 * qJD(5) + t104;
t8 = qJD(3) * qJ(5) + t10;
t7 = -qJD(3) * pkin(4) + qJD(5) - t9;
t6 = t43 * t61 + t58 * t89;
t5 = -qJD(5) * t89 + t92;
t1 = qJD(3) * qJD(5) + t3;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97, -0.2e1 * t116 * t106, t119, -0.2e1 * t97, -t120, 0, -t71 * t119 + t81 * t93, t71 * t120 + t82 * t93, t86, t86 * t71, t6, -t134, t48, t95, -t47, 0, t62 * t42 + t52 * t55 + (-t15 + (qJD(1) * t60 + t53) * t129) * qJD(3), t62 * t43 + t52 * t58 + (t133 * t129 - t16) * qJD(3), -t10 * t55 - t3 * t60 - t9 * t58 + t88, t10 * t16 - t9 * t15 + t128 + t3 * t24 + (t52 + t113) * t104, t6, t48, t134, 0, t47, t95, -t15 * qJD(3) + t11 * t53 + t17 * t55 + t19 * t42 + t5 * t60, -t1 * t60 - t8 * t55 + t7 * t58 + t88, t16 * qJD(3) - t11 * t89 - t17 * t58 - t19 * t43 - t5 * t61, t1 * t24 + t17 * t11 + t7 * t15 + t8 * t16 + t5 * t19 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, t40 * t81 - t41 * t82 + (-t44 * t81 + t45 * t82) * qJD(3), 0, 0, 0, 0, 0, 0, -t47, -t48, t135, t10 * t58 + t3 * t61 - t9 * t55 + t127, 0, 0, 0, 0, 0, 0, -t47, t135, t48, t1 * t61 + t7 * t55 + t8 * t58 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t116 * t84, 0, t105, 0, 0, -t64 * t112, -t64 * t82 * qJD(1) - t74 + (t44 + t121) * qJD(3), 0, 0, t125, t131, t26, -t125, 0, 0, -t53 * t103 - t52 * t89 + t91, t14 * qJD(3) - t103 * t89 + t52 * t53 - t3, (t10 - t13) * t89 + (t14 - t9) * t53 + (-t114 * t43 - t42 * t78) * pkin(3), -t10 * t14 + t9 * t13 + (-t52 * t112 - t114 * t2 + t3 * t78) * pkin(3), t125, t26, -t131, 0, 0, -t125, -t18 * t53 - t126 + t91, -t69 * t42 + t72 * t43 + (-t13 + t8) * t89 + (t7 - t107) * t53, -t17 * t53 + t18 * t89 + (0.2e1 * qJD(5) - t14) * qJD(3) + t3, t1 * t69 + t107 * t8 - t7 * t13 - t17 * t18 + t2 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t25, t132, t10 * t53 + t89 * t9 + t68, 0, 0, 0, 0, 0, 0, t87, t132, -t25, t8 * t53 + (-qJD(5) - t7) * t89 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t26, -t130 - t83, -t8 * qJD(3) + t126 + t2;];
tauc_reg = t4;
