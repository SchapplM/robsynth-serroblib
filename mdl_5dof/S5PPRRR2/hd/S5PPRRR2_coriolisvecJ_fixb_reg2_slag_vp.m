% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:54
% DurationCPUTime: 1.02s
% Computational Cost: add. (1362->165), mult. (3537->238), div. (0->0), fcn. (2715->8), ass. (0->115)
t73 = sin(pkin(9));
t74 = cos(pkin(9));
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t137 = -t77 * t73 + t80 * t74;
t45 = t137 * qJD(1);
t75 = sin(qJ(5));
t76 = sin(qJ(4));
t78 = cos(qJ(5));
t79 = cos(qJ(4));
t55 = t75 * t79 + t76 * t78;
t70 = qJD(4) + qJD(5);
t138 = t70 * t55;
t28 = t138 * qJD(3);
t53 = t73 * t80 + t74 * t77;
t124 = t78 * t79;
t129 = t75 * t76;
t54 = -t124 + t129;
t24 = t54 * t53;
t47 = t137 * qJD(3);
t43 = qJD(1) * t47;
t69 = t79 * qJD(2);
t117 = qJD(4) * t69 + t43 * t79;
t46 = t53 * qJD(1);
t42 = qJD(3) * pkin(6) + t46;
t95 = pkin(7) * qJD(3) + t42;
t88 = t95 * t76;
t10 = -qJD(4) * t88 + t117;
t30 = t69 - t88;
t26 = qJD(4) * pkin(4) + t30;
t136 = (qJD(5) * t26 + t10) * t78;
t107 = t76 * qJD(2);
t31 = t79 * t95 + t107;
t135 = -pkin(7) - pkin(6);
t109 = qJD(4) * t79;
t91 = t70 * t129;
t32 = -qJD(5) * t124 - t109 * t78 + t91;
t134 = t32 * t70;
t68 = -pkin(4) * t79 - pkin(3);
t39 = qJD(3) * t68 - t45;
t51 = t55 * qJD(3);
t133 = t39 * t51;
t48 = t53 * qJD(3);
t44 = qJD(1) * t48;
t132 = t44 * t137;
t111 = qJD(3) * t79;
t101 = t78 * t111;
t112 = qJD(3) * t76;
t102 = t75 * t112;
t49 = -t101 + t102;
t131 = t51 * t49;
t130 = t75 * t31;
t128 = t76 * t42;
t127 = t76 * t43;
t125 = t78 * t31;
t81 = qJD(4) ^ 2;
t122 = t81 * t76;
t121 = t81 * t79;
t58 = t135 * t76;
t59 = t135 * t79;
t37 = t58 * t78 + t59 * t75;
t98 = qJD(4) * t135;
t56 = t76 * t98;
t57 = t79 * t98;
t120 = qJD(5) * t37 + t45 * t54 + t78 * t56 + t75 * t57;
t38 = t58 * t75 - t59 * t78;
t119 = -qJD(5) * t38 + t45 * t55 - t75 * t56 + t78 * t57;
t118 = -t28 * t55 + t32 * t49;
t106 = qJD(3) * qJD(4);
t97 = t79 * t106;
t116 = -qJD(5) * t101 - t78 * t97;
t71 = t76 ^ 2;
t72 = t79 ^ 2;
t115 = t71 - t72;
t114 = t71 + t72;
t113 = qJD(3) * pkin(3);
t110 = qJD(4) * t76;
t108 = t47 * qJD(3);
t82 = qJD(3) ^ 2;
t105 = t76 * t82 * t79;
t104 = pkin(4) * t112;
t103 = pkin(4) * t110;
t100 = -pkin(4) * t70 - t26;
t11 = -qJD(4) * t31 - t127;
t99 = -t75 * t10 + t11 * t78;
t96 = -qJD(5) * t130 + t75 * t11;
t93 = t76 * t97;
t92 = -t46 + t103;
t6 = t26 * t75 + t125;
t27 = qJD(3) * t91 + t116;
t90 = t138 * t51 - t27 * t54;
t34 = t69 - t128;
t35 = t42 * t79 + t107;
t89 = t34 * t76 - t35 * t79;
t87 = t39 * t49 - t96;
t86 = pkin(6) * t81 - t46 * qJD(3) + t44;
t41 = -t45 - t113;
t85 = qJD(4) * (t41 + t45 - t113);
t2 = -qJD(5) * t6 + t99;
t16 = -t110 * t42 + t117;
t17 = -t35 * qJD(4) - t127;
t84 = t16 * t79 - t17 * t76 + (-t34 * t79 - t35 * t76) * qJD(4);
t36 = (t46 + t103) * qJD(3);
t29 = t138 * t70;
t23 = t55 * t53;
t19 = -t49 ^ 2 + t51 ^ 2;
t15 = t51 * t70 - t28;
t14 = -t116 + (-t102 + t49) * t70;
t8 = t30 * t78 - t130;
t7 = -t30 * t75 - t125;
t5 = t26 * t78 - t130;
t4 = t24 * t70 - t55 * t47;
t3 = -t138 * t53 - t47 * t54;
t1 = t96 + t136;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * qJD(3), -t108, 0, t43 * t53 - t45 * t48 + t46 * t47 - t132, 0, 0, 0, 0, 0, 0, -t47 * t110 - t53 * t121 + (-t110 * t137 - t48 * t79) * qJD(3), -t47 * t109 + t53 * t122 + (-t109 * t137 + t48 * t76) * qJD(3), t114 * t108, t41 * t48 - t47 * t89 + t53 * t84 - t132, 0, 0, 0, 0, 0, 0, -t137 * t28 + t4 * t70 + t48 * t49, t137 * t27 - t3 * t70 + t48 * t51, -t23 * t27 + t24 * t28 - t3 * t49 - t4 * t51, -t1 * t24 - t137 * t36 - t2 * t23 + t3 * t6 + t39 * t48 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121, 0, -qJD(4) * t89 + t16 * t76 + t17 * t79, 0, 0, 0, 0, 0, 0, -t29, t134, t90 + t118, t1 * t55 - t138 * t5 - t2 * t54 - t32 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t93, -0.2e1 * t115 * t106, t121, -0.2e1 * t93, -t122, 0, t76 * t85 - t79 * t86, t76 * t86 + t79 * t85, -qJD(3) * t114 * t45 + t84, -t44 * pkin(3) + pkin(6) * t84 - t41 * t46 + t45 * t89, -t27 * t55 - t32 * t51, -t90 + t118, -t134, t138 * t49 + t28 * t54, -t29, 0, t119 * t70 + t138 * t39 + t68 * t28 + t36 * t54 + t49 * t92, -t120 * t70 - t68 * t27 - t39 * t32 + t36 * t55 + t51 * t92, -t1 * t54 - t119 * t51 - t120 * t49 - t138 * t6 - t2 * t55 + t37 * t27 - t38 * t28 + t5 * t32, t1 * t38 + t119 * t5 + t120 * t6 + t2 * t37 + t36 * t68 + t39 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t115 * t82, 0, t105, 0, 0, (-qJD(3) * t41 - t43) * t76, -t41 * t111 + (t34 + t128) * qJD(4) - t117, 0, 0, t131, t19, t14, -t131, t15, 0, -t49 * t104 - t133 - t7 * t70 + (t100 * t75 - t125) * qJD(5) + t99, -t51 * t104 + t8 * t70 + (qJD(5) * t100 - t10) * t78 + t87, (t6 + t7) * t51 + (-t5 + t8) * t49 + (t27 * t78 - t28 * t75 + (-t49 * t78 + t51 * t75) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (-t39 * t112 + t1 * t75 + t2 * t78 + (-t5 * t75 + t6 * t78) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t19, t14, -t131, t15, 0, t6 * t70 - t133 + t2, t5 * t70 - t136 + t87, 0, 0;];
tauc_reg = t9;
