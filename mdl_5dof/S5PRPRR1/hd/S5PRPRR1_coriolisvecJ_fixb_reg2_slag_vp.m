% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:05
% EndTime: 2019-12-05 15:43:09
% DurationCPUTime: 1.02s
% Computational Cost: add. (2009->176), mult. (5147->245), div. (0->0), fcn. (3959->6), ass. (0->106)
t140 = cos(qJ(5));
t100 = cos(qJ(4));
t97 = cos(pkin(9));
t124 = qJD(2) * t97;
t116 = t100 * t124;
t96 = sin(pkin(9));
t99 = sin(qJ(4));
t130 = t99 * t96;
t118 = qJD(2) * t130;
t67 = -t116 + t118;
t78 = t100 * t96 + t99 * t97;
t69 = t78 * qJD(2);
t98 = sin(qJ(5));
t107 = -t140 * t69 + t98 * t67;
t85 = qJD(4) * t116;
t59 = qJD(4) * t118 - t85;
t72 = t78 * qJD(4);
t60 = qJD(2) * t72;
t103 = t107 * qJD(5) - t140 * t60 + t98 * t59;
t95 = qJD(4) + qJD(5);
t138 = t107 * t95;
t147 = t103 - t138;
t112 = qJD(5) * t140;
t122 = qJD(5) * t98;
t106 = -t67 * t112 - t69 * t122 - t140 * t59 - t98 * t60;
t38 = -t140 * t67 - t98 * t69;
t136 = t38 * t95;
t146 = t106 - t136;
t135 = t107 ^ 2;
t137 = t38 ^ 2;
t145 = t135 - t137;
t134 = t38 * t107;
t90 = -t97 * pkin(3) - pkin(2);
t81 = t90 * qJD(2) + qJD(3);
t44 = t67 * pkin(4) + t81;
t144 = t44 * t107;
t123 = qJD(3) * t96;
t115 = qJD(2) * t123;
t119 = qJD(4) * t100;
t129 = pkin(6) + qJ(3);
t82 = t129 * t96;
t92 = t97 * qJD(1);
t62 = -qJD(2) * t82 + t92;
t125 = t100 * t97;
t86 = qJD(3) * t125;
t127 = qJD(2) * t86 + t62 * t119;
t120 = qJ(3) * qJD(2);
t80 = t96 * qJD(1) + t97 * t120;
t63 = pkin(6) * t124 + t80;
t22 = (-qJD(4) * t63 - t115) * t99 + t127;
t16 = -t60 * pkin(7) + t22;
t105 = t78 * qJD(3);
t104 = qJD(2) * t105;
t32 = t100 * t63 + t99 * t62;
t23 = -t32 * qJD(4) - t104;
t17 = t59 * pkin(7) + t23;
t131 = t99 * t63;
t31 = t100 * t62 - t131;
t27 = -t69 * pkin(7) + t31;
t24 = qJD(4) * pkin(4) + t27;
t28 = -t67 * pkin(7) + t32;
t102 = -t24 * t112 + t28 * t122 - t140 * t16 - t98 * t17;
t143 = -t44 * t38 + t102;
t83 = t129 * t97;
t46 = t100 * t83 - t99 * t82;
t142 = t69 ^ 2;
t71 = qJD(4) * t130 - t97 * t119;
t77 = t125 - t130;
t19 = -t77 * t112 + t78 * t122 + t140 * t71 + t98 * t72;
t41 = t140 * t78 + t98 * t77;
t141 = t103 * t41 - t19 * t38;
t139 = t19 * t95;
t133 = t69 * t67;
t132 = t98 * t28;
t128 = -t78 * t60 + t71 * t67;
t126 = t96 ^ 2 + t97 ^ 2;
t121 = t71 * qJD(4);
t117 = t140 * t28;
t114 = t140 * t17 - t98 * t16;
t45 = -t100 * t82 - t99 * t83;
t111 = t126 * qJD(2);
t20 = t41 * qJD(5) + t140 * t72 - t98 * t71;
t40 = -t140 * t77 + t98 * t78;
t110 = t106 * t40 - t107 * t20;
t109 = -t77 * t59 - t69 * t72;
t108 = (-t96 * t120 + t92) * t96 - t80 * t97;
t33 = -t78 * pkin(7) + t45;
t34 = t77 * pkin(7) + t46;
t11 = t140 * t33 - t98 * t34;
t6 = t98 * t24 + t117;
t12 = t140 * t34 + t98 * t33;
t29 = -t82 * t119 + t86 + (-qJD(4) * t83 - t123) * t99;
t2 = -t6 * qJD(5) + t114;
t30 = -t46 * qJD(4) - t105;
t65 = t67 ^ 2;
t61 = t72 * qJD(4);
t52 = -t77 * pkin(4) + t90;
t26 = t71 * pkin(7) + t30;
t25 = -t72 * pkin(7) + t29;
t18 = t20 * t95;
t8 = t140 * t27 - t132;
t7 = -t98 * t27 - t117;
t5 = t140 * t24 - t132;
t4 = -t12 * qJD(5) + t140 * t26 - t98 * t25;
t3 = t11 * qJD(5) + t140 * t25 + t98 * t26;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t121, -t109 + t128, t22 * t78 + t23 * t77 - t31 * t72 - t32 * t71, 0, 0, 0, 0, 0, 0, -t18, t139, t110 + t141, -t102 * t41 - t6 * t19 - t2 * t40 - t5 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t111, (qJ(3) * t111 - t108) * qJD(3), -t59 * t78 - t69 * t71, t109 + t128, -t121, -t60 * t77 + t67 * t72, -t61, 0, t30 * qJD(4) + t90 * t60 + t81 * t72, -t29 * qJD(4) - t90 * t59 - t81 * t71, t22 * t77 - t23 * t78 - t29 * t67 - t30 * t69 + t31 * t71 - t32 * t72 + t45 * t59 - t46 * t60, t22 * t46 + t23 * t45 + t32 * t29 + t31 * t30, t106 * t41 + t107 * t19, -t110 + t141, -t139, -t103 * t40 - t20 * t38, -t18, 0, -t52 * t103 + t44 * t20 + t4 * t95 + (-t38 * t72 + t40 * t60) * pkin(4), t52 * t106 - t44 * t19 - t3 * t95 + (-t107 * t72 + t41 * t60) * pkin(4), t102 * t40 + t103 * t12 - t106 * t11 + t107 * t4 + t5 * t19 - t2 * t41 - t6 * t20 + t3 * t38, -t102 * t12 + t2 * t11 + t6 * t3 + t5 * t4 + (t44 * t72 + t52 * t60) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * qJD(2) ^ 2, t108 * qJD(2), 0, 0, 0, 0, 0, 0, 0.2e1 * t69 * qJD(4), t85 + (-t67 - t118) * qJD(4), -t65 - t142, t31 * t69 + t32 * t67, 0, 0, 0, 0, 0, 0, -t103 - t138, t106 + t136, -t135 - t137, t60 * pkin(4) - t107 * t5 - t6 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, -t65 + t142, t85 + (t67 - t118) * qJD(4), -t133, 0, 0, -t81 * t69 - t104, t99 * t115 + t81 * t67 + (t31 + t131) * qJD(4) - t127, 0, 0, t134, t145, t146, -t134, t147, 0, t69 * pkin(4) * t38 + t144 - t7 * t95 + (-t117 + (-pkin(4) * t95 - t24) * t98) * qJD(5) + t114, t8 * t95 + (t107 * t69 - t95 * t112) * pkin(4) + t143, -t6 * t107 - t8 * t38 + t5 * t38 - t7 * t107 + (-t140 * t106 + t103 * t98 + (-t107 * t98 + t140 * t38) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t140 * t2 - t102 * t98 - t44 * t69 + (t140 * t6 - t5 * t98) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t145, t146, -t134, t147, 0, t6 * t95 + t144 + t2, t5 * t95 + t143, 0, 0;];
tauc_reg = t1;
