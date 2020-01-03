% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:52
% DurationCPUTime: 1.19s
% Computational Cost: add. (726->210), mult. (1719->295), div. (0->0), fcn. (881->4), ass. (0->120)
t111 = qJ(4) * qJD(1);
t71 = sin(qJ(2));
t117 = t71 * qJD(1);
t52 = pkin(6) * t117;
t29 = -t71 * t111 + t52;
t101 = qJD(3) + t29;
t74 = -pkin(2) - pkin(3);
t99 = qJD(2) * t74;
t18 = t99 + t101;
t45 = qJD(5) + t117;
t142 = qJD(5) - t45;
t73 = cos(qJ(2));
t124 = qJD(1) * t73;
t72 = cos(qJ(5));
t115 = t72 * qJD(2);
t70 = sin(qJ(5));
t27 = t70 * t124 - t115;
t109 = qJD(1) * qJD(2);
t141 = -0.2e1 * t109;
t129 = pkin(6) - qJ(4);
t53 = pkin(6) * t124;
t31 = -t73 * t111 + t53;
t65 = qJD(2) * qJ(3);
t24 = -t31 - t65;
t34 = -t73 * pkin(2) - t71 * qJ(3) - pkin(1);
t26 = t73 * pkin(3) - t34;
t87 = t71 * pkin(4) + t73 * pkin(7);
t12 = t87 + t26;
t110 = qJ(4) * qJD(2);
t100 = t73 * t110;
t116 = t71 * qJD(4);
t97 = t73 * t109;
t44 = pkin(6) * t97;
t16 = t44 + (-t100 - t116) * qJD(1);
t19 = qJD(2) * pkin(4) - t24;
t122 = qJD(2) * t73;
t22 = t129 * t122 - t116;
t112 = qJ(3) * qJD(1);
t25 = -qJD(1) * pkin(1) - pkin(2) * t124 - t71 * t112;
t15 = pkin(3) * t124 + qJD(4) - t25;
t6 = t87 * qJD(1) + t15;
t140 = -(qJD(5) * t12 + t22) * t45 + (t19 * qJD(2) - qJD(5) * t6 - t16) * t71;
t98 = t71 * t109;
t42 = qJ(4) * t98;
t63 = qJD(2) * qJD(3);
t123 = qJD(2) * t71;
t81 = pkin(6) * t123 + t73 * qJD(4);
t9 = t81 * qJD(1) - t42 - t63;
t139 = t9 * t70;
t138 = t9 * t72;
t10 = t27 * qJD(5) + t72 * t98;
t137 = t10 * t70;
t136 = t27 * t45;
t118 = t70 * qJD(2);
t28 = t72 * t124 + t118;
t135 = t28 * t45;
t134 = t45 * t71;
t133 = t45 * t72;
t76 = qJD(1) ^ 2;
t132 = t73 * t76;
t75 = qJD(2) ^ 2;
t131 = t75 * t71;
t130 = t75 * t73;
t56 = t71 * qJD(3);
t128 = qJ(3) * t97 + qJD(1) * t56;
t127 = t73 * t65 + t56;
t66 = t71 ^ 2;
t67 = t73 ^ 2;
t126 = t66 - t67;
t125 = qJD(2) * pkin(2);
t64 = -pkin(7) + t74;
t13 = t64 * qJD(2) + t101;
t121 = qJD(5) * t13;
t120 = qJD(5) * t70;
t119 = t19 * qJD(5);
t113 = -qJD(4) - t15;
t108 = t70 * t134;
t107 = t71 * t133;
t106 = t71 * t132;
t105 = t45 * t120;
t104 = t73 * t120;
t103 = qJD(5) * t133;
t88 = t71 * t99;
t14 = t88 + t127;
t8 = qJD(1) * t88 + t128;
t96 = qJD(1) * t14 + t8;
t95 = t113 * t71;
t94 = qJD(1) * t26 + t15;
t93 = qJD(1) * t34 + t25;
t91 = pkin(1) * t141;
t90 = qJD(3) - t125;
t89 = 0.2e1 * t97;
t2 = t72 * t13 + t70 * t6;
t86 = t70 * t13 - t72 * t6;
t36 = t129 * t71;
t80 = pkin(4) * t73 + t64 * t71;
t78 = t80 * qJD(2);
t85 = (-qJD(5) * t36 + t127 + t78) * t45;
t84 = qJD(1) * t67 - t134;
t17 = pkin(2) * t98 - t128;
t23 = pkin(2) * t123 - t127;
t83 = -pkin(6) * t75 - qJD(1) * t23 - t17;
t82 = -t64 * t122 - t19 * t71;
t11 = t28 * qJD(5) - t70 * t98;
t32 = -pkin(6) * t98 + t63;
t33 = t52 + t90;
t35 = t53 + t65;
t77 = t32 * t73 + (t33 * t73 + (-t35 + t53) * t71) * qJD(2);
t69 = qJ(3) + pkin(4);
t60 = 0.2e1 * t63;
t49 = t73 * t112;
t41 = -t66 * t76 - t75;
t37 = t129 * t73;
t30 = pkin(2) * t117 - t49;
t21 = -t71 * t110 + t81;
t20 = t74 * t117 + t49;
t7 = t80 * qJD(1) + t49;
t4 = qJD(1) * t78 + t128;
t3 = t72 * t4;
t1 = [0, 0, 0, t71 * t89, t126 * t141, t130, -t131, 0, -pkin(6) * t130 + t71 * t91, pkin(6) * t131 + t73 * t91, t93 * t123 + t83 * t73, t77, -t93 * t122 + t83 * t71, t77 * pkin(6) + t17 * t34 + t25 * t23, t96 * t71 + (t94 * t73 - t21) * qJD(2), -t96 * t73 + (t94 * t71 + t22) * qJD(2), -t16 * t71 + t9 * t73 + (-t18 * t73 - t24 * t71) * qJD(2) + (t21 * t73 - t22 * t71 + (-t36 * t73 + t37 * t71) * qJD(2)) * qJD(1), t15 * t14 + t16 * t36 + t18 * t22 + t24 * t21 + t8 * t26 - t9 * t37, -t10 * t72 * t73 + (-t115 * t71 - t104) * t28, (t27 * t72 + t28 * t70) * t123 + (t137 - t11 * t72 + (t27 * t70 - t28 * t72) * qJD(5)) * t73, t45 * t104 + t10 * t71 + (-t28 * t73 - t84 * t72) * qJD(2), t73 * t103 + t11 * t71 + (t27 * t73 + t70 * t84) * qJD(2), (t45 + t117) * t122, -t37 * t11 + t21 * t27 + t3 * t71 + (-t121 * t71 + t85) * t72 + t140 * t70 + (-t72 * t119 + t139 + ((t72 * t12 - t70 * t36) * qJD(1) - t86) * qJD(2)) * t73, t37 * t10 + t21 * t28 + (-t85 - (t4 - t121) * t71) * t70 + t140 * t72 + (t70 * t119 + t138 + (-(t70 * t12 + t72 * t36) * qJD(1) - t2) * qJD(2)) * t73; 0, 0, 0, -t106, t126 * t76, 0, 0, 0, t76 * pkin(1) * t71, pkin(1) * t132, (-t25 * t71 + t30 * t73) * qJD(1), ((t35 - t65) * t71 + (-t33 + t90) * t73) * qJD(1), t60 + (t25 * t73 + t30 * t71) * qJD(1), t32 * qJ(3) + t35 * qJD(3) - t25 * t30 + (t35 * t71 + (-t33 - t125) * t73) * qJD(1) * pkin(6), t29 * qJD(2) + t42 + t60 + (t113 * t73 + (-pkin(6) * qJD(2) - t20) * t71) * qJD(1), -t31 * qJD(2) + t44 + ((t20 - t110) * t73 + t95) * qJD(1), 0, -t9 * qJ(3) - t101 * t24 - t15 * t20 + t16 * t74 - t18 * t31, t133 * t28 - t137, (-t10 - t136) * t72 + (-t11 - t135) * t70, -t103 + (-t107 + (t28 - t118) * t73) * qJD(1), t105 + (t108 + (-t27 - t115) * t73) * qJD(1), -t45 * t124, -t69 * t11 - t138 - (-t70 * t31 + t72 * t7) * t45 - t101 * t27 + (-t133 * t64 - t19 * t70) * qJD(5) + (t70 * t82 + t73 * t86) * qJD(1), t69 * t10 + t139 + (t72 * t31 + t70 * t7) * t45 - t101 * t28 + (t70 * t64 * t45 - t19 * t72) * qJD(5) + (t2 * t73 + t82 * t72) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, t41, -t35 * qJD(2) + t25 * t117 + t44, t41, t106, 0, t24 * qJD(2) + t44 + (t95 - t100) * qJD(1), 0, 0, 0, 0, 0, -t103 + qJD(2) * t27 + (-t118 * t73 - t107) * qJD(1), t105 + qJD(2) * t28 + (-t115 * t73 + t108) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0.2e1 * t98, (-t66 - t67) * t76, (-t24 * t73 + (t18 + t99) * t71) * qJD(1) + t128, 0, 0, 0, 0, 0, -t105 + (-t108 + (-t27 + t115) * t73) * qJD(1), -t103 + (-t107 + (-t28 - t118) * t73) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t27, -t27 ^ 2 + t28 ^ 2, t10 - t136, t11 - t135, t97, -t142 * t2 - t70 * t16 + t19 * t28 + t3, t142 * t86 - t72 * t16 - t19 * t27 - t70 * t4;];
tauc_reg = t1;
