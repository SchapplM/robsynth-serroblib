% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:52
% EndTime: 2019-03-09 01:35:55
% DurationCPUTime: 1.01s
% Computational Cost: add. (796->188), mult. (1526->280), div. (0->0), fcn. (876->6), ass. (0->118)
t50 = sin(qJ(5));
t95 = t50 * qJD(1);
t33 = -qJD(6) + t95;
t132 = qJD(6) + t33;
t52 = cos(qJ(5));
t45 = t52 ^ 2;
t49 = sin(qJ(6));
t131 = (qJD(1) * t45 - t33 * t50) * t49;
t53 = -pkin(1) - pkin(2);
t46 = sin(pkin(9));
t92 = qJD(1) * qJD(2);
t78 = t46 * t92;
t31 = t53 * qJD(1) + qJD(2);
t47 = cos(pkin(9));
t93 = qJD(1) * qJ(2);
t18 = t47 * t31 - t46 * t93;
t15 = qJD(1) * pkin(3) + qJD(4) - t18;
t11 = qJD(1) * pkin(7) + t15;
t8 = t52 * qJD(3) + t50 * t11;
t4 = t8 * qJD(5) - t52 * t78;
t130 = t4 * t49;
t51 = cos(qJ(6));
t129 = t4 * t51;
t99 = qJD(6) * t52;
t79 = qJD(1) * t99;
t91 = qJD(1) * qJD(5);
t80 = t51 * t91;
t90 = qJD(5) * qJD(6);
t13 = t49 * t79 + t50 * t80 + t51 * t90;
t128 = t13 * t49;
t127 = t13 * t50;
t81 = t50 * t91;
t14 = t51 * t79 + (-t81 - t90) * t49;
t126 = t14 * t50;
t106 = qJD(1) * t52;
t94 = t51 * qJD(5);
t23 = t49 * t106 + t94;
t125 = t23 * t33;
t85 = t51 * t106;
t96 = t49 * qJD(5);
t24 = t85 - t96;
t124 = t24 * t33;
t123 = t46 * t31;
t55 = qJD(1) ^ 2;
t122 = t46 * t55;
t121 = t47 * t55;
t120 = t49 * t33;
t119 = t49 * t50;
t118 = t50 * t51;
t117 = t51 * t33;
t116 = t51 * t46;
t115 = t52 * t13;
t114 = t52 * t14;
t113 = t52 * t23;
t112 = t52 * t24;
t54 = qJD(5) ^ 2;
t111 = t54 * t52;
t110 = t55 * t50;
t39 = t47 * qJ(2);
t109 = t46 * t53 + t39;
t108 = t50 ^ 2 - t45;
t107 = t54 + t55;
t75 = -t46 * qJ(2) + t47 * t53;
t70 = pkin(3) - t75;
t21 = pkin(7) + t70;
t105 = qJD(5) * t21;
t104 = qJD(5) * t50;
t103 = qJD(5) * t52;
t102 = qJD(6) * t49;
t101 = qJD(6) * t50;
t100 = qJD(6) * t51;
t98 = t46 * qJD(2);
t97 = t47 * qJD(2);
t89 = t50 * t117;
t22 = -qJ(4) + t109;
t88 = t33 * t102;
t87 = t49 * t99;
t86 = t33 * t100;
t84 = 0.2e1 * t92;
t6 = qJD(5) * pkin(8) + t8;
t83 = t21 * t33 + t6;
t82 = t47 * t107;
t77 = t47 * t92;
t76 = t52 * t91;
t19 = t47 * t93 + t123;
t43 = qJD(1) * qJ(4);
t16 = t19 - t43;
t74 = t16 + t98;
t73 = t52 * t86;
t72 = -0.2e1 * t76;
t71 = -qJD(1) + t101;
t69 = -pkin(5) * t52 - pkin(8) * t50;
t68 = -t50 * pkin(5) + t52 * pkin(8);
t9 = t123 - t43 + (t68 + t39) * qJD(1);
t2 = t49 * t9 + t51 * t6;
t67 = t49 * t6 - t51 * t9;
t66 = t18 * t46 - t19 * t47;
t65 = (t33 + t95) * t52;
t63 = t50 * qJD(3) - t52 * t11;
t62 = t74 * qJD(1);
t61 = -t33 * t87 + t45 * t80;
t5 = -qJD(5) * pkin(5) + t63;
t60 = pkin(8) * t103 - t5 * t50;
t42 = qJD(1) * qJD(4);
t26 = -t42 + t77;
t32 = -qJD(4) + t97;
t59 = -t32 * qJD(1) - t21 * t54 - t26;
t58 = qJD(1) * t22 + t16 - t98;
t3 = -t63 * qJD(5) + t50 * t78;
t57 = t5 * qJD(5) + qJD(6) * t9 + t33 * t98 + t3;
t56 = t69 * qJD(5) + t97;
t41 = t54 * t50;
t25 = t69 * qJD(1);
t20 = -qJD(4) + t56;
t17 = t68 + t22;
t12 = t56 * qJD(1) - t42;
t10 = t51 * t12;
t1 = [0, 0, 0, 0, t84, qJ(2) * t84, 0.2e1 * t78, 0.2e1 * t77 ((t47 * t109 - t46 * t75) * qJD(1) - t66) * qJD(2), -0.2e1 * t78, t42 + (-t32 - t97) * qJD(1), t16 * t32 + t26 * t22 + (qJD(1) * t70 + t15) * t98, t50 * t72, 0.2e1 * t108 * t91, t41, t111, 0, -t58 * t103 + t59 * t50, t58 * t104 + t59 * t52, -t51 * t115 + (-t50 * t94 - t87) * t24 (t23 * t51 + t24 * t49) * t104 + (t128 - t14 * t51 + (t23 * t49 - t24 * t51) * qJD(6)) * t52, -t127 + (-t89 + t112) * qJD(5) + t61, -t73 - t126 + (-t113 - t131) * qJD(5), qJD(5) * t65 -(-t17 * t102 + t51 * t20) * t33 + (t83 * t100 - t23 * t105 + t57 * t49 - t10) * t50 + (t23 * t98 - t5 * t100 + t21 * t14 - t130 + (t21 * t120 - (-t21 * t119 + t51 * t17) * qJD(1) + t67) * qJD(5)) * t52 (t17 * t100 + t49 * t20) * t33 + (-t24 * t105 + (-t83 * qJD(6) + t12) * t49 + t57 * t51) * t50 + (t24 * t98 + t5 * t102 - t21 * t13 - t129 + (t21 * t117 + (t21 * t118 + t49 * t17) * qJD(1) + t2) * qJD(5)) * t52; 0, 0, 0, 0, -t55, -t55 * qJ(2), -t122, -t121, t66 * qJD(1), t122, t121, t26 * t46 + (-t16 * t47 + (-t15 - t97) * t46) * qJD(1), 0, 0, 0, 0, 0, t46 * t72 + t50 * t82, 0.2e1 * t46 * t81 + t52 * t82, 0, 0, 0, 0, 0, t46 * t88 + (-(t50 * t100 + t52 * t96) * t33 + t23 * t104 - t114) * t47 + ((-t46 * t119 + t47 * t51) * t33 + (-(t47 * t119 + t116) * qJD(5) - t46 * t23) * t52) * qJD(1), t46 * t86 + ((t49 * t101 - t52 * t94) * t33 + t24 * t104 + t115) * t47 + (-(t50 * t116 + t47 * t49) * t33 + ((-t47 * t118 + t49 * t46) * qJD(5) - t46 * t24) * t52) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t41, 0, 0, 0, 0, 0, t73 - t126 + (-t113 + t131) * qJD(5), t127 + (-t89 - t112) * qJD(5) + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t62, 0, 0, 0, 0, 0, -t41 - t110, -t107 * t52, 0, 0, 0, 0, 0, t114 + t71 * t117 + (-t50 * t23 + t49 * t65) * qJD(5), -t115 - t71 * t120 + (t52 * t117 + (-t24 + t85) * t50) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t110, -t108 * t55, 0, 0, 0, t52 * t62, -t74 * t95, t24 * t117 + t128 (t13 - t125) * t51 + (t14 - t124) * t49, -t86 + (t89 + (-t24 - t96) * t52) * qJD(1), t88 + (-t33 * t119 + (t23 - t94) * t52) * qJD(1), -t33 * t106, pkin(5) * t14 - t129 + (t51 * t25 + t49 * t63) * t33 + t8 * t23 + (pkin(8) * t117 + t5 * t49) * qJD(6) + (t60 * t49 - t52 * t67) * qJD(1), -pkin(5) * t13 + t130 - (t49 * t25 - t51 * t63) * t33 + t8 * t24 + (-pkin(8) * t120 + t5 * t51) * qJD(6) + (-t2 * t52 + t60 * t51) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t23, -t23 ^ 2 + t24 ^ 2, t13 + t125, t14 + t124, -t76, -t132 * t2 + t5 * t24 - t49 * t3 + t10, -t49 * t12 + t132 * t67 - t5 * t23 - t51 * t3;];
tauc_reg  = t1;
