% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:42
% EndTime: 2019-03-09 01:37:45
% DurationCPUTime: 1.13s
% Computational Cost: add. (871->186), mult. (1628->293), div. (0->0), fcn. (944->6), ass. (0->109)
t55 = cos(qJ(5));
t93 = t55 * qJD(1);
t37 = -qJD(6) + t93;
t126 = qJD(6) + t37;
t48 = sin(pkin(9));
t49 = cos(pkin(9));
t32 = t49 * qJD(2) + t48 * qJD(3);
t23 = t32 * qJD(1);
t125 = 0.2e1 * t23;
t92 = qJD(1) * qJD(5);
t124 = -0.2e1 * t92;
t33 = t48 * qJD(2) - t49 * qJD(3);
t24 = t33 * qJD(1);
t53 = sin(qJ(5));
t46 = t53 ^ 2;
t70 = qJD(1) * t46 - t37 * t55;
t52 = sin(qJ(6));
t98 = qJD(6) * t52;
t89 = t53 * t98;
t54 = cos(qJ(6));
t94 = t54 * qJD(5);
t123 = -t37 * t89 - t70 * t94;
t97 = qJD(6) * t54;
t87 = t53 * t97;
t91 = qJD(5) * qJD(6);
t95 = t52 * qJD(5);
t17 = (t55 * t95 + t87) * qJD(1) + t52 * t91;
t122 = 0.2e1 * qJD(1);
t42 = qJD(1) * qJ(2) + qJD(3);
t35 = pkin(3) * qJD(1) + t42;
t51 = -pkin(1) - qJ(3);
t36 = t51 * qJD(1) + qJD(2);
t15 = t48 * t35 + t49 * t36;
t12 = qJD(1) * pkin(7) + t15;
t9 = t53 * qJD(4) + t55 * t12;
t4 = t9 * qJD(5) + t53 * t24;
t121 = t4 * t52;
t120 = t4 * t54;
t82 = t55 * t92;
t102 = qJD(1) * t53;
t86 = t52 * t102;
t16 = -qJD(6) * t86 + (t82 + t91) * t54;
t119 = t16 * t52;
t25 = t86 - t94;
t118 = t25 * t37;
t27 = t54 * t102 + t95;
t117 = t27 * t37;
t57 = qJD(1) ^ 2;
t116 = t48 * t57;
t115 = t49 * t57;
t114 = t52 * t37;
t113 = t52 * t55;
t112 = t53 * t16;
t111 = t53 * t25;
t110 = t54 * t37;
t109 = t54 * t55;
t108 = t55 * t17;
t56 = qJD(5) ^ 2;
t107 = t56 * t53;
t106 = t56 * t55;
t50 = pkin(3) + qJ(2);
t105 = t48 * t50 + t49 * t51;
t104 = -t55 ^ 2 + t46;
t103 = t56 + t57;
t22 = pkin(7) + t105;
t101 = qJD(5) * t22;
t100 = qJD(5) * t53;
t99 = qJD(5) * t55;
t96 = qJD(6) * t55;
t90 = t37 * t98;
t88 = t37 * t97;
t44 = qJD(2) * t122;
t6 = qJD(5) * pkin(8) + t9;
t85 = t22 * t37 + t6;
t84 = t103 * t55;
t83 = t103 * t53;
t80 = t53 * t92;
t79 = t27 * t100 - t16 * t55;
t14 = t49 * t35 - t48 * t36;
t78 = -t48 * t51 + t49 * t50;
t11 = -qJD(1) * pkin(4) - t14;
t77 = -qJD(1) * t11 - t24;
t75 = t37 * t87;
t74 = t49 * t124;
t73 = 0.2e1 * t80;
t72 = pkin(5) * t53 - pkin(8) * t55;
t69 = -t55 * pkin(5) - t53 * pkin(8) - pkin(4);
t7 = t69 * qJD(1) - t14;
t2 = t52 * t7 + t54 * t6;
t71 = t52 * t6 - t54 * t7;
t8 = t55 * qJD(4) - t53 * t12;
t68 = t49 * t109 + t52 * t48;
t67 = t49 * t113 - t54 * t48;
t66 = t48 * t109 - t52 * t49;
t65 = -t48 * t113 - t54 * t49;
t64 = -t22 * t56 + t125;
t63 = qJD(5) * (qJD(1) * (-pkin(4) - t78) + t11 - t33);
t62 = t70 * t52;
t3 = t8 * qJD(5) + t55 * t24;
t5 = -qJD(5) * pkin(5) - t8;
t60 = t5 * qJD(5) + qJD(6) * t7 + t33 * t37 + t3;
t19 = t72 * qJD(5) - t32;
t59 = (-t52 * t96 - t53 * t94) * t37 + t27 * t99 + t112;
t58 = -(t53 * t95 - t54 * t96) * t37 + t25 * t99 + t53 * t17;
t31 = t72 * qJD(1);
t18 = t69 - t78;
t13 = t19 * qJD(1);
t10 = t54 * t13;
t1 = [0, 0, 0, 0, t44, qJ(2) * t44, t44, qJD(3) * t122, t42 * qJD(2) - t36 * qJD(3) + (qJ(2) * qJD(2) - qJD(3) * t51) * qJD(1), t125, -0.2e1 * t24, t105 * t24 + t14 * t32 + t15 * t33 + t23 * t78, t55 * t73, t104 * t124, t106, -t107, 0, t53 * t63 + t55 * t64, -t53 * t64 + t55 * t63, t54 * t112 + (t55 * t94 - t89) * t27 (-t25 * t54 - t27 * t52) * t99 + (-t119 - t17 * t54 + (t25 * t52 - t27 * t54) * qJD(6)) * t53, -t123 + t79, t75 + t108 + (-t62 - t111) * qJD(5) (-t37 - t93) * t100 -(-t18 * t98 + t54 * t19) * t37 + (t101 * t25 + t52 * t60 + t85 * t97 - t10) * t55 + (t5 * t97 + t22 * t17 + t33 * t25 + t121 + (-t22 * t114 + (-t113 * t22 + t54 * t18) * qJD(1) - t71) * qJD(5)) * t53 (t18 * t97 + t52 * t19) * t37 + (t27 * t101 + (-qJD(6) * t85 + t13) * t52 + t60 * t54) * t55 + (-t5 * t98 + t22 * t16 + t33 * t27 + t120 + (-t22 * t110 - (t109 * t22 + t52 * t18) * qJD(1) - t2) * qJD(5)) * t53; 0, 0, 0, 0, -t57, -t57 * qJ(2), -t57, 0 (-qJD(3) - t42) * qJD(1), -t115, t116, -t23 * t48 + t24 * t49 + (-t14 * t49 - t15 * t48) * qJD(1), 0, 0, 0, 0, 0, t48 * t73 - t49 * t84, 0.2e1 * t48 * t82 + t49 * t83, 0, 0, 0, 0, 0, t48 * t90 + t58 * t49 + (t65 * t37 + (-qJD(5) * t67 - t48 * t25) * t53) * qJD(1), t48 * t88 + t59 * t49 + (-t66 * t37 + (-qJD(5) * t68 - t48 * t27) * t53) * qJD(1); 0, 0, 0, 0, 0, 0, 0, -t57 (qJD(2) + t36) * qJD(1), -t116, -t115, t23 * t49 + t24 * t48 + (-t14 * t48 + t15 * t49) * qJD(1), 0, 0, 0, 0, 0, -t48 * t84 + t53 * t74, t48 * t83 + t55 * t74, 0, 0, 0, 0, 0, -t49 * t90 + t58 * t48 + (t67 * t37 + (qJD(5) * t65 + t49 * t25) * t53) * qJD(1), -t49 * t88 + t59 * t48 + (t68 * t37 + (-qJD(5) * t66 + t49 * t27) * t53) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106, 0, 0, 0, 0, 0, t75 - t108 + (-t62 + t111) * qJD(5), t123 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t57 * t55, t104 * t57, 0, 0, 0, t77 * t53, t77 * t55, -t110 * t27 + t119 (t16 + t118) * t54 + (-t17 + t117) * t52, -t88 + (t37 * t109 + (-t27 + t95) * t53) * qJD(1), t90 + (-t37 * t113 + (t25 + t94) * t53) * qJD(1), t37 * t102, -pkin(5) * t17 - t120 + (t54 * t31 - t52 * t8) * t37 - t9 * t25 + (pkin(8) * t110 + t5 * t52) * qJD(6) + (t71 * t53 + (-pkin(8) * t100 - t5 * t55) * t52) * qJD(1), -pkin(5) * t16 + t121 - (t52 * t31 + t54 * t8) * t37 - t9 * t27 + (-pkin(8) * t114 + t5 * t54) * qJD(6) + (-t5 * t109 + (-pkin(8) * t94 + t2) * t53) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t25, -t25 ^ 2 + t27 ^ 2, t16 - t118, -t117 - t17, t80, -t126 * t2 - t5 * t27 - t52 * t3 + t10, t126 * t71 - t52 * t13 + t5 * t25 - t54 * t3;];
tauc_reg  = t1;
