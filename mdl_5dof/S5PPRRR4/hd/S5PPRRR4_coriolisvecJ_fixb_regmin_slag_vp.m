% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:43
% EndTime: 2019-12-05 15:19:48
% DurationCPUTime: 1.14s
% Computational Cost: add. (836->175), mult. (2378->297), div. (0->0), fcn. (2038->12), ass. (0->115)
t52 = sin(pkin(6));
t62 = cos(qJ(3));
t123 = t52 * t62;
t54 = cos(pkin(11));
t55 = cos(pkin(6));
t121 = t54 * t55;
t51 = sin(pkin(11));
t53 = sin(pkin(5));
t59 = sin(qJ(3));
t133 = (t62 * t121 - t51 * t59) * t53;
t56 = cos(pkin(5));
t135 = t56 * t123 + t133;
t45 = t56 * qJD(1) + qJD(2);
t134 = qJD(1) * t133 + t45 * t123;
t111 = qJD(3) * t62;
t112 = qJD(3) * t59;
t114 = qJD(1) * t53;
t95 = t54 * t114;
t84 = t55 * t95;
t94 = t52 * t112;
t96 = t51 * t114;
t17 = t96 * t111 + t84 * t112 + t45 * t94;
t124 = t52 * t59;
t68 = (t59 * t121 + t51 * t62) * t53;
t21 = qJD(1) * t68 + t45 * t124;
t132 = t21 * qJD(3) - t17;
t61 = cos(qJ(4));
t103 = t61 * qJD(3);
t46 = -qJD(5) + t103;
t131 = -qJD(5) - t46;
t20 = -t59 * t96 + (t45 * t52 + t84) * t62;
t101 = qJD(4) * qJD(5);
t57 = sin(qJ(5));
t105 = t57 * qJD(4);
t60 = cos(qJ(5));
t107 = qJD(5) * t60;
t58 = sin(qJ(4));
t90 = t58 * t107;
t29 = (t61 * t105 + t90) * qJD(3) + t57 * t101;
t19 = qJD(3) * pkin(8) + t21;
t30 = t55 * t45 - t52 * t95;
t10 = t61 * t19 + t58 * t30;
t16 = t134 * qJD(3);
t4 = t10 * qJD(4) + t58 * t16;
t130 = t4 * t57;
t129 = t4 * t60;
t104 = t60 * qJD(4);
t88 = t61 * t104;
t113 = qJD(3) * t58;
t89 = t57 * t113;
t28 = qJD(3) * t88 - qJD(5) * t89 + t60 * t101;
t128 = t28 * t57;
t36 = t89 - t104;
t127 = t36 * t46;
t38 = t60 * t113 + t105;
t126 = t38 * t46;
t125 = t46 * t57;
t64 = qJD(3) ^ 2;
t122 = t52 * t64;
t120 = t57 * t61;
t119 = t60 * t46;
t118 = t60 * t61;
t81 = pkin(4) * t58 - pkin(9) * t61;
t41 = t81 * qJD(4);
t117 = t21 - t41;
t49 = t58 ^ 2;
t116 = -t61 ^ 2 + t49;
t115 = qJD(3) * pkin(3);
t110 = qJD(4) * t58;
t109 = qJD(4) * t61;
t108 = qJD(5) * t57;
t102 = qJD(3) * qJD(4);
t98 = t58 * t124;
t97 = t61 * t124;
t93 = t52 * t111;
t92 = t58 * t108;
t91 = t46 * t107;
t86 = t58 * t102;
t18 = -t20 - t115;
t85 = -qJD(3) * t18 - t16;
t83 = t61 * t93;
t82 = t58 * t93;
t43 = -t61 * pkin(4) - t58 * pkin(9) - pkin(3);
t15 = t43 * qJD(3) - t20;
t8 = qJD(4) * pkin(9) + t10;
t1 = t60 * t15 - t57 * t8;
t2 = t57 * t15 + t60 * t8;
t25 = t56 * t124 + t68;
t31 = -t53 * t54 * t52 + t56 * t55;
t14 = t25 * t61 + t31 * t58;
t80 = -t135 * t57 + t14 * t60;
t79 = -t135 * t60 - t14 * t57;
t78 = t58 * t19 - t61 * t30;
t13 = t25 * t58 - t31 * t61;
t77 = qJD(3) * t49 - t46 * t61;
t33 = t58 * t55 + t97;
t76 = -t60 * t123 - t57 * t33;
t75 = t57 * t123 - t60 * t33;
t32 = -t61 * t55 + t98;
t63 = qJD(4) ^ 2;
t73 = pkin(8) * t63 - t132;
t72 = qJD(4) * (t18 + t20 - t115);
t3 = -t78 * qJD(4) + t61 * t16;
t7 = -qJD(4) * pkin(4) + t78;
t66 = qJD(4) * t7 + qJD(5) * t15 - t20 * t46 + t3;
t40 = t81 * qJD(3);
t27 = t33 * qJD(4) + t82;
t26 = -t32 * qJD(4) + t83;
t23 = t25 * qJD(3);
t22 = t135 * qJD(3);
t12 = qJD(3) * t41 + t17;
t11 = t60 * t12;
t6 = -t13 * qJD(4) + t22 * t61;
t5 = t14 * qJD(4) + t22 * t58;
t9 = [0, 0, 0, -t23 * qJD(3), -t22 * qJD(3), 0, 0, 0, 0, 0, -t5 * qJD(4) + (-t110 * t135 - t23 * t61) * qJD(3), -t6 * qJD(4) + (-t109 * t135 + t23 * t58) * qJD(3), 0, 0, 0, 0, 0, -(-t80 * qJD(5) + t23 * t60 - t6 * t57) * t46 + t79 * t86 + t5 * t36 + t13 * t29, (t79 * qJD(5) + t23 * t57 + t6 * t60) * t46 - t80 * t86 + t5 * t38 + t13 * t28; 0, 0, 0, -t59 * t122, -t62 * t122, 0, 0, 0, 0, 0, -t64 * t97 + (-t27 - t82) * qJD(4), t64 * t98 + (-t26 - t83) * qJD(4), 0, 0, 0, 0, 0, -(t75 * qJD(5) - t57 * t26 + t60 * t94) * t46 + t76 * t86 + t27 * t36 + t32 * t29, (t76 * qJD(5) + t60 * t26 + t57 * t94) * t46 + t75 * t86 + t27 * t38 + t32 * t28; 0, 0, 0, t132, (-t134 + t20) * qJD(3), 0.2e1 * t61 * t86, -0.2e1 * t116 * t102, t63 * t61, -t63 * t58, 0, t58 * t72 - t73 * t61, t73 * t58 + t61 * t72, t28 * t60 * t58 + (t88 - t92) * t38, (-t36 * t60 - t38 * t57) * t109 + (-t128 - t29 * t60 + (t36 * t57 - t38 * t60) * qJD(5)) * t58, t46 * t92 - t28 * t61 + (t58 * t38 + t77 * t60) * qJD(4), t46 * t90 + t29 * t61 + (-t58 * t36 - t77 * t57) * qJD(4), (-t46 - t103) * t110, (t43 * t108 + t117 * t60) * t46 + (t8 * t107 - t11 + (qJD(4) * t36 + t91) * pkin(8) + t66 * t57) * t61 + (t7 * t107 + pkin(8) * t29 - t20 * t36 + t130 + (-pkin(8) * t125 + (-pkin(8) * t120 + t60 * t43) * qJD(3) + t1) * qJD(4)) * t58, (t43 * t107 - t117 * t57) * t46 + (qJD(4) * pkin(8) * t38 + (t12 + (-pkin(8) * t46 - t8) * qJD(5)) * t57 + t66 * t60) * t61 + (-t7 * t108 + pkin(8) * t28 - t20 * t38 + t129 + (-pkin(8) * t119 - (pkin(8) * t118 + t57 * t43) * qJD(3) - t2) * qJD(4)) * t58; 0, 0, 0, 0, 0, -t58 * t64 * t61, t116 * t64, 0, 0, 0, t85 * t58, t85 * t61, -t38 * t119 + t128, (t28 + t127) * t60 + (-t29 + t126) * t57, -t91 + (t46 * t118 + (-t38 + t105) * t58) * qJD(3), t46 * t108 + (-t46 * t120 + (t36 + t104) * t58) * qJD(3), t46 * t113, -pkin(4) * t29 - t129 + (t60 * t40 + t57 * t78) * t46 - t10 * t36 + (pkin(9) * t119 + t7 * t57) * qJD(5) + (-t1 * t58 + (-pkin(9) * t110 - t61 * t7) * t57) * qJD(3), -pkin(4) * t28 + t130 - (t57 * t40 - t60 * t78) * t46 - t10 * t38 + (-pkin(9) * t125 + t7 * t60) * qJD(5) + (-t7 * t118 + (-pkin(9) * t104 + t2) * t58) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t36, -t36 ^ 2 + t38 ^ 2, t28 - t127, -t126 - t29, t86, t131 * t2 - t57 * t3 - t7 * t38 + t11, t131 * t1 - t57 * t12 - t60 * t3 + t7 * t36;];
tauc_reg = t9;
