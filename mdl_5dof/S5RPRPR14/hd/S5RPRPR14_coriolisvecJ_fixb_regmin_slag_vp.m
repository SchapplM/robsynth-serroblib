% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:03
% EndTime: 2021-01-15 12:17:10
% DurationCPUTime: 1.14s
% Computational Cost: add. (1265->197), mult. (2809->286), div. (0->0), fcn. (1868->6), ass. (0->115)
t122 = cos(pkin(8));
t76 = sin(pkin(8));
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t91 = -t122 * t78 - t76 * t80;
t146 = t91 * qJD(1);
t148 = qJD(5) - t146;
t79 = cos(qJ(5));
t114 = t79 * qJD(3);
t121 = qJD(1) * t78;
t105 = t76 * t121;
t101 = t122 * t80;
t94 = qJD(1) * t101;
t46 = t94 - t105;
t77 = sin(qJ(5));
t26 = t77 * t46 - t114;
t151 = t148 * t26;
t81 = -pkin(1) - pkin(6);
t56 = t81 * qJD(1) + qJD(2);
t150 = -qJ(4) * qJD(1) + t56;
t55 = qJD(3) * t94;
t36 = qJD(3) * t105 - t55;
t130 = t77 * t36;
t99 = t79 * t148;
t149 = t148 * t99 - t130;
t37 = t146 * qJD(3);
t145 = -qJD(5) + t148;
t141 = pkin(3) * t80;
t106 = qJD(1) * t141;
t109 = qJ(4) * qJD(3);
t115 = t78 * qJD(4);
t119 = qJD(3) * t80;
t25 = t56 * t119 + (-t80 * t109 - t115) * qJD(1);
t120 = qJD(3) * t78;
t113 = t80 * qJD(4);
t90 = t78 * t109 - t113;
t84 = t90 * qJD(1) - t56 * t120;
t3 = -t122 * t84 + t76 * t25;
t64 = t76 * pkin(3) + pkin(7);
t144 = (t46 * pkin(4) - pkin(7) * t146 + qJD(5) * t64 + t106) * t148 + t3;
t53 = pkin(3) * t121 + qJD(1) * qJ(2) + qJD(4);
t13 = -pkin(4) * t146 - t46 * pkin(7) + t53;
t123 = qJ(4) - t81;
t103 = t123 * t80;
t54 = t123 * t78;
t24 = -t76 * t103 - t122 * t54;
t136 = t24 * t36;
t131 = t76 * t78;
t49 = t101 - t131;
t140 = t3 * t49;
t38 = -qJD(3) * t103 - t115;
t88 = t123 * t120 - t113;
t15 = t122 * t38 + t76 * t88;
t65 = t78 * pkin(3) + qJ(2);
t21 = -pkin(4) * t91 - t49 * pkin(7) + t65;
t102 = t122 * t25;
t4 = t76 * t84 + t102;
t96 = qJD(3) * t122;
t45 = -t76 * t119 - t78 * t96;
t39 = t150 * t78;
t132 = t76 * t39;
t40 = t150 * t80;
t35 = qJD(3) * pkin(3) + t40;
t11 = t122 * t35 - t132;
t7 = -qJD(3) * pkin(4) - t11;
t143 = (qJD(5) * t13 + t4) * t91 - (qJD(5) * t21 + t15) * t148 + t7 * t45 + t136 + t140;
t72 = qJD(1) * qJD(2);
t142 = 0.2e1 * t72;
t117 = qJD(5) * t77;
t9 = qJD(5) * t114 - t46 * t117 + t79 * t37;
t139 = t49 * t9;
t138 = t9 * t77;
t137 = t21 * t36;
t28 = t77 * qJD(3) + t79 * t46;
t135 = t28 * t46;
t134 = t46 * t26;
t133 = t91 * t36;
t129 = t77 * t37;
t31 = t79 * t36;
t82 = qJD(3) ^ 2;
t128 = t82 * t78;
t127 = t82 * t80;
t33 = t122 * t39;
t12 = t76 * t35 + t33;
t108 = qJD(1) * qJD(3);
t104 = t80 * t108;
t52 = pkin(3) * t104 + t72;
t126 = t78 ^ 2 - t80 ^ 2;
t83 = qJD(1) ^ 2;
t125 = -t82 - t83;
t124 = t83 * qJ(2);
t118 = qJD(5) * t49;
t116 = t53 * qJD(1);
t57 = pkin(3) * t119 + qJD(2);
t111 = qJ(2) * qJD(3);
t107 = 0.2e1 * qJD(1);
t95 = -qJD(5) * t91 + qJD(1);
t8 = qJD(3) * pkin(7) + t12;
t1 = t79 * t13 - t77 * t8;
t2 = t77 * t13 + t79 * t8;
t93 = -t31 + (t146 * t77 - t117) * t148;
t92 = -t49 * t117 + t79 * t45;
t17 = t122 * t40 - t132;
t86 = t64 * t36 + (t17 + t7) * t148;
t44 = t76 * t120 - t80 * t96;
t85 = t11 * t45 - t12 * t44 - t4 * t91 - t140;
t66 = -t122 * pkin(3) - pkin(4);
t23 = t123 * t101 - t76 * t54;
t18 = -t44 * pkin(4) - t45 * pkin(7) + t57;
t16 = t76 * t40 + t33;
t14 = -t122 * t88 + t76 * t38;
t10 = t28 * qJD(5) + t129;
t6 = -t36 * pkin(4) - t37 * pkin(7) + t52;
t5 = t79 * t6;
t19 = [0, 0, 0, 0, t142, qJ(2) * t142, -0.2e1 * t78 * t104, 0.2e1 * t126 * t108, -t128, -t127, 0, -t81 * t128 + (qJD(2) * t78 + t80 * t111) * t107, -t81 * t127 + (qJD(2) * t80 - t78 * t111) * t107, -t14 * qJD(3) - t146 * t57 - t65 * t36 - t53 * t44 - t52 * t91, -t15 * qJD(3) + t65 * t37 + t53 * t45 + t57 * t46 + t52 * t49, t14 * t46 + t146 * t15 + t23 * t37 + t136 - t85, -t11 * t14 + t12 * t15 + t3 * t23 + t4 * t24 + t52 * t65 + t53 * t57, t79 * t139 + t92 * t28, (-t26 * t79 - t28 * t77) * t45 + (-t10 * t79 - t138 + (t26 * t77 - t28 * t79) * qJD(5)) * t49, t148 * t92 - t28 * t44 - t49 * t31 - t9 * t91, t49 * t130 + t10 * t91 + t26 * t44 + (-t79 * t118 - t77 * t45) * t148, -t148 * t44 + t133, -t1 * t44 + t23 * t10 + t14 * t26 - t5 * t91 + (t18 * t148 - t137 + (-t148 * t24 + t7 * t49 + t8 * t91) * qJD(5)) * t79 + t143 * t77, t14 * t28 + t2 * t44 + t23 * t9 + (-(-qJD(5) * t24 + t18) * t148 + t137 + (-qJD(5) * t8 + t6) * t91 - t7 * t118) * t77 + t143 * t79; 0, 0, 0, 0, -t83, -t124, 0, 0, 0, 0, 0, t125 * t78, t125 * t80, qJD(1) * t146 + t45 * qJD(3), -qJD(1) * t46 + t44 * qJD(3), -t146 * t44 - t49 * t37 - t45 * t46 - t133, t85 - t116, 0, 0, 0, 0, 0, -t91 * t130 - t49 * t10 - t45 * t26 + (t44 * t77 - t95 * t79) * t148, -t91 * t31 - t45 * t28 - t139 + (t44 * t79 + t95 * t77) * t148; 0, 0, 0, 0, 0, 0, t80 * t83 * t78, -t126 * t83, 0, 0, 0, -t80 * t124, t78 * t124, t16 * qJD(3) + t106 * t146 - t53 * t46 - t3, -t102 - t53 * t146 + (t56 * t131 + t17) * qJD(3) + (-t46 * t141 - t76 * t90) * qJD(1), (t12 - t16) * t46 - (-t11 + t17) * t146 + (-t122 * t37 + t36 * t76) * pkin(3), t11 * t16 - t12 * t17 + (-t80 * t116 - t122 * t3 + t4 * t76) * pkin(3), t28 * t99 + t138, (t9 - t151) * t79 + (-t148 * t28 - t10) * t77, -t135 + t149, t93 + t134, -t148 * t46, -t1 * t46 + t66 * t10 - t144 * t79 - t16 * t26 + t86 * t77, t144 * t77 - t16 * t28 + t2 * t46 + t66 * t9 + t86 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 + (t46 - t105) * qJD(3), 0.2e1 * t37, -t146 ^ 2 - t46 ^ 2, t11 * t46 - t12 * t146 + t52, 0, 0, 0, 0, 0, t93 - t134, -t135 - t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t9 + t151, t145 * t28 - t129, -t36, t145 * t2 - t7 * t28 - t77 * t4 + t5, t1 * t145 + t7 * t26 - t79 * t4 - t77 * t6;];
tauc_reg = t19;
