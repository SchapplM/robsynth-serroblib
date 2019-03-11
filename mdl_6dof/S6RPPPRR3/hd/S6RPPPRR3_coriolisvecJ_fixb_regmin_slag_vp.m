% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:59
% EndTime: 2019-03-09 01:34:02
% DurationCPUTime: 1.04s
% Computational Cost: add. (1458->180), mult. (3151->272), div. (0->0), fcn. (2244->8), ass. (0->114)
t78 = sin(pkin(10));
t84 = sin(qJ(5));
t130 = t84 * t78;
t111 = qJD(1) * t130;
t80 = cos(pkin(10));
t123 = qJD(1) * t80;
t86 = cos(qJ(5));
t60 = t86 * t123;
t49 = -t60 + t111;
t115 = qJD(6) - t49;
t85 = cos(qJ(6));
t116 = t85 * qJD(5);
t54 = t86 * t78 + t84 * t80;
t50 = qJD(1) * t54;
t83 = sin(qJ(6));
t28 = -t83 * t50 - t116;
t149 = t115 * t28;
t104 = t115 * t85;
t52 = t54 * qJD(5);
t44 = qJD(1) * t52;
t131 = t83 * t44;
t148 = t104 * t115 - t131;
t53 = -t86 * t80 + t130;
t81 = cos(pkin(9));
t62 = t81 * qJD(2) - qJD(4);
t57 = t62 * qJD(1);
t147 = t53 * t57;
t146 = t54 * t57;
t125 = t78 ^ 2 + t80 ^ 2;
t109 = t125 * t57;
t144 = -qJD(6) + t115;
t114 = qJD(1) * qJ(2);
t87 = -pkin(1) - pkin(2);
t59 = t87 * qJD(1) + qJD(2);
t79 = sin(pkin(9));
t46 = t81 * t114 + t79 * t59;
t34 = -qJD(1) * qJ(4) + t46;
t71 = t80 * qJD(3);
t20 = t71 + (pkin(7) * qJD(1) - t34) * t78;
t25 = t78 * qJD(3) + t80 * t34;
t21 = -pkin(7) * t123 + t25;
t8 = t84 * t20 + t86 * t21;
t4 = t8 * qJD(5) + t146;
t143 = (-t50 * pkin(5) + t115 * pkin(8)) * t115 + t4;
t122 = qJD(6) * t54;
t51 = t53 * qJD(5);
t142 = -(t85 * t122 - t83 * t51) * t115 + t54 * t131;
t45 = -t79 * t114 + t81 * t59;
t31 = qJD(1) * pkin(3) + qJD(4) - t45;
t26 = pkin(4) * t123 + t31;
t11 = -t49 * pkin(5) + t50 * pkin(8) + t26;
t126 = t81 * qJ(2) + t79 * t87;
t55 = -qJ(4) + t126;
t139 = pkin(7) - t55;
t39 = t139 * t78;
t40 = t139 * t80;
t14 = t84 * t39 - t86 * t40;
t108 = -t79 * qJ(2) + t81 * t87;
t102 = pkin(3) - t108;
t48 = t80 * pkin(4) + t102;
t18 = -t53 * pkin(5) + t54 * pkin(8) + t48;
t7 = t86 * t20 - t84 * t21;
t3 = t7 * qJD(5) - t147;
t5 = -qJD(5) * pkin(5) - t7;
t97 = t86 * t39 + t84 * t40;
t9 = t97 * qJD(5) - t53 * t62;
t141 = -(qJD(6) * t18 + t9) * t115 + (qJD(6) * t11 + t3) * t53 + t14 * t44 - t4 * t54 + t5 * t51;
t113 = qJD(1) * qJD(2);
t63 = t79 * t113;
t140 = 0.2e1 * t63;
t121 = qJD(6) * t83;
t58 = qJD(5) * t111;
t43 = -qJD(5) * t60 + t58;
t15 = qJD(6) * t116 + t50 * t121 + t85 * t43;
t138 = t15 * t83;
t137 = t18 * t44;
t30 = t83 * qJD(5) - t85 * t50;
t136 = t30 * t50;
t135 = t50 * t28;
t88 = qJD(1) ^ 2;
t134 = t79 * t88;
t133 = t81 * t88;
t132 = t83 * t43;
t33 = t85 * t44;
t124 = qJD(1) * t53;
t128 = -t81 * t124 + t79 * t52;
t42 = t53 * t79;
t127 = -qJD(5) * t42 - t81 * t50;
t120 = t51 * qJD(5);
t119 = t52 * qJD(5);
t118 = t79 * qJD(1);
t117 = t79 * qJD(2);
t110 = 0.2e1 * t113;
t103 = -0.2e1 * t50;
t101 = -qJD(6) * t81 - t128;
t6 = qJD(5) * pkin(8) + t8;
t1 = t85 * t11 - t83 * t6;
t2 = t83 * t11 + t85 * t6;
t100 = -t53 * t15 - t52 * t30;
t16 = t30 * qJD(6) + t132;
t99 = t53 * t16 + t52 * t28;
t98 = -(-t78 * t34 + t71) * t78 + t25 * t80;
t96 = t45 * t79 - t46 * t81;
t95 = -qJD(6) * t42 + t118;
t94 = -t33 + (t49 * t83 - t121) * t115;
t92 = t54 * t121 + t85 * t51;
t90 = pkin(8) * t44 + (t5 + t7) * t115;
t89 = t115 * t92 + t54 * t33;
t41 = t54 * t79;
t19 = -t52 * pkin(5) - t51 * pkin(8) + t117;
t17 = -t44 * pkin(5) - t43 * pkin(8) + t63;
t12 = t85 * t17;
t10 = t14 * qJD(5) + t54 * t62;
t13 = [0, 0, 0, 0, t110, qJ(2) * t110, t140, t81 * t110 ((-t79 * t108 + t81 * t126) * qJD(1) - t96) * qJD(2), t80 * t140, -0.2e1 * t78 * t63, -0.2e1 * t109, t98 * t62 + t55 * t109 + (qJD(1) * t102 + t31) * t117, -t43 * t54 - t50 * t51, t43 * t53 - t54 * t44 + t51 * t49 - t50 * t52, t120, t119, 0, -t10 * qJD(5) - t26 * t52 - t48 * t44 + (-t49 - t124) * t117, -t9 * qJD(5) + t103 * t117 + t26 * t51 + t48 * t43, -t15 * t85 * t54 + t92 * t30 (-t28 * t85 - t30 * t83) * t51 + (t138 + t16 * t85 + (-t28 * t83 + t30 * t85) * qJD(6)) * t54, t100 + t89, -t142 + t99, -t115 * t52 + t44 * t53, -t1 * t52 + t10 * t28 - t12 * t53 - t97 * t16 + (-t137 + t19 * t115 + (-t115 * t14 - t5 * t54 + t6 * t53) * qJD(6)) * t85 + t141 * t83, t10 * t30 - t97 * t15 + t2 * t52 + (-(-qJD(6) * t14 + t19) * t115 + t137 + (-qJD(6) * t6 + t17) * t53 + t5 * t122) * t83 + t141 * t85; 0, 0, 0, 0, -t88, -t88 * qJ(2), -t134, -t133, t96 * qJD(1), -t80 * t134, t78 * t134, t125 * t133, t79 * t109 + (-t31 * t79 + (-t98 - t117) * t81) * qJD(1), 0, 0, 0, 0, 0, -t127 * qJD(5) + t49 * t118 + t81 * t44, t128 * qJD(5) + t50 * t118 - t81 * t43, 0, 0, 0, 0, 0 -(t83 * t42 - t85 * t81) * t44 + t41 * t16 - (t101 * t83 + t85 * t95) * t115 + t127 * t28 (-t85 * t42 - t83 * t81) * t44 + t41 * t15 - (t101 * t85 - t83 * t95) * t115 + t127 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t120, 0, 0, 0, 0, 0, t142 + t99, -t100 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t88, t98 * qJD(1) + t63, 0, 0, 0, 0, 0, t103 * qJD(5), t58 + (t49 - t60) * qJD(5), 0, 0, 0, 0, 0, t94 + t135, t136 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t49, -t49 ^ 2 + t50 ^ 2, t58 + (-t49 - t60) * qJD(5), 0, 0, t26 * t50 - t146, -t26 * t49 + t147, t30 * t104 + t138 (t15 - t149) * t85 + (-t115 * t30 - t16) * t83, t136 + t148, t94 - t135, t115 * t50, -pkin(5) * t16 + t1 * t50 - t143 * t85 - t8 * t28 + t90 * t83, -pkin(5) * t15 + t143 * t83 - t2 * t50 - t8 * t30 + t90 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t28, -t28 ^ 2 + t30 ^ 2, t15 + t149, t144 * t30 - t132, -t44, t144 * t2 - t83 * t3 - t5 * t30 + t12, t144 * t1 - t83 * t17 + t5 * t28 - t85 * t3;];
tauc_reg  = t13;
