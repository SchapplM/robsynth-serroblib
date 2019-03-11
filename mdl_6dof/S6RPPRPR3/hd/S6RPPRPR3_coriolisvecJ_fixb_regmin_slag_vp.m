% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:58
% EndTime: 2019-03-09 01:45:01
% DurationCPUTime: 1.12s
% Computational Cost: add. (1586->186), mult. (3338->268), div. (0->0), fcn. (2233->8), ass. (0->107)
t60 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t126 = qJ(5) - t60;
t153 = -qJD(1) * t126 + qJD(3);
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t32 = -t81 * qJD(2) + t153 * t83;
t125 = cos(pkin(10));
t77 = sin(pkin(10));
t92 = t125 * t81 + t77 * t83;
t47 = t92 * qJD(1);
t152 = qJD(6) + t47;
t82 = cos(qJ(6));
t117 = t82 * qJD(4);
t103 = t125 * t83;
t124 = qJD(1) * t81;
t49 = qJD(1) * t103 - t77 * t124;
t80 = sin(qJ(6));
t35 = t80 * t49 - t117;
t154 = t152 * t35;
t144 = 0.2e1 * qJD(3);
t62 = sin(pkin(9)) * pkin(1) + qJ(3);
t56 = qJD(1) * t62;
t101 = t152 * t82;
t42 = t49 * qJD(4);
t134 = t80 * t42;
t151 = -t101 * t152 - t134;
t110 = qJD(1) * qJD(5);
t146 = -t83 * qJD(2) - t153 * t81;
t150 = t146 * qJD(4) - t83 * t110;
t149 = -qJD(6) + t152;
t26 = qJD(4) * t32 - t81 * t110;
t3 = -t125 * t150 + t77 * t26;
t66 = t77 * pkin(4) + pkin(8);
t148 = t152 * (t83 * qJD(1) * pkin(4) + t49 * pkin(5) + t47 * pkin(8) + qJD(6) * t66) + t3;
t39 = t82 * t42;
t55 = -t77 * t81 + t103;
t122 = qJD(6) * t80;
t51 = t92 * qJD(4);
t93 = t55 * t122 + t82 * t51;
t147 = -t152 * t93 + t55 * t39;
t104 = t126 * t83;
t38 = -qJD(4) * t104 - t81 * qJD(5);
t118 = t81 * qJD(4);
t90 = -t83 * qJD(5) + t126 * t118;
t14 = t125 * t38 + t77 * t90;
t52 = t126 * t81;
t25 = -t77 * t104 - t125 * t52;
t140 = t25 * t42;
t143 = t3 * t55;
t46 = pkin(4) * t124 + qJD(5) + t56;
t15 = t47 * pkin(5) - t49 * pkin(8) + t46;
t108 = t81 * pkin(4) + t62;
t23 = pkin(5) * t92 - t55 * pkin(8) + t108;
t4 = t125 * t26 + t150 * t77;
t135 = t77 * t146;
t31 = qJD(4) * pkin(4) + t32;
t7 = t125 * t31 + t135;
t5 = -qJD(4) * pkin(5) - t7;
t145 = -(qJD(6) * t23 + t14) * t152 - (qJD(6) * t15 + t4) * t92 - t5 * t51 - t140 + t143;
t43 = qJD(1) * t51;
t17 = qJD(6) * t117 - t49 * t122 - t82 * t43;
t142 = t17 * t80;
t141 = t23 * t42;
t37 = t80 * qJD(4) + t82 * t49;
t139 = t37 * t49;
t138 = t49 * t35;
t137 = t92 * t42;
t136 = t55 * t17;
t133 = t80 * t43;
t84 = qJD(4) ^ 2;
t132 = t84 * t81;
t131 = t84 * t83;
t48 = -qJD(4) * t103 + t77 * t118;
t130 = t17 * t92 - t37 * t48;
t29 = t125 * t146;
t8 = t77 * t31 - t29;
t111 = qJD(1) * qJD(4);
t106 = t83 * t111;
t74 = qJD(3) * qJD(1);
t129 = pkin(4) * t106 + t74;
t128 = t81 ^ 2 - t83 ^ 2;
t85 = qJD(1) ^ 2;
t127 = -t84 - t85;
t123 = qJD(6) * t55;
t121 = t46 * qJD(1);
t120 = t56 * qJD(1);
t115 = t83 * qJD(4);
t114 = pkin(4) * t115 + qJD(3);
t97 = -qJD(6) * t92 - qJD(1);
t6 = qJD(4) * pkin(8) + t8;
t1 = t82 * t15 - t80 * t6;
t2 = t80 * t15 + t82 * t6;
t18 = t37 * qJD(6) - t133;
t95 = t18 * t92 - t48 * t35;
t94 = t39 + (-t47 * t80 - t122) * t152;
t10 = t125 * t32 + t135;
t89 = -t66 * t42 + (t10 + t5) * t152;
t88 = t4 * t92 - t8 * t48 - t7 * t51 - t143;
t86 = (-t82 * t123 + t80 * t51) * t152 - t55 * t134;
t67 = -t125 * pkin(4) - pkin(5);
t24 = t126 * t103 - t77 * t52;
t19 = -t48 * pkin(5) + t51 * pkin(8) + t114;
t16 = t42 * pkin(5) + t43 * pkin(8) + t129;
t13 = -t125 * t90 + t77 * t38;
t12 = t82 * t16;
t9 = t77 * t32 - t29;
t11 = [0, 0, 0, 0, 0, 0.2e1 * t74, t56 * t144, -0.2e1 * t81 * t106, 0.2e1 * t128 * t111, -t132, -t131, 0, t56 * t115 - t60 * t132 + (t62 * t115 + t81 * t144) * qJD(1), -t56 * t118 - t60 * t131 + (-t62 * t118 + t83 * t144) * qJD(1), t13 * t49 - t14 * t47 - t24 * t43 - t140 - t88, t129 * t108 + t46 * t114 - t7 * t13 + t8 * t14 + t3 * t24 + t4 * t25, t82 * t136 - t93 * t37 (t35 * t82 + t37 * t80) * t51 + (-t142 - t18 * t82 + (t35 * t80 - t37 * t82) * qJD(6)) * t55, t130 + t147, t86 - t95, -t152 * t48 + t137, -t1 * t48 + t12 * t92 + t13 * t35 + t24 * t18 + (t19 * t152 + t141 + (-t152 * t25 + t5 * t55 - t6 * t92) * qJD(6)) * t82 + t145 * t80, t13 * t37 + t24 * t17 + t2 * t48 + (-(-qJD(6) * t25 + t19) * t152 - t141 - (-qJD(6) * t6 + t16) * t92 - t5 * t123) * t80 + t145 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t132, -t55 * t42 - t43 * t92 + t51 * t47 - t48 * t49, t3 * t92 + t4 * t55 + t7 * t48 - t8 * t51, 0, 0, 0, 0, 0, t86 + t95, t130 - t147; 0, 0, 0, 0, 0, -t85, -t120, 0, 0, 0, 0, 0, t127 * t81, t127 * t83, t55 * t43 + t48 * t47 + t51 * t49 - t137, t88 - t121, 0, 0, 0, 0, 0, -t92 * t134 - t55 * t18 + t51 * t35 + (t48 * t80 + t82 * t97) * t152, -t92 * t39 - t136 + t51 * t37 + (t48 * t82 - t80 * t97) * t152; 0, 0, 0, 0, 0, 0, 0, t83 * t85 * t81, -t128 * t85, 0, 0, 0, -t83 * t120, t81 * t120 (t8 - t9) * t49 - (-t10 + t7) * t47 + (t125 * t43 - t42 * t77) * pkin(4), -t8 * t10 + t7 * t9 + (-t83 * t121 - t125 * t3 + t4 * t77) * pkin(4), t37 * t101 + t142 (t17 - t154) * t82 + (-t152 * t37 - t18) * t80, -t139 - t151, t94 + t138, -t152 * t49, -t1 * t49 - t148 * t82 + t67 * t18 - t9 * t35 + t89 * t80, t148 * t80 + t67 * t17 + t2 * t49 - t9 * t37 + t89 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 ^ 2 - t49 ^ 2, t8 * t47 + t7 * t49 + t129, 0, 0, 0, 0, 0, t94 - t138, -t139 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t17 + t154, t149 * t37 + t133, t42, t149 * t2 - t5 * t37 - t80 * t4 + t12, t149 * t1 - t80 * t16 + t5 * t35 - t82 * t4;];
tauc_reg  = t11;
