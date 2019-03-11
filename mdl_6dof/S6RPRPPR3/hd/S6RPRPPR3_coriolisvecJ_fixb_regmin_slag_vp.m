% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:28
% DurationCPUTime: 1.34s
% Computational Cost: add. (1040->234), mult. (2288->308), div. (0->0), fcn. (1219->6), ass. (0->131)
t58 = -cos(pkin(9)) * pkin(1) - pkin(2);
t48 = qJD(1) * t58;
t57 = sin(pkin(9)) * pkin(1) + pkin(7);
t47 = t57 * qJD(1);
t91 = sin(qJ(3));
t42 = t91 * t47;
t93 = cos(qJ(3));
t153 = -t93 * qJD(2) + t42;
t171 = qJD(4) + t153;
t23 = -qJD(3) * pkin(3) + t171;
t94 = -pkin(3) - pkin(4);
t121 = qJD(3) * t94;
t133 = qJ(5) * qJD(1);
t21 = -t91 * t133 + t153;
t136 = -qJD(4) - t21;
t10 = t121 - t136;
t92 = cos(qJ(6));
t138 = t92 * qJD(3);
t146 = qJD(1) * t93;
t90 = sin(qJ(6));
t44 = t90 * t146 - t138;
t33 = t91 * qJD(2) + t93 * t47;
t22 = -t93 * t133 + t33;
t83 = qJD(3) * qJ(4);
t15 = -t22 - t83;
t27 = t83 + t33;
t147 = -qJ(5) + t57;
t140 = t91 * qJD(1);
t56 = qJD(6) + t140;
t170 = -qJD(6) + t56;
t169 = (-t153 + t42) * qJD(3);
t11 = qJD(3) * pkin(5) - t15;
t107 = t91 * pkin(5) + t93 * pkin(8);
t37 = -t93 * pkin(3) - t91 * qJ(4) + t58;
t31 = t93 * pkin(4) - t37;
t17 = t107 + t31;
t139 = t91 * qJD(5);
t144 = qJD(3) * t93;
t26 = t147 * t144 - t139;
t134 = qJ(4) * qJD(1);
t28 = -pkin(3) * t146 - t91 * t134 + t48;
t13 = pkin(4) * t146 + qJD(5) - t28;
t5 = t107 * qJD(1) + t13;
t132 = qJ(5) * qJD(3);
t122 = t93 * t132;
t130 = qJD(2) * qJD(3);
t30 = t91 * t130 + t47 * t144;
t8 = (-t122 - t139) * qJD(1) + t30;
t168 = -(qJD(6) * t17 + t26) * t56 + (t11 * qJD(3) - qJD(6) * t5 - t8) * t91;
t143 = qJD(6) * t90;
t125 = t93 * t143;
t158 = t56 * t91;
t85 = t93 ^ 2;
t167 = -(qJD(1) * t85 - t158) * t138 + t56 * t125;
t137 = t93 * qJD(5);
t145 = qJD(3) * t91;
t61 = t93 * t130;
t81 = qJD(3) * qJD(4);
t18 = -t47 * t145 + t61 + t81;
t131 = qJD(1) * qJD(3);
t120 = t91 * t131;
t54 = qJ(5) * t120;
t6 = qJD(1) * t137 - t18 - t54;
t166 = t6 * t90;
t165 = t6 * t92;
t164 = t15 * t93;
t19 = qJD(6) * t44 + t92 * t120;
t163 = t19 * t90;
t141 = t90 * qJD(3);
t45 = t92 * t146 + t141;
t20 = t45 * qJD(6) - t90 * t120;
t162 = t20 * t91;
t161 = t44 * t56;
t160 = t44 * t93;
t159 = t45 * t56;
t157 = t56 * t92;
t95 = qJD(3) ^ 2;
t156 = t95 * t91;
t74 = t95 * t93;
t155 = -t45 * t144 + t19 * t91;
t124 = qJD(6) * t157;
t154 = t90 * t85 * t131 + t93 * t124;
t119 = t93 * t131;
t70 = t91 * qJD(4);
t152 = qJ(4) * t119 + qJD(1) * t70;
t151 = t61 + 0.2e1 * t81;
t150 = t93 * t83 + t70;
t84 = t91 ^ 2;
t149 = t84 - t85;
t82 = -pkin(8) + t94;
t7 = t82 * qJD(3) - t136;
t148 = qJD(6) * t7;
t142 = t11 * qJD(6);
t135 = -qJD(5) - t13;
t129 = t90 * t158;
t128 = t91 * t157;
t96 = qJD(1) ^ 2;
t127 = t91 * t96 * t93;
t126 = t56 * t143;
t117 = t135 * t91;
t116 = qJD(1) * t31 + t13;
t108 = t91 * t121;
t14 = qJD(1) * t108 + t152;
t24 = t108 + t150;
t115 = qJD(1) * t24 + t14;
t114 = qJD(1) * t37 + t28;
t110 = 0.2e1 * t119;
t1 = t92 * t5 - t90 * t7;
t2 = t90 * t5 + t92 * t7;
t106 = t33 * qJD(3) - t30;
t38 = t147 * t91;
t100 = pkin(5) * t93 + t82 * t91;
t98 = t100 * qJD(3);
t105 = (-qJD(6) * t38 + t150 + t98) * t56;
t103 = 0.2e1 * qJD(3) * t48;
t102 = -t11 * t91 - t82 * t144;
t29 = pkin(3) * t120 - t152;
t36 = pkin(3) * t145 - t150;
t101 = -qJD(1) * t36 - t57 * t95 - t29;
t97 = t18 * t93 + t30 * t91 + (t23 * t93 - t27 * t91) * qJD(3);
t89 = qJ(4) + pkin(5);
t65 = t93 * t134;
t53 = -t84 * t96 - t95;
t46 = pkin(3) * t140 - t65;
t39 = t147 * t93;
t34 = t94 * t140 + t65;
t25 = -t91 * t132 + t57 * t145 + t137;
t12 = t100 * qJD(1) + t65;
t4 = qJD(1) * t98 + t152;
t3 = t92 * t4;
t9 = [0, 0, 0, 0, t91 * t110, -0.2e1 * t149 * t131, t74, -t156, 0, t91 * t103 - t57 * t74, t93 * t103 + t57 * t156, t101 * t93 + t114 * t145, t97, t101 * t91 - t114 * t144, t28 * t36 + t29 * t37 + t97 * t57, t115 * t91 + (t116 * t93 - t25) * qJD(3), -t115 * t93 + (t116 * t91 + t26) * qJD(3), t6 * t93 - t8 * t91 + (-t10 * t93 - t15 * t91) * qJD(3) + (t25 * t93 - t26 * t91 + (-t38 * t93 + t39 * t91) * qJD(3)) * qJD(1), t10 * t26 + t13 * t24 + t14 * t31 + t15 * t25 + t8 * t38 - t6 * t39, -t19 * t92 * t93 + (-t91 * t138 - t125) * t45 (t44 * t92 + t45 * t90) * t145 + (t163 - t20 * t92 + (t44 * t90 - t45 * t92) * qJD(6)) * t93, t155 + t167, t162 + (-t129 + t160) * qJD(3) + t154 (t56 + t140) * t144, -t39 * t20 + t25 * t44 + t3 * t91 + (-t91 * t148 + t105) * t92 + t168 * t90 + (-t92 * t142 + t166 + ((t92 * t17 - t90 * t38) * qJD(1) + t1) * qJD(3)) * t93, t39 * t19 + t25 * t45 + (-t105 - (t4 - t148) * t91) * t90 + t168 * t92 + (t90 * t142 + t165 + (-(t90 * t17 + t92 * t38) * qJD(1) - t2) * qJD(3)) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t74, -t156, 0, t74, t18 * t91 - t30 * t93 + (t23 * t91 + t27 * t93) * qJD(3), t74, t156, 0, -t6 * t91 - t8 * t93 + (t10 * t91 - t164) * qJD(3), 0, 0, 0, 0, 0, -t162 + (-t129 - t160) * qJD(3) + t154, t155 - t167; 0, 0, 0, 0, -t127, t149 * t96, 0, 0, 0, -t48 * t140 + t106, -t48 * t146 + t169 - t61 (-t28 * t91 + t46 * t93) * qJD(1) + t106, 0, -t169 + (t28 * t93 + t46 * t91) * qJD(1) + t151, -t30 * pkin(3) + t18 * qJ(4) + t171 * t27 - t23 * t33 - t28 * t46, t54 + (t21 - t42) * qJD(3) + (t135 * t93 - t34 * t91) * qJD(1) + t151, -t22 * qJD(3) + ((t34 - t132) * t93 + t117) * qJD(1) + t30, 0, -t6 * qJ(4) - t10 * t22 - t13 * t34 + t136 * t15 + t8 * t94, t45 * t157 - t163 (-t19 - t161) * t92 + (-t20 - t159) * t90, -t124 + (-t128 + (t45 - t141) * t93) * qJD(1), t126 + (t129 + (-t44 - t138) * t93) * qJD(1), -t56 * t146, -t89 * t20 - t165 - (t92 * t12 - t90 * t22) * t56 + t136 * t44 + (-t11 * t90 - t82 * t157) * qJD(6) + (-t1 * t93 + t102 * t90) * qJD(1), t89 * t19 + t166 + (t90 * t12 + t92 * t22) * t56 + t136 * t45 + (t90 * t82 * t56 - t11 * t92) * qJD(6) + (t102 * t92 + t2 * t93) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, 0, t53, -t27 * qJD(3) + t28 * t140 + t30, t53, t127, 0, t15 * qJD(3) + (t117 - t122) * qJD(1) + t30, 0, 0, 0, 0, 0, -t124 + qJD(3) * t44 + (-t93 * t141 - t128) * qJD(1), t126 + qJD(3) * t45 + (-t93 * t138 + t129) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0.2e1 * t120 (-t84 - t85) * t96 (-t164 + (t10 + t121) * t91) * qJD(1) + t152, 0, 0, 0, 0, 0, -t126 + (-t129 + (-t44 + t138) * t93) * qJD(1), -t124 + (-t128 + (-t45 - t141) * t93) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t44, -t44 ^ 2 + t45 ^ 2, t19 - t161, t20 - t159, t119, t11 * t45 + t170 * t2 - t90 * t8 + t3, t170 * t1 - t11 * t44 - t90 * t4 - t92 * t8;];
tauc_reg  = t9;
