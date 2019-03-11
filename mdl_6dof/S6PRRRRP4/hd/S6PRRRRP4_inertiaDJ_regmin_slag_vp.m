% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:52
% EndTime: 2019-03-09 00:17:58
% DurationCPUTime: 2.41s
% Computational Cost: add. (2770->293), mult. (7243->477), div. (0->0), fcn. (6626->10), ass. (0->142)
t95 = cos(qJ(3));
t161 = qJD(3) * t95;
t92 = sin(qJ(4));
t142 = t92 * t161;
t94 = cos(qJ(4));
t159 = qJD(4) * t94;
t93 = sin(qJ(3));
t187 = t93 * t159 + t142;
t186 = -0.4e1 * t93;
t176 = cos(qJ(5));
t172 = t93 * t94;
t177 = pkin(8) * t92;
t126 = -t95 * pkin(3) - t93 * pkin(9);
t61 = -pkin(2) + t126;
t56 = t94 * t61;
t35 = -pkin(10) * t172 + t56 + (-pkin(4) - t177) * t95;
t171 = t94 * t95;
t76 = pkin(8) * t171;
t167 = t92 * t61 + t76;
t173 = t92 * t93;
t41 = -pkin(10) * t173 + t167;
t91 = sin(qJ(5));
t185 = t176 * t41 + t91 * t35;
t163 = sin(pkin(6));
t124 = t163 * cos(qJ(2));
t117 = qJD(2) * t124;
t123 = t163 * sin(qJ(2));
t164 = cos(pkin(6));
t47 = t93 * t123 - t164 * t95;
t180 = t47 * qJD(3) - t95 * t117;
t184 = -qJD(4) * t124 - t180;
t106 = t95 * t123 + t164 * t93;
t116 = qJD(2) * t123;
t183 = t106 * qJD(4) - t116;
t133 = t176 * qJD(5);
t182 = t176 * qJD(4) + t133;
t89 = t94 ^ 2;
t166 = t92 ^ 2 - t89;
t132 = t166 * qJD(4);
t181 = qJD(4) + qJD(5);
t158 = qJD(4) * t95;
t148 = t92 * t158;
t162 = qJD(3) * t94;
t110 = t93 * t162 + t148;
t125 = pkin(3) * t93 - pkin(9) * t95;
t59 = t125 * qJD(3);
t29 = t110 * pkin(8) - t61 * t159 - t92 * t59;
t85 = qJD(3) * t93;
t143 = t92 * t85;
t168 = pkin(8) * t143 + t94 * t59;
t178 = pkin(4) * t93;
t18 = (-pkin(10) * t171 + t178) * qJD(3) + (-t76 + (pkin(10) * t93 - t61) * t92) * qJD(4) + t168;
t22 = -pkin(10) * t187 - t29;
t6 = -qJD(5) * t185 + t176 * t18 - t91 * t22;
t96 = 2 * qJD(6);
t179 = -pkin(10) - pkin(9);
t175 = t91 * t92;
t174 = t91 * t94;
t60 = pkin(4) * t173 + t93 * pkin(8);
t88 = t93 ^ 2;
t165 = -t95 ^ 2 + t88;
t160 = qJD(4) * t92;
t157 = qJD(5) * t91;
t156 = t95 * qJD(6);
t155 = -0.2e1 * pkin(2) * qJD(3);
t154 = -0.2e1 * pkin(3) * qJD(4);
t153 = t91 * t173;
t83 = pkin(8) * t161;
t44 = t187 * pkin(4) + t83;
t152 = pkin(5) * t85;
t151 = pkin(4) * t160;
t150 = pkin(4) * t157;
t149 = t93 * t160;
t146 = t94 * t158;
t145 = t47 * t160;
t144 = t47 * t159;
t141 = t92 * t159;
t140 = t93 * t161;
t139 = t94 * t161;
t82 = -t94 * pkin(4) - pkin(3);
t138 = t176 * t94;
t137 = t179 * qJD(4);
t135 = qJD(3) * t176;
t131 = t165 * qJD(3);
t130 = 0.2e1 * t140;
t129 = t92 * t139;
t128 = t179 * t176;
t127 = t95 * t135;
t37 = t175 * t181 - t182 * t94;
t39 = t106 * qJD(3) + t93 * t117;
t58 = t176 * t92 + t174;
t122 = t47 * t37 - t39 * t58;
t121 = t92 * t128;
t120 = qJD(4) * t128;
t119 = qJD(3) * t124;
t65 = t179 * t94;
t25 = -qJD(5) * t121 - t92 * t120 - t137 * t174 - t65 * t157;
t43 = t179 * t175 - t176 * t65;
t115 = t25 * t95 + t43 * t85;
t26 = -t65 * t133 - t94 * t120 + (qJD(5) * t179 + t137) * t175;
t42 = -t91 * t65 - t121;
t114 = t26 * t95 - t42 * t85;
t113 = t176 * t35 - t91 * t41;
t5 = -t35 * t133 + t41 * t157 - t176 * t22 - t91 * t18;
t79 = qJ(6) * t85;
t108 = -t5 + t79;
t102 = -t106 * t92 - t94 * t124;
t101 = t91 * t102;
t40 = t106 * t94 - t92 * t124;
t97 = -t183 * t92 + t184 * t94;
t98 = t183 * t94 + t184 * t92;
t2 = qJD(5) * t101 + t40 * t133 + t176 * t98 + t91 * t97;
t99 = t176 * t102;
t20 = t91 * t40 - t99;
t24 = t92 * t127 - t91 * t149 - qJD(5) * t153 + (t91 * t161 + t182 * t93) * t94;
t45 = t58 * t93;
t107 = t2 * t95 - t20 * t85 + t47 * t24 + t39 * t45;
t1 = -qJD(5) * t99 + t40 * t157 - t176 * t97 + t91 * t98;
t21 = t176 * t40 + t101;
t38 = t181 * t58;
t23 = -t94 * t127 + t91 * t142 + t38 * t93;
t46 = t93 * t138 - t153;
t105 = t1 * t95 + t21 * t85 + t47 * t23 - t39 * t46;
t103 = t95 * t150 + t6;
t84 = pkin(4) * t133;
t81 = -t176 * pkin(4) - pkin(5);
t78 = t91 * pkin(4) + qJ(6);
t77 = -0.2e1 * t150;
t72 = t84 + qJD(6);
t71 = -0.2e1 * t140;
t57 = -t138 + t175;
t33 = t57 * pkin(5) - t58 * qJ(6) + t82;
t30 = -qJD(4) * t167 + t168;
t27 = t45 * pkin(5) - t46 * qJ(6) + t60;
t17 = t95 * pkin(5) - t113;
t16 = -t95 * qJ(6) + t185;
t9 = t38 * pkin(5) + t37 * qJ(6) - t58 * qJD(6) + t151;
t8 = t47 * t38 + t39 * t57;
t7 = t24 * pkin(5) + t23 * qJ(6) - t46 * qJD(6) + t44;
t4 = -t152 - t6;
t3 = t108 - t156;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21 * t1 + 0.2e1 * t20 * t2 + 0.2e1 * t47 * t39; 0, 0, -t116, -t117, 0, 0, 0, 0, 0, -t95 * t116 - t93 * t119, t93 * t116 - t95 * t119, 0, 0, 0, 0, 0, t102 * t85 + t142 * t47 + t144 * t93 + t173 * t39 + t95 * t98, t139 * t47 - t145 * t93 + t172 * t39 - t40 * t85 + t95 * t97, 0, 0, 0, 0, 0, t107, -t105, t107, t1 * t45 + t2 * t46 - t20 * t23 - t21 * t24, t105, -t1 * t16 + t2 * t17 + t20 * t4 + t21 * t3 + t39 * t27 + t47 * t7; 0, 0, 0, 0, t130, -0.2e1 * t131, 0, 0, 0, t93 * t155, t95 * t155, 0.2e1 * t89 * t140 - 0.2e1 * t88 * t141, t129 * t186 + 0.2e1 * t132 * t88, 0.2e1 * t93 * t148 + 0.2e1 * t165 * t162, -0.2e1 * t131 * t92 + 0.2e1 * t146 * t93, t71, 0.2e1 * t56 * t85 - 0.2e1 * t30 * t95 + 0.2e1 * (t140 * t92 + t159 * t88) * pkin(8), -0.2e1 * t29 * t95 - 0.2e1 * t167 * t85 + 0.2e1 * (t130 * t94 - t160 * t88) * pkin(8), -0.2e1 * t46 * t23, 0.2e1 * t23 * t45 - 0.2e1 * t46 * t24, 0.2e1 * t23 * t95 + 0.2e1 * t46 * t85, 0.2e1 * t24 * t95 - 0.2e1 * t45 * t85, t71, 0.2e1 * t113 * t85 + 0.2e1 * t60 * t24 + 0.2e1 * t44 * t45 - 0.2e1 * t6 * t95, -0.2e1 * t185 * t85 - 0.2e1 * t60 * t23 + 0.2e1 * t44 * t46 - 0.2e1 * t5 * t95, -0.2e1 * t17 * t85 + 0.2e1 * t27 * t24 + 0.2e1 * t4 * t95 + 0.2e1 * t7 * t45, -0.2e1 * t16 * t24 - 0.2e1 * t17 * t23 - 0.2e1 * t3 * t45 + 0.2e1 * t4 * t46, 0.2e1 * t16 * t85 + 0.2e1 * t27 * t23 - 0.2e1 * t3 * t95 - 0.2e1 * t7 * t46, 0.2e1 * t16 * t3 + 0.2e1 * t17 * t4 + 0.2e1 * t27 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t180, 0, 0, 0, 0, 0, -t39 * t94 + t145, t39 * t92 + t144, 0, 0, 0, 0, 0, t8, -t122, t8, t1 * t57 + t2 * t58 - t20 * t37 - t21 * t38, t122, -t1 * t43 + t2 * t42 + t20 * t26 - t21 * t25 + t39 * t33 + t47 * t9; 0, 0, 0, 0, 0, 0, t161, -t85, 0, -t83, pkin(8) * t85, -t93 * t132 + t129, t141 * t186 - t166 * t161, t143 - t146, t110, 0 (pkin(9) * t171 + (-pkin(3) * t94 + t177) * t93) * qJD(4) + (t126 * t92 - t76) * qJD(3) (pkin(8) * t172 + t125 * t92) * qJD(4) + (t126 * t94 + t177 * t95) * qJD(3), -t23 * t58 - t46 * t37, t23 * t57 - t58 * t24 + t37 * t45 - t46 * t38, t37 * t95 + t58 * t85, t38 * t95 - t57 * t85, 0, t151 * t45 + t82 * t24 + t60 * t38 + t44 * t57 + t114, t151 * t46 - t82 * t23 - t60 * t37 + t44 * t58 - t115, t33 * t24 + t27 * t38 + t9 * t45 + t7 * t57 + t114, -t16 * t38 - t17 * t37 - t42 * t23 - t43 * t24 + t25 * t45 + t26 * t46 - t3 * t57 + t4 * t58, t33 * t23 + t27 * t37 - t9 * t46 - t7 * t58 + t115, -t16 * t25 + t17 * t26 + t27 * t9 + t3 * t43 + t7 * t33 + t4 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t141, -0.2e1 * t132, 0, 0, 0, t92 * t154, t94 * t154, -0.2e1 * t58 * t37, 0.2e1 * t37 * t57 - 0.2e1 * t58 * t38, 0, 0, 0, 0.2e1 * t151 * t57 + 0.2e1 * t82 * t38, 0.2e1 * t151 * t58 - 0.2e1 * t82 * t37, 0.2e1 * t33 * t38 + 0.2e1 * t9 * t57, 0.2e1 * t25 * t57 + 0.2e1 * t26 * t58 - 0.2e1 * t42 * t37 - 0.2e1 * t43 * t38, 0.2e1 * t33 * t37 - 0.2e1 * t9 * t58, -0.2e1 * t43 * t25 + 0.2e1 * t42 * t26 + 0.2e1 * t33 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -t1 * t78 + t150 * t20 + t2 * t81 + t21 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139 - t149, -t187, t85, t30, t29, 0, 0, -t23, -t24, t85, t135 * t178 + t103 (t133 * t95 - t85 * t91) * pkin(4) + t5 (pkin(5) - t81) * t85 + t103, t150 * t46 - t81 * t23 - t78 * t24 - t72 * t45, t78 * t85 + (-qJD(6) - t72) * t95 + t108, t150 * t17 + t16 * t72 + t3 * t78 + t4 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t160, 0, -pkin(9) * t159, pkin(9) * t160, 0, 0, -t37, -t38, 0, -t26, t25, -t26, t150 * t58 - t81 * t37 - t78 * t38 - t72 * t57, -t25, t150 * t42 - t25 * t78 + t26 * t81 + t43 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t84, t77, 0, 0.2e1 * t72, 0.2e1 * t150 * t81 + 0.2e1 * t78 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -t2 * pkin(5) - t1 * qJ(6) + t21 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, t85, t6, t5, t6 + 0.2e1 * t152, pkin(5) * t23 - t24 * qJ(6) - t45 * qJD(6), -t5 + 0.2e1 * t79 - 0.2e1 * t156, -t4 * pkin(5) + t3 * qJ(6) + t16 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, -t26, t25, -t26, pkin(5) * t37 - t38 * qJ(6) - t57 * qJD(6), -t25, -t26 * pkin(5) - t25 * qJ(6) + t43 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t84, -t150, 0, t96 + t84, -pkin(5) * t150 + t72 * qJ(6) + t78 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, qJ(6) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t23, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
