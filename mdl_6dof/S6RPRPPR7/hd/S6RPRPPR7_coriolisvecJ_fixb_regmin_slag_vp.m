% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:31
% EndTime: 2019-03-09 02:57:34
% DurationCPUTime: 1.43s
% Computational Cost: add. (1944->251), mult. (4217->331), div. (0->0), fcn. (2803->6), ass. (0->145)
t105 = cos(qJ(3));
t160 = cos(pkin(9));
t133 = t160 * t105;
t125 = qJD(1) * t133;
t103 = sin(qJ(3));
t159 = sin(pkin(9));
t132 = t159 * t103;
t82 = qJD(1) * t132;
t71 = t125 - t82;
t65 = qJD(6) + t71;
t191 = qJD(6) - t65;
t190 = t65 ^ 2;
t66 = t71 ^ 2;
t111 = t160 * t103 + t159 * t105;
t67 = t111 * qJD(1);
t189 = -t67 ^ 2 - t66;
t147 = qJD(1) * t103;
t106 = -pkin(1) - pkin(7);
t81 = t106 * qJD(1) + qJD(2);
t63 = -qJ(4) * t147 + t103 * t81;
t52 = t159 * t63;
t146 = qJD(1) * t105;
t64 = -qJ(4) * t146 + t105 * t81;
t34 = t160 * t64 - t52;
t151 = -qJD(5) + t34;
t140 = t105 * qJD(4);
t145 = qJD(3) * t103;
t42 = -t81 * t145 + (qJ(4) * t145 - t140) * qJD(1);
t141 = t103 * qJD(4);
t144 = qJD(3) * t105;
t43 = t81 * t144 + (-qJ(4) * t144 - t141) * qJD(1);
t16 = t159 * t42 + t160 * t43;
t55 = qJD(3) * pkin(3) + t64;
t26 = t160 * t55 - t52;
t136 = t160 * t63;
t27 = t159 * t55 + t136;
t130 = qJD(3) * t159;
t131 = qJD(3) * t160;
t69 = t103 * t130 - t105 * t131;
t70 = -t103 * t131 - t105 * t130;
t188 = -t111 * t16 - t26 * t70 + t27 * t69;
t11 = -qJD(3) * qJD(5) - t16;
t118 = qJD(5) - t26;
t21 = -qJD(3) * pkin(4) + t118;
t23 = -qJD(3) * qJ(5) - t27;
t187 = t11 * t111 + t21 * t70 - t23 * t69;
t181 = t67 * pkin(5);
t10 = -t23 - t181;
t33 = t159 * t64 + t136;
t139 = qJD(1) * qJD(3);
t59 = t111 * t139;
t92 = -t160 * pkin(3) - pkin(4);
t86 = -pkin(8) + t92;
t186 = -t86 * t59 + (t10 - t33 + t181) * t65;
t104 = cos(qJ(6));
t102 = sin(qJ(6));
t46 = t104 * qJD(3) + t102 * t67;
t80 = qJD(3) * t125;
t58 = qJD(3) * t82 - t80;
t25 = t46 * qJD(6) + t104 * t58;
t98 = qJD(1) * qJD(2);
t184 = 0.2e1 * t98;
t183 = pkin(4) + pkin(8);
t182 = t58 * pkin(4);
t180 = t71 * pkin(5);
t15 = t159 * t43 - t160 * t42;
t153 = qJ(4) - t106;
t78 = t153 * t103;
t79 = t153 * t105;
t40 = -t159 * t78 + t160 * t79;
t178 = t15 * t40;
t74 = -t132 + t133;
t177 = t15 * t74;
t161 = t103 * pkin(3) + qJ(2);
t126 = -t74 * qJ(5) + t161;
t22 = t111 * t183 + t126;
t174 = t22 * t59;
t142 = t102 * qJD(3);
t44 = -t104 * t67 + t142;
t170 = t44 * t65;
t169 = t46 * t65;
t168 = t46 * t67;
t167 = t67 * t44;
t166 = t74 * t59;
t135 = t105 * t139;
t165 = pkin(3) * t135 + t98;
t164 = t102 * t59;
t163 = t102 * t65;
t51 = t104 * t59;
t143 = qJD(6) * t104;
t24 = -qJD(6) * t142 - t102 * t58 + t67 * t143;
t162 = t24 * t104;
t157 = qJD(6) * t111;
t107 = qJD(3) ^ 2;
t156 = t107 * t103;
t155 = t107 * t105;
t108 = qJD(1) ^ 2;
t154 = t108 * qJ(2);
t152 = pkin(3) * t144 + qJD(2);
t150 = t180 - t151;
t149 = t103 ^ 2 - t105 ^ 2;
t148 = -t107 - t108;
t138 = 0.2e1 * qJD(1);
t77 = pkin(3) * t147 + qJD(1) * qJ(2) + qJD(4);
t137 = t111 * t143;
t134 = pkin(3) * t146 + t67 * qJ(5);
t128 = qJD(6) * t74 + qJD(1);
t127 = t59 * qJ(5) + t165;
t123 = t111 * t24 - t69 * t46;
t122 = -t111 * t59 - t65 * t69;
t60 = t153 * t145 - t140;
t61 = -qJD(3) * t79 - t141;
t29 = t159 * t61 - t160 * t60;
t120 = -t71 * qJ(5) + t77;
t12 = t183 * t67 + t120;
t8 = -t183 * qJD(3) + t118 + t180;
t2 = t102 * t8 + t104 * t12;
t121 = t102 * t12 - t104 * t8;
t119 = -t163 * t65 - t51;
t4 = t58 * pkin(5) - t11;
t117 = t4 + (-qJD(6) * t86 + t183 * t71 + t134) * t65;
t31 = t74 * pkin(5) + t40;
t116 = -t10 * t69 + t111 * t4 + t31 * t59;
t115 = -t71 * qJD(5) + t127;
t28 = t67 * pkin(4) + t120;
t114 = t28 * t71 + t15;
t113 = -t70 * qJ(5) - t74 * qJD(5) + t152;
t112 = -t104 * t190 + t164;
t30 = t159 * t60 + t160 * t61;
t41 = -t159 * t79 - t160 * t78;
t110 = t111 * t58 + t69 * t67 - t70 * t71 + t166;
t109 = t29 * t71 - t30 * t67 - t40 * t59 + t41 * t58 + t177;
t88 = t159 * pkin(3) + qJ(5);
t36 = pkin(4) * t111 + t126;
t35 = t71 * pkin(4) + t134;
t32 = -pkin(5) * t111 + t41;
t20 = -t69 * pkin(4) + t113;
t14 = t69 * pkin(5) + t30;
t13 = t70 * pkin(5) + t29;
t9 = t115 - t182;
t7 = -t183 * t69 + t113;
t6 = -t59 * pkin(5) + t15;
t5 = t104 * t6;
t3 = -t183 * t58 + t115;
t1 = [0, 0, 0, 0, t184, qJ(2) * t184, -0.2e1 * t103 * t135, 0.2e1 * t149 * t139, -t156, -t155, 0, -t106 * t156 + (qJ(2) * t144 + qJD(2) * t103) * t138, -t106 * t155 + (-qJ(2) * t145 + qJD(2) * t105) * t138, t109 + t188, t77 * t152 + t16 * t41 + t165 * t161 - t26 * t29 + t27 * t30 + t178, t109 + t187, t29 * qJD(3) - t111 * t9 - t20 * t67 + t28 * t69 + t36 * t58, t30 * qJD(3) - t20 * t71 - t28 * t70 + t36 * t59 - t9 * t74, -t11 * t41 + t28 * t20 + t21 * t29 - t23 * t30 + t9 * t36 + t178, t102 * t123 + t46 * t137 (t102 * t44 - t104 * t46) * t69 - (t102 * t25 - t162 + (t102 * t46 + t104 * t44) * qJD(6)) * t111, t102 * t122 + t137 * t65 + t24 * t74 + t46 * t70, t104 * t122 - t157 * t163 - t25 * t74 - t44 * t70, t65 * t70 - t166, -t121 * t70 + t14 * t44 + t32 * t25 + t5 * t74 + (-t3 * t74 - t7 * t65 + t174) * t102 + (t13 * t65 - t116) * t104 + ((-t102 * t31 - t104 * t22) * t65 - t2 * t74 + t10 * t102 * t111) * qJD(6), t14 * t46 - t2 * t70 + t32 * t24 + (-(qJD(6) * t31 + t7) * t65 + t174 - (qJD(6) * t8 + t3) * t74 + t10 * t157) * t104 + (-(-qJD(6) * t22 + t13) * t65 - (-qJD(6) * t12 + t6) * t74 + t116) * t102; 0, 0, 0, 0, -t108, -t154, 0, 0, 0, 0, 0, t148 * t103, t148 * t105, t110, -t77 * qJD(1) - t177 - t188, t110, qJD(1) * t67 - t70 * qJD(3), qJD(1) * t71 - t69 * qJD(3), -t28 * qJD(1) - t177 - t187, 0, 0, 0, 0, 0, t74 * t51 + t111 * t25 - t69 * t44 + (t102 * t128 - t104 * t70) * t65, -t74 * t164 + (t102 * t70 + t104 * t128) * t65 + t123; 0, 0, 0, 0, 0, 0, t105 * t108 * t103, -t149 * t108, 0, 0, 0, -t105 * t154, t103 * t154 (t27 - t33) * t71 + (-t26 + t34) * t67 + (t159 * t58 + t160 * t59) * pkin(3), t26 * t33 - t27 * t34 + (-t77 * t146 - t160 * t15 + t159 * t16) * pkin(3), t88 * t58 - t92 * t59 + (-t23 - t33) * t71 + (t21 + t151) * t67, -t33 * qJD(3) + t35 * t67 + t114, -t28 * t67 + t35 * t71 + (0.2e1 * qJD(5) - t34) * qJD(3) + t16, -t11 * t88 + t15 * t92 + t151 * t23 - t21 * t33 - t28 * t35, -t163 * t46 + t162 (-t25 - t169) * t104 + (-t24 + t170) * t102, t119 + t168, t112 - t167, t65 * t67, t117 * t102 + t104 * t186 - t121 * t67 + t150 * t44 + t88 * t25, -t102 * t186 + t117 * t104 + t150 * t46 - t2 * t67 + t88 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t26 * t71 + t27 * t67 + t165, t189, -t80 + (t82 - t71) * qJD(3), qJD(3) * t67 + t59, -t182 - t23 * t67 + (-qJD(5) - t21) * t71 + t127, 0, 0, 0, 0, 0, t112 + t167, t102 * t190 + t168 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 * t67, -t66 - t107, t23 * qJD(3) + t114, 0, 0, 0, 0, 0, -qJD(3) * t44 + t119, -qJD(3) * t46 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, t24 + t170, t169 - t25, -t59, -t10 * t46 - t102 * t3 - t191 * t2 + t5, t10 * t44 - t102 * t6 - t104 * t3 + t191 * t121;];
tauc_reg  = t1;
