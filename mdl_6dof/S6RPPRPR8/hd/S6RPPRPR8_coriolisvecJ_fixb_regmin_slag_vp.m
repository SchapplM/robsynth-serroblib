% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:10
% EndTime: 2019-03-09 01:56:14
% DurationCPUTime: 1.34s
% Computational Cost: add. (1700->235), mult. (3824->301), div. (0->0), fcn. (2701->6), ass. (0->128)
t96 = cos(qJ(4));
t139 = qJD(1) * t96;
t92 = cos(pkin(9));
t128 = t92 * t139;
t156 = sin(qJ(4));
t125 = qJD(1) * t156;
t91 = sin(pkin(9));
t78 = t91 * t125;
t63 = -t78 + t128;
t59 = qJD(6) + t63;
t172 = qJD(6) - t59;
t69 = t156 * t92 + t96 * t91;
t140 = qJD(1) * t69;
t171 = qJD(4) * t140;
t93 = -pkin(1) - qJ(3);
t166 = t93 * qJD(1);
t77 = qJD(2) + t166;
t123 = -pkin(7) * qJD(1) + t77;
t55 = t123 * t91;
t56 = t123 * t92;
t144 = -t156 * t55 + t96 * t56;
t170 = qJD(5) - t144;
t169 = t59 ^ 2;
t68 = t156 * t91 - t96 * t92;
t34 = t68 * t171;
t124 = qJD(4) * t156;
t137 = qJD(4) * t96;
t65 = -t92 * t124 - t91 * t137;
t168 = t63 * t65 + t34;
t57 = t65 * qJD(4);
t167 = -qJD(1) * t140 + t57;
t143 = t91 ^ 2 + t92 ^ 2;
t165 = t143 * qJD(3);
t129 = t63 * pkin(5) + t170;
t158 = t140 * pkin(5);
t31 = t156 * t56 + t96 * t55;
t26 = -qJD(4) * qJ(5) - t31;
t13 = -t26 - t158;
t160 = pkin(4) + pkin(8);
t164 = t160 * t171 + (t13 - t31 + t158) * t59;
t94 = sin(qJ(6));
t95 = cos(qJ(6));
t39 = t95 * qJD(4) + t140 * t94;
t76 = qJD(4) * t78;
t54 = qJD(4) * t128 - t76;
t22 = t39 * qJD(6) - t95 * t54;
t157 = -pkin(7) + t93;
t70 = t157 * t91;
t71 = t157 * t92;
t24 = t69 * qJD(3) + t70 * t124 - t71 * t137;
t163 = (t63 + t128) * qJD(4) - t76;
t162 = t140 ^ 2;
t161 = t63 ^ 2;
t89 = qJD(1) * qJD(2);
t127 = 0.2e1 * t89;
t159 = t54 * pkin(4);
t86 = t91 * pkin(3);
t82 = qJ(2) + t86;
t115 = t68 * qJ(5) + t82;
t20 = t160 * t69 + t115;
t155 = t20 * t171;
t131 = t94 * qJD(4);
t135 = qJD(6) * t95;
t21 = -qJD(6) * t131 + t135 * t140 + t94 * t54;
t154 = t21 * t95;
t90 = qJD(1) * qJ(2);
t85 = qJD(3) + t90;
t72 = qJD(1) * t86 + t85;
t110 = -t63 * qJ(5) + t72;
t27 = pkin(4) * t140 + t110;
t153 = t27 * t63;
t37 = -t140 * t95 + t131;
t152 = t37 * t59;
t151 = t39 * t59;
t150 = t39 * t140;
t149 = t140 * t37;
t148 = t63 * t140;
t146 = t69 * t94;
t145 = t94 * t171;
t45 = t95 * t171;
t142 = t140 * qJ(5);
t136 = qJD(6) * t69;
t134 = t24 * qJD(4);
t106 = t156 * t71 + t96 * t70;
t25 = -t68 * qJD(3) + t106 * qJD(4);
t133 = t25 * qJD(4);
t66 = -t91 * t124 + t92 * t137;
t132 = t66 * qJD(4);
t126 = qJD(3) * t139;
t35 = t156 * t70 - t96 * t71;
t122 = qJ(5) * t171 + t89;
t121 = qJD(1) * t143;
t120 = t59 * t94;
t114 = qJD(3) * t125;
t119 = -t92 * t114 - t55 * t124 - t91 * t126 + t56 * t137;
t118 = -qJD(6) * t68 + qJD(1);
t14 = t140 * t160 + t110;
t8 = -t160 * qJD(4) + t129;
t2 = t95 * t14 + t94 * t8;
t113 = t94 * t14 - t95 * t8;
t111 = qJD(1) * t63 + t132;
t108 = -t59 * t120 - t45;
t107 = t69 * t135 + t94 * t66;
t105 = -t63 * qJD(5) + t122;
t9 = -qJD(4) * qJD(5) - t119;
t3 = -t54 * pkin(5) - t9;
t104 = t3 + (t59 * t160 + t142) * t59;
t103 = -t65 * qJ(5) + t68 * qJD(5) + qJD(2);
t28 = -t68 * pkin(5) + t35;
t102 = t13 * t66 + t171 * t28 + t3 * t69;
t101 = -t169 * t95 + t145;
t12 = -t91 * t114 + t56 * t124 + t92 * t126 + t55 * t137;
t23 = -qJD(4) * pkin(4) + t170;
t100 = t12 * t68 - t23 * t65 - t26 * t66 - t9 * t69;
t99 = t31 * qJD(4) - t12;
t98 = qJD(1) ^ 2;
t33 = t63 * pkin(4) + t142;
t32 = t69 * pkin(4) + t115;
t29 = -t69 * pkin(5) + t106;
t19 = t66 * pkin(4) + t103;
t15 = t105 + t159;
t11 = t65 * pkin(5) + t25;
t10 = -t66 * pkin(5) - t24;
t7 = t160 * t66 + t103;
t6 = t160 * t54 + t105;
t5 = -pkin(5) * t171 + t12;
t4 = t95 * t5;
t1 = [0, 0, 0, 0, t127, qJ(2) * t127, t91 * t127, t92 * t127, 0.2e1 * qJD(3) * t121 (t85 + t90) * qJD(2) + (-t77 - t166) * t165, t168, -t140 * t65 + t171 * t69 + t68 * t54 - t63 * t66, t57, -t132, 0, 0.2e1 * t140 * qJD(2) + t82 * t54 + t72 * t66 - t133, t134 - t82 * t171 + t72 * t65 + (-qJD(1) * t68 + t63) * qJD(2), -t106 * t54 + t140 * t24 - t171 * t35 + t25 * t63 - t100, -t140 * t19 - t15 * t69 - t27 * t66 - t32 * t54 + t133, t15 * t68 + t171 * t32 - t19 * t63 - t27 * t65 - t134, -t106 * t9 + t12 * t35 + t15 * t32 + t27 * t19 + t23 * t25 + t26 * t24, t107 * t39 + t21 * t146 (-t37 * t94 + t39 * t95) * t66 + (t154 - t22 * t94 + (-t37 * t95 - t39 * t94) * qJD(6)) * t69, t107 * t59 - t69 * t145 - t21 * t68 + t39 * t65, -t69 * t45 + t22 * t68 - t37 * t65 + (-t94 * t136 + t95 * t66) * t59, t59 * t65 + t34, -t113 * t65 + t10 * t37 + t29 * t22 - t4 * t68 + (-t7 * t59 + t6 * t68 + t155) * t94 + (t11 * t59 - t102) * t95 + ((-t95 * t20 - t94 * t28) * t59 + t2 * t68 + t13 * t146) * qJD(6), t10 * t39 - t2 * t65 + t29 * t21 + (-(qJD(6) * t28 + t7) * t59 + t155 + (qJD(6) * t8 + t6) * t68 + t13 * t136) * t95 + (-(-qJD(6) * t20 + t11) * t59 + (-qJD(6) * t14 + t5) * t68 + t102) * t94; 0, 0, 0, 0, -t98, -t98 * qJ(2), -t98 * t91, -t98 * t92, 0 (-t85 - t165) * qJD(1), 0, 0, 0, 0, 0, t167, -t111, -t140 * t66 - t69 * t54 - t168, -t167, t111, -t27 * qJD(1) + t100, 0, 0, 0, 0, 0, -t68 * t45 + t69 * t22 + t66 * t37 + (t118 * t94 - t65 * t95) * t59, t68 * t145 + t69 * t21 + t66 * t39 + (t118 * t95 + t65 * t94) * t59; 0, 0, 0, 0, 0, 0, 0, 0, -t143 * t98, t77 * t121 + t89, 0, 0, 0, 0, 0, t163, -0.2e1 * t171, -t161 - t162, -t163, 0.2e1 * t171, t159 - t26 * t140 + (-qJD(5) - t23) * t63 + t122, 0, 0, 0, 0, 0, t101 + t149, t169 * t94 + t150 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t161 - t162, 0, t76 + (t63 - t128) * qJD(4), 0, -t72 * t63 + t99, qJD(4) * t144 + t140 * t72 - t119, pkin(4) * t171 - qJ(5) * t54 + (-t26 - t31) * t63 + (t23 - t170) * t140, t140 * t33 + t153 - t99, -t27 * t140 + t33 * t63 + (0.2e1 * qJD(5) - t144) * qJD(4) + t119, -t12 * pkin(4) - t9 * qJ(5) - t170 * t26 - t23 * t31 - t27 * t33, -t39 * t120 + t154 (-t22 - t151) * t95 + (-t21 + t152) * t94, t108 + t150, t101 - t149, t59 * t140, qJ(5) * t22 + t104 * t94 - t113 * t140 + t129 * t37 + t164 * t95, qJ(5) * t21 + t104 * t95 + t129 * t39 - t140 * t2 - t164 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -qJD(4) ^ 2 - t161, t26 * qJD(4) + t12 + t153, 0, 0, 0, 0, 0, -qJD(4) * t37 + t108, -qJD(4) * t39 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t21 + t152, t151 - t22, -t171, -t13 * t39 - t172 * t2 - t94 * t6 + t4, t172 * t113 + t13 * t37 - t94 * t5 - t95 * t6;];
tauc_reg  = t1;
