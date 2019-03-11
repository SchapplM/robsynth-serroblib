% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:44
% EndTime: 2019-03-08 19:37:47
% DurationCPUTime: 1.45s
% Computational Cost: add. (1193->235), mult. (3088->333), div. (0->0), fcn. (2295->10), ass. (0->145)
t88 = sin(pkin(6));
t155 = qJD(1) * t88;
t96 = cos(qJ(2));
t131 = t96 * t155;
t66 = qJD(2) * pkin(2) + t131;
t93 = sin(qJ(2));
t136 = t93 * t155;
t89 = cos(pkin(11));
t68 = t89 * t136;
t87 = sin(pkin(11));
t35 = t87 * t66 + t68;
t30 = qJD(2) * pkin(8) + t35;
t90 = cos(pkin(6));
t76 = t90 * qJD(1) + qJD(3);
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t159 = -t92 * t30 + t95 * t76;
t184 = -qJD(5) + t159;
t15 = -qJD(4) * pkin(4) - t184;
t94 = cos(qJ(6));
t144 = t94 * qJD(4);
t153 = qJD(2) * t95;
t91 = sin(qJ(6));
t62 = -t91 * t153 + t144;
t46 = (t87 * t93 - t89 * t96) * t88;
t142 = qJD(4) * qJ(5);
t21 = t95 * t30 + t92 * t76;
t16 = -t142 - t21;
t145 = t92 * qJD(5);
t152 = qJD(4) * t92;
t42 = t87 * t131 + t68;
t183 = -pkin(4) * t152 + t145 + t42;
t146 = t92 * qJD(2);
t78 = qJD(6) + t146;
t182 = -qJD(6) + t78;
t137 = pkin(5) * t146;
t143 = t137 - t184;
t149 = qJD(6) * t94;
t132 = t95 * t149;
t147 = t91 * qJD(4);
t86 = t95 ^ 2;
t154 = qJD(2) * t86;
t170 = t78 * t92;
t181 = -(t154 - t170) * t147 - t78 * t132;
t148 = t46 * qJD(4);
t47 = (t87 * t96 + t89 * t93) * t88;
t43 = qJD(2) * t47;
t28 = t47 * t95 + t90 * t92;
t44 = qJD(2) * t46;
t8 = t28 * qJD(4) - t44 * t92;
t180 = (t92 * t148 - t43 * t95) * qJD(2) - t8 * qJD(4);
t9 = -t47 * t152 + (qJD(4) * t90 - t44) * t95;
t179 = (t95 * t148 + t43 * t92) * qJD(2) - t9 * qJD(4);
t129 = t95 * t142;
t161 = t129 + t183;
t79 = t87 * pkin(2) + pkin(8);
t98 = qJD(4) ^ 2;
t168 = t79 * t98;
t101 = qJD(1) * t47 - t145;
t141 = qJD(2) * qJD(4);
t127 = t92 * t141;
t77 = pkin(4) * t127;
t19 = t77 + (t101 - t129) * qJD(2);
t178 = t161 * qJD(2) - t168 - t19;
t97 = -pkin(4) - pkin(9);
t151 = qJD(4) * t95;
t37 = qJD(1) * t44;
t138 = -t76 * t151 + t30 * t152 + t95 * t37;
t3 = (qJD(5) - t137) * qJD(4) - t138;
t177 = t3 * t91;
t176 = t3 * t94;
t175 = t89 * pkin(2);
t174 = pkin(5) + t79;
t60 = t94 * t153 + t147;
t38 = -t60 * qJD(6) + t91 * t127;
t173 = t38 * t94;
t172 = t60 * t78;
t171 = t62 * t78;
t169 = t78 * t97;
t99 = qJD(2) ^ 2;
t167 = t88 * t99;
t39 = t62 * qJD(6) - t94 * t127;
t166 = t92 * t39;
t165 = t95 * t60;
t164 = t98 * t92;
t163 = t62 * t151 + t38 * t92;
t120 = pkin(9) * t92 - qJ(5) * t95;
t103 = t120 * qJD(4);
t162 = -t103 + t183;
t150 = qJD(6) * t91;
t134 = t78 * t150;
t140 = t94 * t170;
t160 = qJD(4) * t140 + t95 * t134;
t85 = t92 ^ 2;
t158 = t85 - t86;
t139 = t92 * t99 * t95;
t7 = t30 * t151 + t76 * t152 - t92 * t37;
t135 = t94 * t154;
t133 = t78 * t149;
t57 = t174 * t95;
t12 = t77 + (t103 + t101) * qJD(2);
t126 = t95 * t141;
t5 = pkin(5) * t126 + t7;
t128 = -t91 * t12 + t94 * t5;
t125 = -t92 * qJ(5) - pkin(3);
t67 = t87 * t136;
t34 = t89 * t66 - t67;
t10 = t97 * qJD(4) + t143;
t102 = t97 * t95 + t125;
t17 = t102 * qJD(2) - t34;
t1 = t94 * t10 - t91 * t17;
t2 = t91 * t10 + t94 * t17;
t119 = t15 * t92 - t16 * t95;
t27 = t47 * t92 - t90 * t95;
t118 = t27 * t94 - t46 * t91;
t117 = t27 * t91 + t46 * t94;
t48 = t102 - t175;
t56 = t174 * t92;
t116 = t94 * t48 + t91 * t56;
t112 = t78 * t91;
t111 = t21 * qJD(4) - t7;
t110 = -t95 * pkin(4) + t125;
t81 = pkin(5) * t153;
t11 = -t16 + t81;
t109 = t11 * t92 + t97 * t151;
t36 = qJD(1) * t43;
t106 = qJD(2) * t42 - t168 - t36;
t22 = t110 * qJD(2) - t34;
t45 = t89 * t131 - t67;
t55 = t110 - t175;
t105 = qJD(4) * (-qJD(2) * t55 - t22 - t45);
t29 = -qJD(2) * pkin(3) - t34;
t104 = qJD(4) * (qJD(2) * (-pkin(3) - t175) + t29 + t45);
t6 = -qJD(4) * qJD(5) + t138;
t100 = -t6 * t95 + t7 * t92 + (t15 * t95 + t16 * t92) * qJD(4);
t84 = t98 * t95;
t83 = pkin(4) * t146;
t73 = t94 * t126;
t63 = -qJ(5) * t153 + t83;
t52 = qJD(4) * t57;
t51 = t174 * t152;
t49 = t120 * qJD(2) + t83;
t18 = t22 * t146;
t14 = t81 + t21;
t4 = [0, 0, -t93 * t167, -t96 * t167, -t34 * t43 - t35 * t44 + t36 * t46 - t37 * t47, 0, 0, 0, 0, 0, t180, t179 (t8 * t92 + t9 * t95 + (t27 * t95 - t28 * t92) * qJD(4)) * qJD(2), -t180, -t179, t15 * t8 - t16 * t9 + t19 * t46 + t22 * t43 + t7 * t27 - t6 * t28, 0, 0, 0, 0, 0 (-qJD(6) * t117 - t43 * t91 + t8 * t94) * t78 + t118 * t126 + t9 * t60 + t28 * t39 -(qJD(6) * t118 + t43 * t94 + t8 * t91) * t78 - t117 * t126 + t9 * t62 + t28 * t38; 0, 0, 0, 0, t34 * t42 - t35 * t45 + (-t36 * t89 - t37 * t87) * pkin(2), 0.2e1 * t92 * t126, -0.2e1 * t158 * t141, t84, -t164, 0, t92 * t104 + t106 * t95, t95 * t104 - t106 * t92 (-t85 - t86) * t45 * qJD(2) + t100, t92 * t105 - t178 * t95, t95 * t105 + t178 * t92, t100 * t79 - t119 * t45 - t161 * t22 + t19 * t55, -t38 * t91 * t95 + (t92 * t147 - t132) * t62 (-t60 * t91 + t62 * t94) * t152 + (-t173 + t39 * t91 + (t60 * t94 + t62 * t91) * qJD(6)) * t95, t163 + t181, -t166 + (-t135 - t165) * qJD(4) + t160 (t78 + t146) * t151, t57 * t39 - t51 * t60 + (-t11 * t144 + t128) * t92 + ((-t45 * t92 + t52) * t94 + t162 * t91) * t78 + (-t116 * t78 - t2 * t92) * qJD(6) + (-t11 * t150 + t176 - t45 * t60 + ((-t91 * t48 + t94 * t56) * qJD(2) + t1) * qJD(4)) * t95, t57 * t38 - t51 * t62 + (-(qJD(6) * t10 + t12) * t92 + (-qJD(6) * t56 + t162) * t78) * t94 + (-(-qJD(6) * t48 + t52) * t78 + (t11 * qJD(4) + qJD(6) * t17 + t45 * t78 - t5) * t92) * t91 + (-t11 * t149 - t177 - t45 * t62 + (-t116 * qJD(2) - t2) * qJD(4)) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, -t84, 0, t164, t84, t119 * qJD(4) - t6 * t92 - t7 * t95, 0, 0, 0, 0, 0, t166 + (-t135 + t165) * qJD(4) + t160, t163 - t181; 0, 0, 0, 0, 0, -t139, t158 * t99, 0, 0, 0, -t29 * t146 + t111, qJD(4) * t159 - t29 * t153 + t138, 0, -t63 * t153 - t111 + t18 (0.2e1 * qJD(5) - t159) * qJD(4) + (t22 * t95 + t63 * t92) * qJD(2) - t138, -t7 * pkin(4) - t6 * qJ(5) - t15 * t21 + t184 * t16 - t22 * t63, -t62 * t112 + t173 (-t39 - t171) * t94 + (-t38 + t172) * t91, -t134 + t73 + (-t91 * t170 - t95 * t62) * qJD(2), -t133 + (-t140 + (t60 - t147) * t95) * qJD(2), -t78 * t153, qJ(5) * t39 + t177 - (t94 * t14 - t91 * t49) * t78 + t143 * t60 + (t11 * t94 - t91 * t169) * qJD(6) + (-t1 * t95 + t109 * t94) * qJD(2), qJ(5) * t38 + t176 + (t91 * t14 + t94 * t49) * t78 + t143 * t62 + (-t11 * t91 - t94 * t169) * qJD(6) + (-t109 * t91 + t2 * t95) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t85 * t99 - t98, t16 * qJD(4) + t18 + t7, 0, 0, 0, 0, 0, -qJD(4) * t60 - t112 * t78 + t73, -t133 - qJD(4) * t62 + (-t95 * t147 - t140) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t38 + t172, -t39 + t171, t126, -t11 * t62 + t182 * t2 + t128, t182 * t1 + t11 * t60 - t94 * t12 - t91 * t5;];
tauc_reg  = t4;
