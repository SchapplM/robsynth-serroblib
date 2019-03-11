% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:57
% EndTime: 2019-03-08 21:11:03
% DurationCPUTime: 1.88s
% Computational Cost: add. (1188->293), mult. (2977->385), div. (0->0), fcn. (1918->8), ass. (0->158)
t94 = sin(qJ(3));
t155 = t94 * qJD(2);
t68 = qJD(6) + t155;
t194 = qJD(6) - t68;
t89 = sin(pkin(6));
t167 = qJD(1) * t89;
t95 = sin(qJ(2));
t134 = t95 * t167;
t51 = qJD(2) * pkin(8) + t134;
t41 = t94 * t51;
t90 = cos(pkin(6));
t166 = qJD(1) * t90;
t97 = cos(qJ(3));
t64 = t97 * t166;
t176 = t41 - t64;
t193 = qJD(4) + t176;
t25 = -qJD(3) * pkin(3) + t193;
t99 = -pkin(3) - pkin(4);
t132 = qJD(3) * t99;
t149 = qJ(5) * qJD(2);
t23 = -t94 * t149 + t176;
t152 = -qJD(4) - t23;
t14 = t132 - t152;
t165 = qJD(2) * t89;
t98 = cos(qJ(2));
t138 = t98 * t165;
t179 = t89 * t95;
t144 = t94 * t179;
t21 = -qJD(3) * t144 + (qJD(3) * t90 + t138) * t97;
t120 = t94 * t138;
t40 = t97 * t179 + t90 * t94;
t22 = t40 * qJD(3) + t120;
t39 = -t90 * t97 + t144;
t192 = (t21 * t97 + t22 * t94 + (t39 * t97 - t40 * t94) * qJD(3)) * qJD(2);
t177 = pkin(8) - qJ(5);
t30 = t94 * t166 + t97 * t51;
t24 = -t97 * t149 + t30;
t86 = qJD(3) * qJ(4);
t17 = -t24 - t86;
t26 = t86 + t30;
t85 = -pkin(9) + t99;
t106 = pkin(5) * t97 + t85 * t94;
t104 = t106 * qJD(3);
t77 = t94 * qJD(4);
t173 = t97 * t86 + t77;
t55 = t177 * t94;
t191 = (-qJD(6) * t55 + t104 + t134 + t173) * t68;
t190 = (-t176 + t41) * qJD(3);
t15 = qJD(3) * pkin(5) - t17;
t115 = pkin(5) * t94 + pkin(9) * t97;
t53 = -t97 * pkin(3) - t94 * qJ(4) - pkin(2);
t43 = t97 * pkin(4) - t53;
t31 = t115 + t43;
t154 = t94 * qJD(5);
t161 = qJD(3) * t97;
t36 = t177 * t161 - t154;
t63 = t98 * t167;
t131 = qJD(1) * t165;
t116 = t98 * t131;
t162 = qJD(3) * t94;
t130 = qJD(1) * t162;
t12 = t94 * t116 + t90 * t130 + t51 * t161;
t148 = qJ(5) * qJD(3);
t133 = t97 * t148;
t8 = (-t133 - t154) * qJD(2) + t12;
t163 = qJD(2) * t97;
t150 = qJ(4) * qJD(2);
t91 = qJD(2) * pkin(2);
t52 = -t63 - t91;
t32 = -pkin(3) * t163 - t94 * t150 + t52;
t19 = pkin(4) * t163 + qJD(5) - t32;
t9 = t115 * qJD(2) + t19;
t189 = -(qJD(6) * t31 + t36) * t68 + (qJD(3) * t15 - qJD(6) * t9 + t68 * t63 - t8) * t94;
t175 = qJD(3) * t64 + t97 * t116;
t84 = qJD(3) * qJD(4);
t10 = -t51 * t162 + t175 + t84;
t160 = qJD(5) * t97;
t147 = qJD(2) * qJD(3);
t129 = t94 * t147;
t66 = qJ(5) * t129;
t3 = qJD(2) * t160 - t10 - t66;
t93 = sin(qJ(6));
t187 = t3 * t93;
t96 = cos(qJ(6));
t186 = t3 * t96;
t185 = t17 * t97;
t158 = qJD(6) * t93;
t136 = t97 * t158;
t153 = t96 * qJD(3);
t27 = qJD(2) * t136 - qJD(6) * t153 + t96 * t129;
t184 = t27 * t93;
t44 = t93 * t163 - t153;
t183 = t44 * t68;
t156 = t93 * qJD(3);
t45 = t96 * t163 + t156;
t182 = t45 * t68;
t181 = t68 * t94;
t180 = t68 * t96;
t178 = t89 * t98;
t128 = t97 * t147;
t174 = qJ(4) * t128 + qJD(2) * t77;
t87 = t94 ^ 2;
t88 = t97 ^ 2;
t172 = t87 - t88;
t171 = t87 + t88;
t100 = qJD(3) ^ 2;
t170 = t100 * t94;
t169 = t100 * t97;
t101 = qJD(2) ^ 2;
t168 = t101 * t89;
t164 = qJD(2) * t95;
t11 = t85 * qJD(3) - t152;
t159 = qJD(6) * t11;
t157 = qJD(6) * t96;
t151 = -qJD(5) - t19;
t146 = t93 * t181;
t145 = t94 * t180;
t143 = 0.2e1 * t84 + t175;
t142 = pkin(3) * t162;
t141 = t95 * t168;
t140 = t94 * t101 * t97;
t139 = t89 * t164;
t137 = t68 * t158;
t135 = t68 * t157;
t18 = (t134 + t142) * qJD(2) - t174;
t127 = -pkin(8) * t100 - t18;
t125 = t52 - t91;
t124 = t151 * t94;
t123 = qJD(2) * t53 + t32;
t119 = 0.2e1 * t128;
t117 = t94 * t132;
t2 = t11 * t96 + t9 * t93;
t114 = t11 * t93 - t96 * t9;
t111 = qJD(2) * t88 - t181;
t110 = t30 * qJD(3) - t12;
t109 = t96 * t178 - t39 * t93;
t108 = t93 * t178 + t39 * t96;
t107 = -t15 * t94 - t85 * t161;
t28 = t45 * qJD(6) - t93 * t129;
t102 = t10 * t97 + t12 * t94 + (t25 * t97 - t26 * t94) * qJD(3);
t6 = -t94 * t141 + (t97 * t138 + t21) * qJD(3);
t5 = -t97 * t141 + (-t22 - t120) * qJD(3);
t92 = qJ(4) + pkin(5);
t72 = t97 * t150;
t65 = -t87 * t101 - t100;
t56 = t177 * t97;
t49 = t95 * t94 * t131;
t47 = t130 * t178;
t46 = pkin(3) * t155 - t72;
t37 = t142 - t173;
t35 = pkin(8) * t162 - t94 * t148 + t160;
t34 = t99 * t155 + t72;
t33 = t117 + t173;
t20 = t106 * qJD(2) + t72;
t13 = (t117 - t134) * qJD(2) + t174;
t7 = (t104 - t134) * qJD(2) + t174;
t4 = t96 * t7;
t1 = [0, 0, -t141, -t98 * t168, 0, 0, 0, 0, 0, t5, -t6, t5, t192, t6, t10 * t40 + t12 * t39 + t26 * t21 + t25 * t22 + (t32 * t164 - t18 * t98) * t89, t6, -t5, -t192, t14 * t22 - t17 * t21 - t3 * t40 + t8 * t39 + (t13 * t98 - t164 * t19) * t89, 0, 0, 0, 0, 0 (-qJD(6) * t108 - t139 * t96 - t22 * t93) * t68 + t109 * t128 - t21 * t44 - t40 * t28 -(qJD(6) * t109 - t139 * t93 + t22 * t96) * t68 - t108 * t128 - t21 * t45 + t40 * t27; 0, 0, 0, 0, t94 * t119, -0.2e1 * t172 * t147, t169, -t170, 0, -pkin(8) * t169 + t125 * t162 + t47, pkin(8) * t170 + (t125 + t63) * t161, t47 + t123 * t162 + ((-t37 + t134) * qJD(2) + t127) * t97, -t171 * t116 + t102, t49 + (-qJD(2) * t37 + t127) * t94 + (-t123 - t63) * t161, t18 * t53 + t32 * t37 + (-t32 * t95 + (-t25 * t94 - t26 * t97) * t98) * t167 + t102 * pkin(8), t49 + (qJD(2) * t33 + t13) * t94 + (-t35 + (qJD(2) * t43 + t19 - t63) * t97) * qJD(3), -t13 * t97 - t47 + (t19 * t94 + t36) * qJD(3) + (t43 * t162 + (-t33 - t134) * t97) * qJD(2), t3 * t97 - t8 * t94 + (-t14 * t97 - t17 * t94) * qJD(3) + (t35 * t97 - t36 * t94 + (-t55 * t97 + t56 * t94) * qJD(3) + t171 * t63) * qJD(2), t13 * t43 + t14 * t36 + t17 * t35 + t19 * t33 - t3 * t56 + t8 * t55 + (t19 * t95 + (-t14 * t94 + t185) * t98) * t167, -t27 * t96 * t97 + (-t153 * t94 - t136) * t45 (t44 * t96 + t45 * t93) * t162 + (t184 - t28 * t96 + (t44 * t93 - t45 * t96) * qJD(6)) * t97, t68 * t136 + t27 * t94 + (-t111 * t96 - t45 * t97) * qJD(3), t97 * t135 + t28 * t94 + (t111 * t93 + t44 * t97) * qJD(3) (t68 + t155) * t161, -t56 * t28 + t35 * t44 + t4 * t94 + (-t159 * t94 + t191) * t96 + t189 * t93 + (t44 * t63 - t15 * t157 + t187 + ((t31 * t96 - t55 * t93) * qJD(2) - t114) * qJD(3)) * t97, t56 * t27 + t35 * t45 + (-(t7 - t159) * t94 - t191) * t93 + t189 * t96 + (t45 * t63 + t15 * t158 + t186 + (-(t31 * t93 + t55 * t96) * qJD(2) - t2) * qJD(3)) * t97; 0, 0, 0, 0, -t140, t172 * t101, 0, 0, 0, -t52 * t155 + t110, -t52 * t163 - t175 + t190 (-t32 * t94 + t46 * t97) * qJD(2) + t110, 0, -t190 + (t32 * t97 + t46 * t94) * qJD(2) + t143, -t12 * pkin(3) + t10 * qJ(4) + t193 * t26 - t25 * t30 - t32 * t46, t66 + (t23 - t41) * qJD(3) + (t151 * t97 - t34 * t94) * qJD(2) + t143, -t24 * qJD(3) + ((t34 - t148) * t97 + t124) * qJD(2) + t12, 0, -t3 * qJ(4) - t14 * t24 + t152 * t17 - t19 * t34 + t8 * t99, t180 * t45 - t184 (-t27 - t183) * t96 + (-t28 - t182) * t93, -t135 + (-t145 + (t45 - t156) * t97) * qJD(2), t137 + (t146 + (-t44 - t153) * t97) * qJD(2), -t68 * t163, -t92 * t28 - t186 - (t20 * t96 - t24 * t93) * t68 + t152 * t44 + (-t15 * t93 - t180 * t85) * qJD(6) + (t107 * t93 + t114 * t97) * qJD(2), t92 * t27 + t187 + (t20 * t93 + t24 * t96) * t68 + t152 * t45 + (t68 * t85 * t93 - t15 * t96) * qJD(6) + (t107 * t96 + t2 * t97) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, 0, t65, -qJD(3) * t26 + t155 * t32 + t12, t65, t140, 0, t17 * qJD(3) + (t124 - t133) * qJD(2) + t12, 0, 0, 0, 0, 0, -t135 + qJD(3) * t44 + (-t156 * t97 - t145) * qJD(2), t137 + qJD(3) * t45 + (-t153 * t97 + t146) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, 0.2e1 * t129, -t171 * t101 (-t134 - t185 + (t14 + t132) * t94) * qJD(2) + t174, 0, 0, 0, 0, 0, -t137 + (-t146 + (-t44 + t153) * t97) * qJD(2), -t135 + (-t145 + (-t45 - t156) * t97) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t44, -t44 ^ 2 + t45 ^ 2, t27 - t183, t28 - t182, t128, t15 * t45 - t194 * t2 - t93 * t8 + t4, t194 * t114 - t15 * t44 - t93 * t7 - t96 * t8;];
tauc_reg  = t1;
