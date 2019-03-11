% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:39
% EndTime: 2019-03-09 02:18:45
% DurationCPUTime: 1.96s
% Computational Cost: add. (3734->240), mult. (9264->331), div. (0->0), fcn. (7408->10), ass. (0->148)
t116 = cos(qJ(6));
t154 = qJD(6) * t116;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t111 = cos(pkin(11));
t118 = cos(qJ(4));
t163 = t118 * t111;
t109 = sin(pkin(11));
t115 = sin(qJ(4));
t164 = t115 * t109;
t129 = -t163 + t164;
t86 = t129 * qJD(1);
t94 = t118 * t109 + t115 * t111;
t87 = t94 * qJD(1);
t60 = t114 * t87 + t117 * t86;
t205 = t116 * t60;
t213 = t154 + t205;
t108 = qJD(4) + qJD(5);
t166 = t60 * t108;
t156 = qJD(5) * t117;
t157 = qJD(5) * t114;
t149 = qJD(1) * t164;
t153 = t118 * qJD(4);
t158 = qJD(1) * t111;
t97 = t153 * t158;
t81 = -qJD(4) * t149 + t97;
t89 = t94 * qJD(4);
t82 = qJD(1) * t89;
t37 = -t114 * t82 + t117 * t81 - t86 * t156 - t87 * t157;
t212 = t37 + t166;
t202 = -qJD(6) - t60;
t211 = qJD(6) + t202;
t113 = sin(qJ(6));
t132 = -t114 * t86 + t117 * t87;
t155 = qJD(6) * t113;
t23 = t108 * t154 + t116 * t37 - t132 * t155;
t56 = t113 * t108 + t116 * t132;
t24 = qJD(6) * t56 + t113 * t37;
t54 = -t116 * t108 + t113 * t132;
t210 = -t113 * t24 + t23 * t116 - t213 * t54;
t21 = t23 * t113;
t209 = t213 * t56 + t21;
t38 = t132 * qJD(5) + t114 * t81 + t117 * t82;
t34 = t113 * t38;
t57 = t202 * t154;
t177 = t34 - t57;
t183 = t56 * t132;
t208 = -t202 * t205 + t177 - t183;
t105 = t111 * qJD(2);
t101 = sin(pkin(10)) * pkin(1) + qJ(3);
t96 = t101 * qJD(1);
t69 = t105 + (-pkin(7) * qJD(1) - t96) * t109;
t77 = t109 * qJD(2) + t111 * t96;
t70 = pkin(7) * t158 + t77;
t131 = -t115 * t69 - t118 * t70;
t45 = -t86 * pkin(8) - t131;
t173 = t114 * t45;
t194 = -t115 * t70 + t118 * t69;
t44 = -t87 * pkin(8) + t194;
t43 = qJD(4) * pkin(4) + t44;
t16 = t117 * t43 - t173;
t14 = -t108 * pkin(5) - t16;
t207 = t14 * t60;
t95 = -cos(pkin(10)) * pkin(1) - t111 * pkin(3) - pkin(2);
t84 = t95 * qJD(1) + qJD(3);
t64 = t86 * pkin(4) + t84;
t206 = t64 * t60;
t204 = t132 * t60;
t203 = t113 * t202;
t167 = t132 * t108;
t201 = -t38 + t167;
t199 = t132 ^ 2 - t60 ^ 2;
t40 = pkin(5) * t132 + t60 * pkin(9);
t182 = t132 * t54;
t196 = t202 * t132;
t36 = t116 * t38;
t195 = -t155 * t202 - t36;
t193 = qJD(3) * t86;
t169 = t117 * t45;
t17 = t114 * t43 + t169;
t15 = t108 * pkin(9) + t17;
t26 = t60 * pkin(5) - pkin(9) * t132 + t64;
t133 = t113 * t15 - t116 * t26;
t192 = t132 * t133 + t14 * t155;
t31 = -t82 * pkin(8) + t194 * qJD(4) - t193;
t126 = t94 * qJD(3);
t125 = qJD(1) * t126;
t32 = -t81 * pkin(8) + t131 * qJD(4) - t125;
t145 = t114 * t31 - t117 * t32;
t3 = t17 * qJD(5) + t145;
t5 = t113 * t26 + t116 * t15;
t191 = t3 * t113 + t5 * t132 + t14 * t154;
t190 = -t64 * t132 - t145;
t66 = -t114 * t129 + t117 * t94;
t181 = t66 * t38;
t65 = t114 * t94 + t117 * t129;
t88 = t129 * qJD(4);
t46 = -t65 * qJD(5) - t114 * t89 - t117 * t88;
t135 = -t202 * t46 + t181;
t150 = t66 * t155;
t189 = -t116 * t135 - t150 * t202;
t144 = t114 * t32 - t45 * t157;
t2 = (qJD(5) * t43 + t31) * t117 + t144;
t180 = pkin(7) + t101;
t90 = t180 * t109;
t91 = t180 * t111;
t52 = -t94 * pkin(8) - t115 * t91 - t118 * t90;
t130 = t115 * t90 - t118 * t91;
t53 = -pkin(8) * t129 - t130;
t28 = t114 * t52 + t117 * t53;
t68 = pkin(4) * t129 + t95;
t33 = t65 * pkin(5) - t66 * pkin(9) + t68;
t27 = t114 * t53 - t117 * t52;
t124 = -t90 * t153 + qJD(3) * t163 + (-qJD(3) * t109 - qJD(4) * t91) * t115;
t48 = -t89 * pkin(8) + t124;
t120 = t130 * qJD(4) - t126;
t49 = t88 * pkin(8) + t120;
t6 = -t27 * qJD(5) + t114 * t49 + t117 * t48;
t188 = t14 * t46 + (qJD(6) * t33 + t6) * t202 - (qJD(6) * t26 + t2) * t65 - t28 * t38 + t3 * t66;
t187 = t87 * pkin(4);
t185 = t14 * t66;
t184 = t33 * t38;
t47 = t66 * qJD(5) - t114 * t88 + t117 * t89;
t179 = t23 * t65 + t56 * t47;
t175 = t113 * t56;
t168 = t46 * t108;
t162 = t88 * qJD(4);
t159 = t109 ^ 2 + t111 ^ 2;
t148 = -pkin(4) * t108 - t43;
t102 = t114 * pkin(4) + pkin(9);
t140 = qJD(6) * t102 + t187 + t40;
t139 = qJD(1) * t159;
t18 = t114 * t44 + t169;
t137 = pkin(4) * t157 - t18;
t136 = -t65 * t24 - t47 * t54;
t134 = t109 * (-t109 * t96 + t105) - t111 * t77;
t128 = t203 * t60 - t195;
t19 = t117 * t44 - t173;
t122 = -t102 * t38 + t207 - (-pkin(4) * t156 + t19) * t202;
t121 = -t135 * t113 + t66 * t57;
t103 = -t117 * pkin(4) - pkin(5);
t83 = t89 * qJD(4);
t42 = t47 * t108;
t11 = t89 * pkin(4) + t47 * pkin(5) - t46 * pkin(9);
t9 = t82 * pkin(4) + t38 * pkin(5) - t37 * pkin(9);
t8 = t116 * t9;
t7 = t28 * qJD(5) + t114 * t48 - t117 * t49;
t1 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t139 (t101 * t139 - t134) * qJD(3), t81 * t94 - t87 * t88, -t129 * t81 - t94 * t82 + t88 * t86 - t87 * t89, -t162, -t83, 0, t120 * qJD(4) + t95 * t82 + t84 * t89, -t124 * qJD(4) + t95 * t81 - t84 * t88, t132 * t46 + t37 * t66, -t132 * t47 - t37 * t65 - t46 * t60 - t181, t168, -t42, 0, -t7 * t108 + t68 * t38 + t64 * t47 + (t60 * t89 + t65 * t82) * pkin(4), -t6 * t108 + t68 * t37 + t64 * t46 + (t132 * t89 + t66 * t82) * pkin(4), -t56 * t150 + (t23 * t66 + t46 * t56) * t116 (-t116 * t54 - t175) * t46 + (-t21 - t116 * t24 + (t113 * t54 - t116 * t56) * qJD(6)) * t66, t179 - t189, t121 + t136, -t202 * t47 + t38 * t65, t27 * t24 - t133 * t47 + t7 * t54 + t8 * t65 + (-t11 * t202 + t184 + (-t15 * t65 + t202 * t28 + t185) * qJD(6)) * t116 + t188 * t113, t27 * t23 - t5 * t47 + t7 * t56 + ((-qJD(6) * t28 + t11) * t202 - t184 - (-qJD(6) * t15 + t9) * t65 - qJD(6) * t185) * t113 + t188 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t162, 0, 0, 0, 0, 0, -t42, -t168, 0, 0, 0, 0, 0, t121 - t136, t179 + t189; 0, 0, 0, 0, 0, 0, -t159 * qJD(1) ^ 2, t134 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t87 * qJD(4), t97 + (-t86 - t149) * qJD(4), 0, 0, 0, 0, 0, t38 + t167, t37 - t166, 0, 0, 0, 0, 0, t128 - t182, -t116 * t202 ^ 2 - t183 - t34; 0, 0, 0, 0, 0, 0, 0, 0, t87 * t86, -t86 ^ 2 + t87 ^ 2, t97 + (t86 - t149) * qJD(4), 0, 0, -t84 * t87 - t125, t84 * t86 + t193, t204, t199, t212, t201, 0, -t60 * t187 + t18 * t108 + (t114 * t148 - t169) * qJD(5) + t190, -t132 * t187 + t19 * t108 + t206 + (qJD(5) * t148 - t31) * t117 - t144, t209, t175 * t202 + t210, t208, t128 + t182, t196, t103 * t24 + t137 * t54 + (t140 * t202 - t3) * t116 + t122 * t113 + t192, t103 * t23 + t116 * t122 + t137 * t56 - t140 * t203 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t199, t212, t201, 0 (-qJD(5) + t108) * t17 + t190, t16 * t108 - t2 + t206, t209, t203 * t56 + t210, t208, -t202 * t203 + t182 + t36, t196, -pkin(5) * t24 - t3 * t116 + (-t113 * t16 + t116 * t40) * t202 - t17 * t54 + t113 * t207 - t177 * pkin(9) + t192, -pkin(5) * t23 - (t113 * t40 + t116 * t16) * t202 - t17 * t56 + t14 * t205 + t195 * pkin(9) + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, -t202 * t54 + t23, -t202 * t56 - t24, t38, -t113 * t2 - t14 * t56 - t211 * t5 + t8, -t113 * t9 - t116 * t2 + t211 * t133 + t14 * t54;];
tauc_reg  = t1;
