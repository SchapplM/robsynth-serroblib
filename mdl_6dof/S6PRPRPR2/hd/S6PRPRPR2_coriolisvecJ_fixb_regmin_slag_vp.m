% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:54
% EndTime: 2019-03-08 19:33:00
% DurationCPUTime: 2.08s
% Computational Cost: add. (2156->285), mult. (5696->442), div. (0->0), fcn. (4523->12), ass. (0->163)
t128 = cos(pkin(11));
t132 = sin(qJ(2));
t126 = sin(pkin(6));
t185 = qJD(1) * t126;
t173 = t132 * t185;
t110 = t128 * t173;
t125 = sin(pkin(11));
t135 = cos(qJ(2));
t172 = t135 * t185;
t72 = t125 * t172 + t110;
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t159 = pkin(4) * t131 - qJ(5) * t134;
t89 = qJD(4) * t159 - t131 * qJD(5);
t214 = t72 - t89;
t182 = qJD(2) * t134;
t116 = -qJD(6) + t182;
t213 = qJD(6) + t116;
t127 = cos(pkin(12));
t124 = sin(pkin(12));
t193 = t124 * t134;
t109 = t125 * t173;
t75 = t128 * t172 - t109;
t118 = pkin(2) * t125 + pkin(8);
t180 = qJD(4) * t131;
t170 = t118 * t180;
t97 = t124 * t170;
t206 = t127 * t214 - t193 * t75 - t97;
t190 = t127 * t134;
t212 = -t124 * t214 - t190 * t75;
t129 = cos(pkin(6));
t115 = qJD(1) * t129 + qJD(3);
t108 = qJD(2) * pkin(2) + t172;
t68 = t125 * t108 + t110;
t65 = qJD(2) * pkin(8) + t68;
t211 = t115 * t134 - t131 * t65;
t181 = qJD(4) * t124;
t183 = qJD(2) * t131;
t100 = t127 * t183 + t181;
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t176 = t127 * qJD(4);
t98 = t124 * t183 - t176;
t56 = t100 * t133 - t130 * t98;
t210 = qJD(6) * t56;
t81 = (t125 * t132 - t128 * t135) * t126;
t209 = pkin(2) * t128;
t208 = pkin(9) + qJ(5);
t74 = qJD(2) * t81;
t70 = qJD(1) * t74;
t16 = -t134 * t70 + (qJD(5) + t211) * qJD(4);
t82 = (t125 * t135 + t128 * t132) * t126;
t40 = (qJD(1) * t82 + t89) * qJD(2);
t7 = t124 * t40 + t127 * t16;
t150 = pkin(5) * t131 - pkin(9) * t190;
t140 = t150 * qJD(4);
t207 = -t140 + t206;
t191 = t127 * t131;
t205 = (-pkin(9) * t193 - t118 * t191) * qJD(4) + t212;
t164 = t127 * t170;
t204 = t164 - t212;
t44 = t131 * t115 + t134 * t65;
t39 = qJD(4) * qJ(5) + t44;
t149 = -pkin(4) * t134 - qJ(5) * t131 - pkin(3);
t67 = t108 * t128 - t109;
t49 = qJD(2) * t149 - t67;
t11 = t124 * t49 + t127 * t39;
t175 = qJD(2) * qJD(4);
t168 = t131 * t175;
t103 = t124 * t133 + t127 * t130;
t143 = t103 * t134;
t139 = qJD(4) * t143;
t177 = qJD(6) * t133;
t178 = qJD(6) * t131;
t194 = t124 * t130;
t48 = t177 * t191 - t178 * t194 + t139;
t83 = t103 * t131;
t203 = t48 * t116 - t83 * t168;
t105 = t159 * qJD(2);
t21 = t124 * t105 + t127 * t211;
t189 = t133 * t127;
t102 = -t189 + t194;
t142 = t102 * t134;
t202 = qJD(2) * t142 - t102 * qJD(6);
t201 = -qJD(2) * t143 + t103 * qJD(6);
t200 = t116 * t56;
t179 = qJD(4) * t134;
t18 = t115 * t180 - t131 * t70 + t65 * t179;
t199 = t124 * t18;
t198 = t127 * t18;
t165 = -qJD(4) * pkin(4) + qJD(5);
t37 = t165 - t211;
t197 = t131 * t37;
t96 = t149 - t209;
t61 = t118 * t190 + t124 * t96;
t169 = t134 * t175;
t196 = t169 * t189 - t98 * t177;
t137 = qJD(2) ^ 2;
t192 = t126 * t137;
t136 = qJD(4) ^ 2;
t188 = t136 * t131;
t187 = t136 * t134;
t122 = t131 ^ 2;
t186 = -t134 ^ 2 + t122;
t184 = qJD(2) * t122;
t6 = -t124 * t16 + t127 * t40;
t4 = qJD(2) * t140 + t6;
t163 = t124 * t169;
t5 = -pkin(9) * t163 + t7;
t174 = -t130 * t5 + t133 * t4;
t171 = t124 * t182;
t167 = pkin(5) * t124 + t118;
t10 = -t124 * t39 + t127 * t49;
t20 = t127 * t105 - t124 * t211;
t23 = (-qJD(6) * t100 - t163) * t130 + t196;
t166 = -t134 * t23 + t56 * t180;
t162 = -t124 * t6 + t127 * t7;
t161 = t130 * t4 + t133 * t5;
t8 = -pkin(5) * t182 - pkin(9) * t100 + t10;
t9 = -pkin(9) * t98 + t11;
t160 = t130 * t9 - t133 * t8;
t2 = t130 * t8 + t133 * t9;
t158 = -t10 * t124 + t11 * t127;
t63 = t129 * t131 + t134 * t82;
t27 = -t124 * t63 + t127 * t81;
t28 = t124 * t81 + t127 * t63;
t157 = -t130 * t28 + t133 * t27;
t156 = t130 * t27 + t133 * t28;
t86 = t127 * t96;
t46 = -pkin(9) * t191 + t86 + (-t118 * t124 - pkin(5)) * t134;
t51 = -pkin(9) * t124 * t131 + t61;
t155 = -t130 * t51 + t133 * t46;
t154 = t130 * t46 + t133 * t51;
t152 = t129 * t134 - t131 * t82;
t24 = qJD(2) * t139 + t210;
t90 = t133 * t98;
t54 = t100 * t130 + t90;
t148 = t134 * t24 - t180 * t54;
t73 = qJD(2) * t82;
t69 = qJD(1) * t73;
t147 = qJD(2) * t72 - t118 * t136 - t69;
t64 = -qJD(2) * pkin(3) - t67;
t146 = qJD(4) * (qJD(2) * (-pkin(3) - t209) + t64 + t75);
t113 = t208 * t124;
t145 = pkin(9) * t171 + qJD(5) * t127 - qJD(6) * t113 - t21;
t114 = t208 * t127;
t144 = qJD(2) * t150 + qJD(5) * t124 + qJD(6) * t114 + t20;
t47 = -qJD(4) * t142 - t103 * t178;
t84 = t102 * t131;
t141 = t116 * t47 + t168 * t84;
t138 = -qJ(5) * t180 + (t165 - t37) * t134;
t120 = -pkin(5) * t127 - pkin(4);
t88 = t167 * t131;
t80 = t167 * t179;
t60 = -t118 * t193 + t86;
t31 = pkin(5) * t171 + t44;
t26 = qJD(4) * t152 - t74 * t134;
t25 = qJD(4) * t63 - t74 * t131;
t22 = pkin(5) * t98 + t37;
t14 = pkin(5) * t163 + t18;
t13 = t124 * t73 + t127 * t26;
t12 = -t124 * t26 + t127 * t73;
t1 = [0, 0, -t132 * t192, -t135 * t192, -t67 * t73 - t68 * t74 + t69 * t81 - t70 * t82, 0, 0, 0, 0, 0, -t25 * qJD(4) + (-t134 * t73 + t180 * t81) * qJD(2), -t26 * qJD(4) + (t131 * t73 + t179 * t81) * qJD(2), t25 * t98 + (-t12 * t134 + (t131 * t27 - t152 * t193) * qJD(4)) * qJD(2), t100 * t25 + (t13 * t134 + (-t131 * t28 - t152 * t190) * qJD(4)) * qJD(2), -t100 * t12 - t13 * t98 + (-t124 * t28 - t127 * t27) * t169, t10 * t12 + t11 * t13 - t152 * t18 + t25 * t37 + t27 * t6 + t28 * t7, 0, 0, 0, 0, 0 -(-qJD(6) * t156 + t133 * t12 - t130 * t13) * t116 + t157 * t168 + t25 * t54 - t152 * t24 (qJD(6) * t157 + t130 * t12 + t133 * t13) * t116 - t156 * t168 + t25 * t56 - t152 * t23; 0, 0, 0, 0, t67 * t72 - t68 * t75 + (-t125 * t70 - t128 * t69) * pkin(2), 0.2e1 * t134 * t168, -0.2e1 * t186 * t175, t187, -t188, 0, t131 * t146 + t134 * t147, -t131 * t147 + t134 * t146 (t199 - t75 * t98 + (qJD(2) * t60 + t10) * qJD(4)) * t131 + (-t6 + (t118 * t98 + t124 * t37) * qJD(4) + (t97 + t206) * qJD(2)) * t134 (-t100 * t75 + t198 + (-qJD(2) * t61 - t11) * qJD(4)) * t131 + (t7 + (t100 * t118 + t127 * t37) * qJD(4) + (t164 - t204) * qJD(2)) * t134, t204 * t98 + (-t124 * t7 - t127 * t6) * t131 + t206 * t100 + (-t10 * t127 - t11 * t124 + (-t124 * t61 - t127 * t60) * qJD(2)) * t179, -t75 * t197 + t6 * t60 + t61 * t7 + (t131 * t18 + t179 * t37) * t118 - t204 * t11 - t206 * t10, -t23 * t84 + t47 * t56, -t23 * t83 + t24 * t84 - t47 * t54 - t48 * t56, -t141 + t166, t148 + t203 (-t116 - t182) * t180, -t174 * t134 + t80 * t54 + t88 * t24 + t14 * t83 + t22 * t48 + (t205 * t130 + t207 * t133) * t116 + (t116 * t154 + t134 * t2) * qJD(6) + (-t75 * t54 + (qJD(2) * t155 - t160) * qJD(4)) * t131, t161 * t134 + t80 * t56 + t88 * t23 - t14 * t84 + t22 * t47 + (-t207 * t130 + t205 * t133) * t116 + (t116 * t155 - t134 * t160) * qJD(6) + (-t75 * t56 + (-qJD(2) * t154 - t2) * qJD(4)) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, -t187 (-t124 * t184 + t131 * t98) * qJD(4) (t100 * t131 - t127 * t184) * qJD(4) (t100 * t124 - t127 * t98) * t179, -t134 * t18 + t162 * t131 + (t134 * t158 + t197) * qJD(4), 0, 0, 0, 0, 0, -t148 + t203, t141 + t166; 0, 0, 0, 0, 0, -t131 * t137 * t134, t186 * t137, 0, 0, 0, qJD(4) * t44 - t183 * t64 - t18 (-qJD(2) * t64 + t70) * t134, -t198 - t44 * t98 + (-t10 * t131 + t124 * t138 + t134 * t20) * qJD(2), -t100 * t44 + t199 + (t11 * t131 + t127 * t138 - t134 * t21) * qJD(2), t100 * t20 + t21 * t98 + (-qJD(5) * t98 + t10 * t182 + t7) * t127 + (qJD(5) * t100 + t11 * t182 - t6) * t124, -pkin(4) * t18 + qJ(5) * t162 + qJD(5) * t158 - t10 * t20 - t11 * t21 - t37 * t44, t23 * t103 + t202 * t56, -t102 * t23 - t103 * t24 - t201 * t56 - t202 * t54, -t202 * t116 + (qJD(4) * t103 - t56) * t183, t201 * t116 + (-qJD(4) * t102 + t54) * t183, t116 * t183, t14 * t102 + t120 * t24 - t31 * t54 + t201 * t22 + (t130 * t145 + t133 * t144) * t116 + ((-t113 * t133 - t114 * t130) * qJD(4) + t160) * t183, t14 * t103 + t120 * t23 - t31 * t56 + t202 * t22 + (-t130 * t144 + t133 * t145) * t116 + (-(-t113 * t130 + t114 * t133) * qJD(4) + t2) * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t100 + t181) * t182 (t98 + t176) * t182, -t100 ^ 2 - t98 ^ 2, t10 * t100 + t11 * t98 + t18, 0, 0, 0, 0, 0, t24 - t200, t90 * t116 + (-t163 + (-qJD(6) + t116) * t100) * t130 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, -t116 * t54 + t23, -t103 * t169 - t200 - t210, t168, -t2 * t213 - t22 * t56 + t174, t160 * t213 + t22 * t54 - t161;];
tauc_reg  = t1;
