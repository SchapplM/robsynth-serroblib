% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:39
% EndTime: 2019-03-08 19:37:46
% DurationCPUTime: 2.51s
% Computational Cost: add. (2821->308), mult. (6968->419), div. (0->0), fcn. (5098->10), ass. (0->185)
t108 = sin(qJ(4));
t171 = qJD(5) * t108;
t174 = qJD(4) * t108;
t104 = sin(pkin(6));
t103 = sin(pkin(11));
t105 = cos(pkin(11));
t109 = sin(qJ(2));
t112 = cos(qJ(2));
t133 = t103 * t112 + t105 * t109;
t59 = t133 * t104;
t53 = qJD(1) * t59;
t227 = pkin(4) * t174 - t171 - t53;
t111 = cos(qJ(4));
t139 = pkin(9) * t108 - qJ(5) * t111;
t121 = t139 * qJD(4);
t226 = t121 + t227;
t180 = qJD(1) * t104;
t161 = t109 * t180;
t160 = t112 * t180;
t78 = qJD(2) * pkin(2) + t160;
t46 = t103 * t78 + t105 * t161;
t42 = qJD(2) * pkin(8) + t46;
t106 = cos(pkin(6));
t90 = qJD(1) * t106 + qJD(3);
t80 = t111 * t90;
t33 = t108 * t42 - t80;
t225 = -qJD(5) - t33;
t107 = sin(qJ(6));
t151 = pkin(5) * qJD(2) + t42;
t195 = t108 * t90;
t58 = (t103 * t109 - t105 * t112) * t104;
t55 = qJD(2) * t58;
t48 = qJD(1) * t55;
t197 = t108 * t48;
t11 = -t197 + (t151 * t111 + t195) * qJD(4);
t110 = cos(qJ(6));
t113 = -pkin(4) - pkin(9);
t182 = t151 * t108 + qJD(5) - t80;
t19 = t113 * qJD(4) + t182;
t152 = -qJ(5) * t108 - pkin(3);
t120 = t113 * t111 + t152;
t79 = t103 * t161;
t45 = t105 * t78 - t79;
t30 = qJD(2) * t120 - t45;
t137 = t107 * t30 - t110 * t19;
t119 = t53 - t171;
t167 = qJD(2) * qJD(4);
t154 = t108 * t167;
t91 = pkin(4) * t154;
t25 = t91 + (t121 + t119) * qJD(2);
t1 = -t137 * qJD(6) + t107 * t11 + t110 * t25;
t178 = qJD(2) * t108;
t92 = qJD(6) + t178;
t224 = t137 * t92 + t1;
t6 = t107 * t19 + t110 * t30;
t2 = -qJD(6) * t6 - t107 * t25 + t110 * t11;
t223 = t6 * t92 + t2;
t28 = -qJD(4) * pkin(4) - t225;
t34 = t111 * t42 + t195;
t29 = -qJD(4) * qJ(5) - t34;
t175 = qJD(4) * t107;
t177 = qJD(2) * t111;
t73 = t110 * t177 + t175;
t131 = t73 * t92;
t49 = qJD(6) * t73 - t107 * t154;
t222 = t49 - t131;
t155 = t107 * t177;
t173 = qJD(4) * t110;
t75 = -t155 + t173;
t204 = t75 * t92;
t169 = qJD(6) * t110;
t50 = qJD(4) * t169 - qJD(6) * t155 - t110 * t154;
t221 = -t50 + t204;
t168 = qJD(6) * t111;
t156 = t110 * t168;
t102 = t111 ^ 2;
t179 = qJD(2) * t102;
t220 = -(-t108 * t92 + t179) * t175 - t92 * t156;
t40 = t106 * t108 + t111 * t59;
t17 = t40 * qJD(4) - t108 * t55;
t54 = qJD(2) * t59;
t219 = qJD(2) * (-t111 * t54 + t58 * t174) - qJD(4) * t17;
t172 = qJD(4) * t111;
t18 = -t59 * t174 + (qJD(4) * t106 - t55) * t111;
t218 = qJD(2) * (t108 * t54 + t58 * t172) - qJD(4) * t18;
t114 = qJD(4) ^ 2;
t94 = pkin(2) * t103 + pkin(8);
t189 = t114 * t94;
t162 = qJ(5) * t172;
t203 = t162 - t227;
t32 = t91 + (t119 - t162) * qJD(2);
t217 = t203 * qJD(2) - t189 - t32;
t140 = -t107 * t137 - t110 * t6;
t216 = -qJD(6) * t140 + t1 * t107 + t110 * t2;
t213 = pkin(5) + t94;
t184 = t108 * t110;
t210 = pkin(2) * t105;
t60 = t120 - t210;
t69 = t213 * t108;
t37 = t107 * t69 + t110 * t60;
t56 = t105 * t160 - t79;
t70 = t213 * t111;
t64 = qJD(4) * t70;
t212 = t37 * qJD(6) + t226 * t107 - t110 * t64 + t56 * t184;
t185 = t107 * t108;
t36 = -t107 * t60 + t110 * t69;
t211 = -t36 * qJD(6) - t107 * t64 - t226 * t110 + t56 * t185;
t166 = t111 * t48 - t90 * t172 + t42 * t174;
t9 = (-pkin(5) * t178 + qJD(5)) * qJD(4) - t166;
t209 = t107 * t9;
t208 = t110 * t9;
t14 = qJD(4) * t34 - t197;
t39 = -t106 * t111 + t108 * t59;
t207 = t14 * t39;
t47 = qJD(1) * t54;
t206 = t47 * t58;
t205 = t75 * t73;
t202 = -t49 * t108 + t75 * t172;
t157 = t107 * t168;
t158 = t108 * t173;
t201 = (t157 + t158) * t92;
t200 = t107 * t50;
t199 = t107 * t75;
t198 = t108 * t14;
t196 = t108 * t50;
t194 = t110 * t49;
t193 = t110 * t73;
t192 = t111 * t14;
t191 = t111 * t73;
t190 = t113 * t92;
t115 = qJD(2) ^ 2;
t186 = t104 * t115;
t96 = pkin(5) * t177;
t24 = -t29 + t96;
t183 = t24 * qJD(4);
t101 = t108 ^ 2;
t181 = t101 - t102;
t170 = qJD(6) * t107;
t165 = t92 * t184;
t164 = t73 * t174;
t163 = t92 * t169;
t159 = t110 * t179;
t153 = t111 * t167;
t129 = -pkin(4) * t111 + t152;
t35 = t129 * qJD(2) - t45;
t150 = t35 * t178 - t197;
t146 = t75 * t158;
t144 = t108 * t153;
t143 = (-t101 - t102) * t56 * qJD(2);
t141 = t107 * t6 - t110 * t137;
t20 = -t107 * t58 + t110 * t39;
t21 = t107 * t39 + t110 * t58;
t136 = t108 * t28 - t111 * t29;
t135 = t108 * t33 + t111 * t34;
t130 = t107 * t92;
t127 = qJD(2) * t53 - t189 - t47;
t68 = t129 - t210;
t126 = qJD(4) * (-qJD(2) * t68 - t35 - t56);
t41 = -qJD(2) * pkin(3) - t45;
t95 = -pkin(3) - t210;
t125 = qJD(4) * (qJD(2) * t95 + t41 + t56);
t122 = t108 * t24 + t113 * t172;
t12 = -qJD(4) * qJD(5) + t166;
t118 = t198 - t111 * t12 + (t108 * t29 + t111 * t28) * qJD(4);
t117 = t198 - t111 * t166 + (-t108 * t34 + t111 * t33) * qJD(4);
t116 = (t108 * t17 + t111 * t18 + (-t108 * t40 + t111 * t39) * qJD(4)) * qJD(2);
t100 = t114 * t111;
t99 = t114 * t108;
t98 = pkin(4) * t178;
t89 = t108 * t115 * t111;
t86 = t110 * t153;
t83 = -0.2e1 * t144;
t82 = 0.2e1 * t144;
t81 = t181 * t115;
t76 = -qJ(5) * t177 + t98;
t66 = -0.2e1 * t181 * t167;
t63 = t213 * t174;
t61 = t139 * qJD(2) + t98;
t57 = t107 * t164;
t27 = t34 + t96;
t16 = t107 * t27 + t110 * t61;
t15 = -t107 * t61 + t110 * t27;
t4 = t20 * qJD(6) + t107 * t17 + t110 * t54;
t3 = -t21 * qJD(6) - t107 * t54 + t110 * t17;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t186, -t112 * t186, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * qJD(2), t55 * qJD(2), 0, -t45 * t54 - t46 * t55 - t48 * t59 + t206, 0, 0, 0, 0, 0, 0, t219, t218, t116, -t166 * t40 + t17 * t33 + t18 * t34 + t41 * t54 + t206 + t207, 0, 0, 0, 0, 0, 0, t116, -t219, -t218, -t12 * t40 + t17 * t28 - t18 * t29 + t32 * t58 + t35 * t54 + t207, 0, 0, 0, 0, 0, 0, t153 * t20 + t18 * t73 + t3 * t92 + t40 * t50, -t153 * t21 + t18 * t75 - t4 * t92 - t40 * t49, t20 * t49 - t21 * t50 - t3 * t75 - t4 * t73, t1 * t21 - t137 * t3 + t18 * t24 + t2 * t20 + t4 * t6 + t40 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t133 * t180 + t53) * qJD(2) (qJD(1) * t58 + t56) * qJD(2), 0, t45 * t53 - t46 * t56 + (-t103 * t48 - t105 * t47) * pkin(2), t82, t66, t100, t83, -t99, 0, t108 * t125 + t111 * t127, -t108 * t127 + t111 * t125, t143 + t117, t117 * t94 - t135 * t56 - t41 * t53 + t47 * t95, 0, -t100, t99, t82, t66, t83, t143 + t118, t108 * t126 - t217 * t111, t217 * t108 + t111 * t126, t118 * t94 - t136 * t56 - t203 * t35 + t32 * t68, -t75 * t156 + (t111 * t49 + t75 * t174) * t107, t146 - t57 + (t200 + t194 + (t193 + t199) * qJD(6)) * t111, t202 + t220, -t73 * t157 + (t111 * t50 - t164) * t110, -t196 + (-t159 - t191) * qJD(4) + t201 (t92 + t178) * t172, t50 * t70 - t63 * t73 - t212 * t92 + (-t173 * t24 + t2) * t108 + (-t24 * t170 + t208 - t56 * t73 + (qJD(2) * t36 - t137) * qJD(4)) * t111, -t49 * t70 - t63 * t75 + t211 * t92 + (t175 * t24 - t1) * t108 + (-t24 * t169 - t209 - t56 * t75 + (-qJD(2) * t37 - t6) * qJD(4)) * t111, t36 * t49 - t37 * t50 + t212 * t75 + t211 * t73 - t140 * t174 + (qJD(6) * t141 - t1 * t110 + t107 * t2) * t111, t1 * t37 + t2 * t36 + t70 * t9 - t211 * t6 + t212 * t137 + (-t111 * t56 - t63) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t100, 0, t135 * qJD(4) - t108 * t166 - t192, 0, 0, 0, 0, 0, 0, 0, t99, t100, t136 * qJD(4) - t108 * t12 - t192, 0, 0, 0, 0, 0, 0, t196 + (-t159 + t191) * qJD(4) + t201, t202 - t220, -t146 - t57 + (t200 - t194 + (t193 - t199) * qJD(6)) * t111 (qJD(4) * t141 + t9) * t108 + (t183 - t216) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t81, 0, t89, 0, 0 (-qJD(2) * t41 + t48) * t108, -qJD(4) * t33 - t41 * t177 + t166, 0, 0, 0, 0, 0, -t89, t81, t89, 0, -t76 * t177 + t150 (0.2e1 * qJD(5) + t33) * qJD(4) + (t108 * t76 + t111 * t35) * qJD(2) - t166, -pkin(4) * t14 - qJ(5) * t12 + t225 * t29 - t28 * t34 - t35 * t76, -t130 * t75 - t194 (-t50 - t204) * t110 + (t49 + t131) * t107, -t92 * t170 + t86 + (-t111 * t75 - t185 * t92) * qJD(2), t110 * t131 + t200, -t163 + (-t165 + (t73 - t175) * t111) * qJD(2), -t92 * t177, qJ(5) * t50 + t209 - t15 * t92 + t182 * t73 + (-t107 * t190 + t110 * t24) * qJD(6) + (t110 * t122 + t111 * t137) * qJD(2), -qJ(5) * t49 + t208 + t16 * t92 + t182 * t75 + (-t107 * t24 - t110 * t190) * qJD(6) + (-t107 * t122 + t111 * t6) * qJD(2), t15 * t75 + t16 * t73 + (-t6 * t178 + t113 * t49 - t2 + (-t113 * t73 - t6) * qJD(6)) * t110 + (-t137 * t178 - t113 * t50 - t1 + (t113 * t75 - t137) * qJD(6)) * t107, qJ(5) * t9 + t216 * t113 + t137 * t15 - t16 * t6 + t182 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t101 * t115 - t114 (t29 + t34) * qJD(4) + t150, 0, 0, 0, 0, 0, 0, -qJD(4) * t73 - t130 * t92 + t86, -t163 - qJD(4) * t75 + (-t107 * t172 - t165) * qJD(2), t221 * t107 + t222 * t110, t224 * t107 + t223 * t110 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t73 ^ 2 + t75 ^ 2, -t222, -t205, t221, t153, -t24 * t75 + t223, t24 * t73 - t224, 0, 0;];
tauc_reg  = t5;
