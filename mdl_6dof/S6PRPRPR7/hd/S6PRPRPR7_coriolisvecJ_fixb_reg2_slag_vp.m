% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:47
% EndTime: 2019-03-08 19:53:53
% DurationCPUTime: 2.12s
% Computational Cost: add. (2351->310), mult. (5068->419), div. (0->0), fcn. (3219->8), ass. (0->182)
t94 = cos(pkin(6));
t184 = qJD(1) * t94;
t96 = sin(qJ(4));
t159 = t96 * t184;
t102 = -pkin(2) - pkin(8);
t100 = cos(qJ(2));
t93 = sin(pkin(6));
t185 = qJD(1) * t93;
t149 = t100 * t185;
t125 = qJD(3) - t149;
t47 = t102 * qJD(2) + t125;
t99 = cos(qJ(4));
t33 = -t99 * t47 + t159;
t227 = qJD(5) + t33;
t151 = t99 * t184;
t22 = t151 + (-pkin(5) * qJD(2) + t47) * t96;
t97 = sin(qJ(2));
t153 = t97 * t185;
t137 = qJD(2) * t153;
t62 = t99 * t137;
t13 = qJD(4) * t22 - t62;
t113 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t99;
t165 = qJD(2) * qJD(4);
t147 = t96 * t165;
t148 = t99 * t165;
t194 = pkin(4) * t148 + qJ(5) * t147;
t17 = (t113 + t149) * qJD(2) + t194;
t101 = -pkin(4) - pkin(9);
t166 = qJD(4) * t101;
t182 = qJD(2) * t99;
t162 = pkin(5) * t182;
t171 = t162 + t227;
t18 = t171 + t166;
t145 = -qJ(5) * t99 + qJ(3);
t115 = pkin(9) * t96 + t145;
t183 = qJD(2) * t96;
t134 = pkin(4) * t183 + t153;
t29 = qJD(2) * t115 + t134;
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t5 = t18 * t98 - t29 * t95;
t1 = qJD(6) * t5 + t95 * t13 + t98 * t17;
t84 = qJD(6) + t182;
t226 = -t5 * t84 + t1;
t6 = t18 * t95 + t29 * t98;
t2 = -qJD(6) * t6 + t98 * t13 - t95 * t17;
t225 = t6 * t84 + t2;
t91 = t96 ^ 2;
t92 = t99 ^ 2;
t224 = (t91 + t92) * t137;
t167 = qJD(4) * qJ(5);
t179 = qJD(4) * t99;
t223 = pkin(4) * t179 + t96 * t167;
t25 = -qJD(4) * pkin(4) + t227;
t34 = t47 * t96 + t151;
t26 = -t34 - t167;
t152 = t98 * t183;
t54 = qJD(4) * t95 - t152;
t118 = t54 * t84;
t177 = qJD(6) * t95;
t35 = qJD(4) * t177 - qJD(6) * t152 - t95 * t148;
t222 = t35 - t118;
t180 = qJD(4) * t98;
t56 = t95 * t183 + t180;
t203 = t56 * t84;
t36 = qJD(6) * t56 - t98 * t148;
t221 = -t36 + t203;
t202 = t93 * t97;
t158 = qJD(2) * t202;
t142 = t99 * t158;
t104 = qJD(2) ^ 2;
t186 = t104 * t93;
t154 = t100 * t186;
t190 = t100 * t93;
t161 = t96 * t190;
t28 = -qJD(4) * t161 + t94 * t179 - t142;
t220 = qJD(4) * (-t28 + t142) + t96 * t154;
t139 = t96 * t158;
t45 = t99 * t190 + t94 * t96;
t27 = qJD(4) * t45 - t139;
t219 = qJD(4) * (-t27 + t139) - t99 * t154;
t163 = qJD(4) * t159 - t96 * t137 - t47 * t179;
t16 = t34 * qJD(4) - t62;
t210 = t16 * t99;
t106 = (t33 * t96 + t34 * t99) * qJD(4) - t163 * t96 - t210;
t14 = -qJD(4) * qJD(5) + t163;
t107 = (t25 * t96 - t26 * t99) * qJD(4) - t14 * t96 - t210;
t131 = t5 * t95 - t6 * t98;
t218 = -qJD(6) * t131 + t1 * t95 + t2 * t98;
t201 = t95 * t99;
t90 = t96 * pkin(4);
t49 = t115 + t90;
t196 = pkin(5) - t102;
t65 = t196 * t99;
t23 = -t49 * t95 + t65 * t98;
t37 = t113 + t223;
t146 = qJD(4) * t196;
t50 = t96 * t146;
t215 = t23 * qJD(6) + t98 * t37 - t95 * t50 - (t100 * t98 - t97 * t201) * t185;
t197 = t98 * t99;
t24 = t49 * t98 + t65 * t95;
t214 = t24 * qJD(6) + t95 * t37 + t98 * t50 + (-t100 * t95 - t97 * t197) * t185;
t11 = (qJD(5) - t162) * qJD(4) - t163;
t213 = t11 * t95;
t212 = t11 * t98;
t211 = t16 * t45;
t209 = t35 * t99;
t208 = t36 * t95;
t207 = t36 * t99;
t124 = qJD(3) + t149;
t52 = t124 * qJD(2);
t206 = t52 * t97;
t205 = t54 * t96;
t204 = t56 * t54;
t200 = t96 * t35;
t199 = t96 * t98;
t198 = t98 * t35;
t40 = t145 * qJD(2) + t134;
t173 = t40 * qJD(2);
t195 = t99 * t173 - t62;
t57 = pkin(4) * t182 + qJ(5) * t183;
t193 = t91 - t92;
t192 = qJD(2) * pkin(2);
t169 = qJD(2) * qJ(3);
t58 = t153 + t169;
t191 = t100 * t58;
t189 = t101 * t84;
t103 = qJD(4) ^ 2;
t188 = t103 * t96;
t187 = t103 * t99;
t181 = qJD(4) * t96;
t178 = qJD(5) * t99;
t176 = qJD(6) * t98;
t175 = t102 * t103;
t19 = t22 + t167;
t174 = t19 * qJD(4);
t172 = t58 * qJD(2);
t170 = t103 + t104;
t168 = qJD(2) * t100;
t164 = t84 * t201;
t160 = t97 * t186;
t157 = t84 * t177;
t156 = t96 * t177;
t155 = t84 * t176;
t150 = t93 * t168;
t144 = t84 + t182;
t138 = qJD(6) * t99 + qJD(2);
t136 = t99 * t147;
t135 = -t58 + t153;
t132 = t5 * t98 + t6 * t95;
t123 = t144 * t96;
t122 = -qJD(2) * t91 + t84 * t99;
t121 = t52 * qJ(3) + t58 * qJD(3);
t120 = t84 * t138;
t31 = t98 * t202 + t45 * t95;
t30 = -t95 * t202 + t45 * t98;
t63 = t145 + t90;
t112 = -qJD(2) * t63 + t153 - t40;
t111 = t135 - t169;
t20 = (t124 - t178) * qJD(2) + t194;
t43 = qJD(3) - t178 + t223;
t109 = t175 - t20 + (-t43 + t149) * qJD(2);
t108 = t125 * qJD(2) - t175 + t52;
t46 = t94 * t99 - t161;
t105 = (t27 * t96 + t28 * t99 + (-t45 * t96 - t46 * t99) * qJD(4)) * qJD(2);
t82 = t99 * t104 * t96;
t76 = t95 * t147;
t73 = -0.2e1 * t136;
t72 = 0.2e1 * t136;
t71 = t170 * t99;
t70 = t170 * t96;
t69 = t193 * t104;
t64 = t196 * t96;
t53 = t125 - t192;
t51 = t99 * t146;
t48 = 0.2e1 * t193 * t165;
t44 = pkin(9) * t182 + t57;
t10 = t22 * t95 + t44 * t98;
t9 = t22 * t98 - t44 * t95;
t8 = qJD(6) * t30 + t98 * t150 + t28 * t95;
t7 = -qJD(6) * t31 - t95 * t150 + t28 * t98;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t154, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t154 (t206 + (-t135 * t100 + t53 * t97) * qJD(2)) * t93, 0, 0, 0, 0, 0, 0, t220, -t219, t105, -t163 * t46 + t211 - t27 * t34 + t28 * t33 + (t58 * t168 + t206) * t93, 0, 0, 0, 0, 0, 0, t105, -t220, t219, -t14 * t46 + t211 + t25 * t28 + t26 * t27 + (t40 * t168 + t20 * t97) * t93, 0, 0, 0, 0, 0, 0, -t30 * t147 - t27 * t54 + t36 * t46 + t7 * t84, t31 * t147 - t27 * t56 - t35 * t46 - t8 * t84, t30 * t35 - t31 * t36 - t54 * t8 - t56 * t7, t1 * t31 + t11 * t46 - t19 * t27 + t2 * t30 + t5 * t7 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3) (-t191 + (-t53 - t192) * t97) * t185 + t121, t73, t48, -t188, t72, -t187, 0, t108 * t96 - t111 * t179, t108 * t99 + t111 * t181, -t106 + t224 (-t191 + (t33 * t99 - t34 * t96) * t97) * t185 + t106 * t102 + t121, 0, t188, t187, t73, t48, t72, -t107 + t224, t109 * t96 + t112 * t179, t109 * t99 - t112 * t181, t20 * t63 + t40 * t43 + (-t100 * t40 + (t25 * t99 + t26 * t96) * t97) * t185 + t107 * t102, -t95 * t200 + (t96 * t176 + t95 * t179) * t56 (-t54 * t95 + t56 * t98) * t179 + (-t198 - t208 + (-t54 * t98 - t56 * t95) * qJD(6)) * t96, t96 * t155 - t209 + (t122 * t95 - t56 * t96) * qJD(4), -t36 * t199 + (-t98 * t179 + t156) * t54, -t84 * t156 - t207 + (t122 * t98 + t205) * qJD(4), -qJD(4) * t123, -t36 * t64 - t51 * t54 + (-t98 * t174 + t2) * t99 - t214 * t84 + (-t54 * t153 + t19 * t177 - t212 + (-qJD(2) * t23 - t5) * qJD(4)) * t96, t35 * t64 - t51 * t56 + (t95 * t174 - t1) * t99 - t215 * t84 + (-t56 * t153 + t19 * t176 + t213 + (qJD(2) * t24 + t6) * qJD(4)) * t96, t23 * t35 - t24 * t36 + t214 * t56 - t215 * t54 - t131 * t179 + (-qJD(6) * t132 + t1 * t98 - t2 * t95) * t96, t1 * t24 - t11 * t64 + t2 * t23 + t215 * t6 - t214 * t5 + (-t153 * t96 - t51) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t135 * qJD(2), 0, 0, 0, 0, 0, 0, -t70, -t71, 0, t106 - t172, 0, 0, 0, 0, 0, 0, 0, t70, t71, t107 - t173, 0, 0, 0, 0, 0, 0, t96 * t36 + t95 * t120 + (t144 * t199 + t54 * t99) * qJD(4), -t200 + t98 * t120 + (-t123 * t95 + t56 * t99) * qJD(4) (t138 * t54 - t181 * t56 - t209) * t98 + (-t138 * t56 - t181 * t54 + t207) * t95, t131 * qJD(2) + (qJD(4) * t132 + t11) * t96 + (t174 - t218) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t69, 0, -t82, 0, 0, -t99 * t172 + t62, -qJD(4) * t33 + t96 * t172 + t163, 0, 0, 0, 0, 0, t82, -t69, -t82, 0, t57 * t183 + t195 (0.2e1 * qJD(5) + t33) * qJD(4) + (-t40 * t96 + t57 * t99) * qJD(2) - t163, -pkin(4) * t16 - qJ(5) * t14 - t227 * t26 - t25 * t34 - t40 * t57, -t203 * t95 - t198 (-t36 - t203) * t98 + (t35 + t118) * t95, -t157 + (-t164 + (t56 - t180) * t96) * qJD(2), t118 * t98 + t208, -t155 + t76 + (-t84 * t197 - t205) * qJD(2), t84 * t183, qJ(5) * t36 + t213 - t84 * t9 + t171 * t54 + (-t95 * t189 + t19 * t98) * qJD(6) + (t19 * t197 + (-t98 * t166 + t5) * t96) * qJD(2), -qJ(5) * t35 + t10 * t84 + t212 + t171 * t56 + (-t98 * t189 - t19 * t95) * qJD(6) + (-t6 * t96 + (t96 * t166 - t19 * t99) * t95) * qJD(2), t10 * t54 + t56 * t9 + (-t6 * t182 + t101 * t35 - t2 + (-t101 * t54 - t6) * qJD(6)) * t98 + (t5 * t182 - t101 * t36 - t1 + (t101 * t56 + t5) * qJD(6)) * t95, qJ(5) * t11 - t10 * t6 + t101 * t218 + t171 * t19 - t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t104 * t92 - t103 (t26 + t34) * qJD(4) + t195, 0, 0, 0, 0, 0, 0, -t157 - qJD(4) * t54 + (-t96 * t180 - t164) * qJD(2), -t84 ^ 2 * t98 - qJD(4) * t56 + t76, t221 * t95 + t222 * t98, t225 * t98 + t226 * t95 - t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, -t54 ^ 2 + t56 ^ 2, -t222, -t204, t221, -t147, -t19 * t56 + t225, t19 * t54 - t226, 0, 0;];
tauc_reg  = t3;
