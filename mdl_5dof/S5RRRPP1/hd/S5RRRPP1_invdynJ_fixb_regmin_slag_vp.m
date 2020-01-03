% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:00
% EndTime: 2019-12-31 20:50:04
% DurationCPUTime: 1.48s
% Computational Cost: add. (2566->270), mult. (3853->318), div. (0->0), fcn. (2359->12), ass. (0->161)
t222 = qJ(4) + pkin(7);
t130 = sin(pkin(8));
t131 = cos(pkin(8));
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t169 = qJD(3) * t222;
t151 = -t133 * qJD(4) - t136 * t169;
t137 = cos(qJ(2));
t189 = qJD(1) * t137;
t178 = pkin(1) * t189;
t117 = t136 * qJD(4);
t75 = -t133 * t169 + t117;
t197 = t131 * t136;
t198 = t130 * t133;
t81 = -t197 + t198;
t205 = t130 * t151 + t131 * t75 + t81 * t178;
t129 = qJ(1) + qJ(2);
t119 = cos(t129);
t106 = g(2) * t119;
t134 = sin(qJ(2));
t210 = t134 * pkin(1);
t179 = qJD(1) * t210;
t207 = t137 * pkin(1);
t191 = -qJD(2) * t179 + qJDD(1) * t207;
t123 = qJDD(1) + qJDD(2);
t212 = t123 * pkin(2);
t221 = -t191 - t212 + t106;
t124 = qJD(1) + qJD(2);
t82 = t130 * t136 + t131 * t133;
t69 = t82 * t124;
t64 = t69 ^ 2;
t177 = t124 * t197;
t67 = t124 * t198 - t177;
t220 = -t67 ^ 2 - t64;
t118 = sin(t129);
t219 = g(1) * t119 + g(2) * t118;
t107 = g(1) * t118;
t218 = t106 - t107;
t168 = t222 * t124 + t179;
t59 = t168 * t136;
t203 = t130 * t59;
t58 = t168 * t133;
t57 = qJD(3) * pkin(3) - t58;
t34 = t131 * t57 - t203;
t32 = -qJD(3) * pkin(4) + qJD(5) - t34;
t53 = t131 * t59;
t35 = t130 * t57 + t53;
t33 = qJD(3) * qJ(5) + t35;
t76 = t82 * qJD(3);
t185 = qJD(3) * t136;
t172 = t131 * t185;
t186 = qJD(3) * t133;
t173 = t130 * t186;
t77 = t172 - t173;
t184 = qJDD(1) * t134;
t187 = qJD(2) * t137;
t72 = t123 * pkin(7) + (qJD(1) * t187 + t184) * pkin(1);
t154 = qJ(4) * t123 + qJD(4) * t124 + t72;
t158 = qJD(3) * t168;
t25 = qJDD(3) * pkin(3) - t154 * t133 - t136 * t158;
t28 = -t133 * t158 + t154 * t136;
t11 = t130 * t25 + t131 * t28;
t181 = qJDD(3) * qJ(5) + t11;
t8 = qJD(3) * qJD(5) + t181;
t10 = -t130 * t28 + t131 * t25;
t176 = -qJDD(5) + t10;
t202 = qJDD(3) * pkin(4);
t9 = -t176 - t202;
t217 = t32 * t77 - t33 * t76 - t8 * t81 + t9 * t82;
t216 = -t10 * t82 - t11 * t81 - t34 * t77 - t35 * t76;
t149 = t82 * t123 - t124 * t173;
t125 = qJ(3) + pkin(8);
t115 = sin(t125);
t116 = cos(t125);
t208 = t136 * pkin(3);
t110 = pkin(2) + t208;
t66 = -t110 * t124 + qJD(4) - t178;
t29 = t67 * pkin(4) - t69 * qJ(5) + t66;
t215 = -g(3) * t116 + t115 * t219 - t29 * t69 + t176;
t214 = pkin(3) * t133;
t213 = g(3) * t136;
t211 = t124 * pkin(2);
t135 = sin(qJ(1));
t209 = t135 * pkin(1);
t206 = t130 * t75 - t131 * t151 - t82 * t178;
t85 = -t178 - t211;
t204 = t136 * t107 + t85 * t186;
t201 = t115 * t119;
t200 = t116 * t119;
t199 = t119 * t222;
t196 = t133 * t123;
t195 = t136 * t123;
t109 = pkin(7) + t210;
t194 = -qJ(4) - t109;
t38 = -t131 * t58 - t203;
t193 = qJD(5) - t38;
t127 = t133 ^ 2;
t190 = -t136 ^ 2 + t127;
t188 = qJD(2) * t134;
t183 = t221 * t133 + t85 * t185;
t91 = t119 * t110;
t182 = pkin(4) * t200 + qJ(5) * t201 + t91;
t180 = pkin(1) * t187;
t113 = pkin(3) * t186;
t175 = t124 * t188;
t174 = t124 * t186;
t170 = t222 * t133;
t167 = t118 * t222 + t91;
t166 = t194 * t133;
t165 = qJD(3) * t194;
t164 = t124 * t179;
t163 = -t191 + t218;
t36 = t76 * pkin(4) - t77 * qJ(5) - t82 * qJD(5) + t113;
t162 = -t36 + t179;
t161 = t130 * t196 - t131 * t195;
t159 = -t116 * pkin(4) - t115 * qJ(5);
t157 = -t118 * t110 + t199;
t40 = t124 * t76 + t161;
t41 = t124 * t172 + t149;
t45 = pkin(3) * t174 - t110 * t123 + qJDD(4) - t191;
t140 = t40 * pkin(4) - t41 * qJ(5) + t45;
t7 = -t69 * qJD(5) + t140;
t156 = -g(2) * t200 + t116 * t107 + t29 * t76 + t7 * t81;
t155 = -g(2) * t201 + t115 * t107 - t29 * t77 - t7 * t82;
t48 = t81 * pkin(4) - t82 * qJ(5) - t110;
t150 = -t124 * t85 + t219 - t72;
t142 = (-qJD(4) - t180) * t133 + t136 * t165;
t50 = t133 * t165 + t136 * t180 + t117;
t26 = t130 * t50 - t131 * t142;
t27 = t130 * t142 + t131 * t50;
t120 = t136 * qJ(4);
t80 = t136 * t109 + t120;
t46 = t130 * t80 - t131 * t166;
t47 = t130 * t166 + t131 * t80;
t148 = t26 * t69 - t27 * t67 - t47 * t40 + t46 * t41 - t219;
t139 = qJD(3) ^ 2;
t147 = pkin(7) * t139 - t164 - t212;
t111 = -pkin(2) - t207;
t146 = pkin(1) * t175 + t109 * t139 + t111 * t123;
t145 = -pkin(7) * qJDD(3) + (t178 - t211) * qJD(3);
t144 = -qJDD(3) * t109 + (t111 * t124 - t180) * qJD(3);
t143 = (-g(1) * (-t110 + t159) - g(2) * t222) * t118;
t96 = t136 * pkin(7) + t120;
t55 = t130 * t96 + t131 * t170;
t56 = -t130 * t170 + t131 * t96;
t141 = -t205 * t67 + t206 * t69 - t56 * t40 + t55 * t41 - t219;
t138 = cos(qJ(1));
t122 = t124 ^ 2;
t121 = t138 * pkin(1);
t114 = pkin(1) * t188;
t104 = -t131 * pkin(3) - pkin(4);
t102 = t130 * pkin(3) + qJ(5);
t93 = qJDD(3) * t136 - t139 * t133;
t92 = qJDD(3) * t133 + t139 * t136;
t73 = t127 * t123 + 0.2e1 * t136 * t174;
t52 = -0.2e1 * t190 * t124 * qJD(3) + 0.2e1 * t133 * t195;
t44 = t48 - t207;
t39 = t69 * pkin(4) + t67 * qJ(5) + t124 * t214;
t37 = -t130 * t58 + t53;
t30 = t114 + t36;
t1 = [qJDD(1), g(1) * t135 - g(2) * t138, g(1) * t138 + g(2) * t135, t123, (t123 * t137 - t175) * pkin(1) - t163, ((-qJDD(1) - t123) * t134 + (-qJD(1) - t124) * t187) * pkin(1) + t219, t73, t52, t92, t93, 0, t144 * t133 + (-t146 - t221) * t136 + t204, t144 * t136 + (t146 - t107) * t133 + t183, t148 + t216, t11 * t47 + t35 * t27 - t10 * t46 - t34 * t26 + t45 * (-t110 - t207) + t66 * (t114 + t113) - g(1) * (t157 - t209) - g(2) * (t121 + t167), -t26 * qJD(3) - t46 * qJDD(3) + t30 * t67 + t44 * t40 + t156, t148 + t217, t27 * qJD(3) + t47 * qJDD(3) - t30 * t69 - t44 * t41 + t155, t8 * t47 + t33 * t27 + t7 * t44 + t29 * t30 + t9 * t46 + t32 * t26 - g(1) * (t199 - t209) - g(2) * (t121 + t182) + t143; 0, 0, 0, t123, -t163 + t164, (-t184 + (-qJD(2) + t124) * t189) * pkin(1) + t219, t73, t52, t92, t93, 0, t145 * t133 + (-t147 - t221) * t136 + t204, t145 * t136 + (t147 - t107) * t133 + t183, t141 + t216, t11 * t56 - t10 * t55 - t45 * t110 - g(1) * t157 - g(2) * t167 + (-t179 + t113) * t66 + t205 * t35 - t206 * t34, -t206 * qJD(3) - t55 * qJDD(3) - t162 * t67 + t48 * t40 + t156, t141 + t217, t205 * qJD(3) + t56 * qJDD(3) + t162 * t69 - t48 * t41 + t155, -g(1) * t199 - g(2) * t182 - t162 * t29 + t205 * t33 + t206 * t32 + t7 * t48 + t9 * t55 + t8 * t56 + t143; 0, 0, 0, 0, 0, 0, -t133 * t122 * t136, t190 * t122, t196, t195, qJDD(3), t150 * t133 - t213, g(3) * t133 + t150 * t136, (t35 - t37) * t69 + (-t34 + t38) * t67 + (-t130 * t40 - t131 * t41) * pkin(3), t34 * t37 - t35 * t38 + (-t213 + t10 * t131 + t11 * t130 + (-t124 * t66 + t219) * t133) * pkin(3), t37 * qJD(3) - t39 * t67 + (pkin(4) - t104) * qJDD(3) + t215, -t102 * t40 + t104 * t41 + (t33 - t37) * t69 + (t32 - t193) * t67, -g(3) * t115 + t102 * qJDD(3) - t29 * t67 + t39 * t69 - t219 * t116 + (0.2e1 * qJD(5) - t38) * qJD(3) + t181, t8 * t102 + t9 * t104 - t29 * t39 - t32 * t37 - g(3) * (-t159 + t208) + t193 * t33 + t219 * (pkin(4) * t115 - qJ(5) * t116 + t214); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t34 * t69 + t35 * t67 + t218 + t45, 0.2e1 * t69 * qJD(3) + t161, t220, (t67 - t177) * qJD(3) - t149, t33 * t67 + (-qJD(5) - t32) * t69 + t140 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67 - qJDD(3), (t67 + t177) * qJD(3) + t149, -t64 - t139, -t33 * qJD(3) - t202 - t215;];
tau_reg = t1;
