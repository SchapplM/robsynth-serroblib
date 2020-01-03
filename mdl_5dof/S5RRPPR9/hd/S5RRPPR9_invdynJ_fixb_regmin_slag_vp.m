% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:50
% EndTime: 2019-12-31 19:41:55
% DurationCPUTime: 2.02s
% Computational Cost: add. (1088->310), mult. (2263->373), div. (0->0), fcn. (1258->6), ass. (0->175)
t209 = pkin(2) + pkin(3);
t154 = t209 * qJD(2);
t104 = sin(qJ(2));
t175 = qJD(1) * t104;
t74 = pkin(6) * t175;
t38 = -qJ(4) * t175 + t74;
t158 = qJD(3) + t38;
t21 = -t154 + t158;
t219 = t104 * t209;
t107 = cos(qJ(2));
t166 = qJD(1) * qJD(2);
t149 = t107 * t166;
t165 = t104 * qJDD(1);
t218 = t149 + t165;
t150 = t104 * t166;
t164 = t107 * qJDD(1);
t217 = -t150 + t164;
t174 = qJD(1) * t107;
t151 = t209 * qJDD(2);
t57 = qJD(5) + t175;
t216 = t57 ^ 2;
t201 = t107 * pkin(2) + t104 * qJ(3);
t215 = -pkin(1) - t201;
t203 = pkin(6) - qJ(4);
t75 = pkin(6) * t174;
t40 = -qJ(4) * t174 + t75;
t97 = qJD(2) * qJ(3);
t28 = -t40 - t97;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t214 = g(1) * t108 + g(2) * t105;
t183 = t104 * t108;
t184 = t104 * t105;
t90 = g(3) * t107;
t213 = -g(1) * t183 - g(2) * t184 + t90;
t56 = pkin(6) * t149;
t70 = pkin(6) * t165;
t156 = qJDD(3) + t56 + t70;
t170 = qJD(4) * t104;
t112 = -qJ(4) * t218 - qJD(1) * t170 + t156;
t135 = pkin(4) * t104 + pkin(7) * t107;
t29 = -qJD(1) * pkin(1) - pkin(2) * t174 - qJ(3) * t175;
t19 = pkin(3) * t174 + qJD(4) - t29;
t14 = qJD(1) * t135 + t19;
t96 = -pkin(7) - t209;
t148 = -qJD(5) * t14 - t96 * qJDD(2) - t112;
t86 = t107 * pkin(3);
t161 = t86 + t201;
t34 = pkin(1) + t161;
t16 = t135 + t34;
t22 = qJD(2) * pkin(4) - t28;
t171 = qJD(2) * t107;
t26 = t203 * t171 - t170;
t35 = qJDD(5) + t218;
t45 = t203 * t104;
t172 = qJD(2) * t104;
t125 = pkin(6) * t172 + qJD(4) * t107;
t71 = pkin(6) * t164;
t94 = qJDD(2) * qJ(3);
t95 = qJD(2) * qJD(3);
t162 = t71 + t94 + t95;
t51 = qJ(4) * t150;
t9 = qJ(4) * t164 + qJD(1) * t125 - t162 - t51;
t7 = qJDD(2) * pkin(4) - t9;
t212 = -(qJD(5) * t16 + t26) * t57 + (t22 * qJD(2) + t148) * t104 - t45 * t35 - t7 * t107;
t190 = pkin(6) * qJDD(2);
t211 = (qJD(1) * t215 + t29) * qJD(2) - t190;
t124 = pkin(4) * t107 + t96 * t104;
t207 = g(3) * t104;
t67 = qJ(3) * t174;
t210 = t107 * t214 + (qJD(1) * t124 + qJD(5) * t96 + t67) * t57 - t7 + t207;
t92 = g(1) * t105;
t208 = g(2) * t108;
t206 = t16 * t35;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t167 = t106 * qJD(2);
t36 = t103 * t174 - t167;
t205 = t36 * t57;
t173 = qJD(2) * t103;
t37 = t106 * t174 + t173;
t204 = t37 * t57;
t79 = t104 * qJD(3);
t202 = qJ(3) * t171 + t79;
t98 = t104 ^ 2;
t99 = t107 ^ 2;
t199 = t98 - t99;
t198 = t98 + t99;
t168 = qJD(5) * t107;
t153 = t103 * t168;
t10 = qJD(1) * t153 - qJDD(2) * t103 + (-qJD(2) * qJD(5) - t217) * t106;
t197 = t10 * t103;
t196 = t103 * t35;
t195 = t103 * t57;
t194 = t106 * t35;
t193 = t106 * t57;
t192 = t107 * t36;
t191 = t107 * t37;
t189 = qJ(3) * t107;
t188 = qJD(2) * t28;
t187 = qJD(2) * t36;
t186 = qJD(2) * t37;
t100 = qJDD(1) * pkin(1);
t185 = qJDD(2) * pkin(2);
t111 = qJD(1) ^ 2;
t182 = t104 * t111;
t181 = t105 * t106;
t180 = t105 * t107;
t179 = t106 * t108;
t178 = t107 * t108;
t176 = -qJD(4) - t19;
t169 = qJD(5) * t106;
t163 = t71 + 0.2e1 * t94 + 0.2e1 * t95;
t160 = t104 * t195;
t159 = t104 * t193;
t157 = t107 * t182;
t152 = t92 - t208;
t17 = t96 * qJD(2) + t158;
t119 = t124 * qJD(2);
t136 = pkin(2) * t164 + t218 * qJ(3) + qJD(1) * t79 + t100;
t123 = pkin(3) * t164 + qJDD(4) + t136;
t2 = qJD(1) * t119 + qJDD(1) * t135 + t123;
t147 = qJD(5) * t17 - t2;
t146 = t108 * pkin(1) + pkin(2) * t178 + t105 * pkin(6) + qJ(3) * t183;
t145 = t70 + t213;
t144 = qJD(1) * t34 + t19;
t140 = t176 * t104;
t139 = -qJD(2) * pkin(2) + qJD(3);
t138 = 0.2e1 * t149;
t137 = t104 * t154;
t110 = qJD(2) ^ 2;
t134 = pkin(6) * t110 + t208;
t132 = -qJDD(3) - t145;
t41 = t139 + t74;
t43 = t75 + t97;
t131 = t104 * t43 - t107 * t41;
t129 = t148 - t90;
t23 = t156 - t185;
t128 = -t57 * t169 - t196;
t127 = qJD(5) * t195 - t194;
t126 = -0.2e1 * pkin(1) * t166 - t190;
t121 = -t134 + 0.2e1 * t100;
t18 = -t137 + t202;
t5 = -qJD(1) * t137 + t123;
t118 = -qJD(1) * t18 - qJDD(1) * t34 + t208 - t5;
t117 = -qJ(4) * t165 - t132 + t56;
t116 = -t96 * t35 + (-t22 + t40) * t57;
t12 = pkin(2) * t150 - t136;
t27 = pkin(2) * t172 - t202;
t115 = -qJD(1) * t27 - qJDD(1) * t215 - t12 - t134;
t20 = -pkin(6) * t150 + t162;
t114 = -qJD(2) * t131 + t23 * t104 + t20 * t107;
t11 = qJD(5) * t37 - t106 * qJDD(2) + t217 * t103;
t102 = qJ(3) + pkin(4);
t88 = t108 * pkin(6);
t61 = g(1) * t180;
t60 = g(1) * t184;
t55 = qJ(3) * t178;
t53 = qJ(3) * t180;
t50 = -t111 * t98 - t110;
t46 = t203 * t107;
t44 = qJDD(2) + t157;
t39 = pkin(2) * t175 - t67;
t33 = -t103 * t105 + t104 * t179;
t32 = -t103 * t183 - t181;
t31 = -t103 * t108 - t104 * t181;
t30 = t103 * t184 - t179;
t25 = -qJ(4) * t172 + t125;
t24 = -t209 * t175 + t67;
t13 = t119 + t202;
t8 = -t151 + t112;
t4 = t103 * t14 + t106 * t17;
t3 = -t103 * t17 + t106 * t14;
t1 = t106 * t2;
t6 = [qJDD(1), t152, t214, qJDD(1) * t98 + t104 * t138, 0.2e1 * t104 * t164 - 0.2e1 * t199 * t166, qJDD(2) * t104 + t107 * t110, qJDD(2) * t107 - t104 * t110, 0, t104 * t126 + t107 * t121 + t61, -t104 * t121 + t107 * t126 - t60, t104 * t211 + t115 * t107 + t61, t198 * qJDD(1) * pkin(6) + t114 - t214, t115 * t104 - t107 * t211 + t60, t114 * pkin(6) - g(1) * t88 - g(2) * t146 + t29 * t27 + (t12 - t92) * t215, qJDD(2) * t46 + t60 + (t144 * t107 - t25) * qJD(2) - t118 * t104, qJDD(2) * t45 - t61 + (t144 * t104 + t26) * qJD(2) + t118 * t107, (-qJD(2) * t21 - qJDD(1) * t46 + t9 + (-qJD(2) * t45 + t25) * qJD(1)) * t107 + (-t188 - qJDD(1) * t45 - t8 + (qJD(2) * t46 - t26) * qJD(1)) * t104 + t214, t8 * t45 + t21 * t26 - t9 * t46 + t28 * t25 + t5 * t34 + t19 * t18 - g(1) * (-qJ(4) * t108 + t88) - g(2) * (pkin(3) * t178 + t146) + (-g(1) * (t215 - t86) + g(2) * qJ(4)) * t105, -t37 * t153 + (-t10 * t107 - t37 * t172) * t106, (t103 * t37 + t106 * t36) * t172 + (t197 - t106 * t11 + (t103 * t36 - t106 * t37) * qJD(5)) * t107, (t57 * t167 + t10) * t104 + (t127 - t186) * t107, (-t57 * t173 + t11) * t104 + (-t128 + t187) * t107, t104 * t35 + t57 * t171, t3 * t171 - g(1) * t31 - g(2) * t33 + t1 * t104 - t46 * t11 + t25 * t36 + (t13 * t57 + t206 + (-t17 * t104 - t22 * t107 - t45 * t57) * qJD(5)) * t106 + t212 * t103, -t4 * t171 - g(1) * t30 - g(2) * t32 + t46 * t10 + t25 * t37 + (-(-qJD(5) * t45 + t13) * t57 - t206 + t147 * t104 + t22 * t168) * t103 + t212 * t106; 0, 0, 0, -t157, t199 * t111, t165, t164, qJDD(2), pkin(1) * t182 - t145, t207 - t71 + (pkin(1) * t111 + t214) * t107, 0.2e1 * t185 + (-t104 * t29 + t107 * t39) * qJD(1) + t132, (-pkin(2) * t104 + t189) * qJDD(1) + ((t43 - t97) * t104 + (t139 - t41) * t107) * qJD(1), (qJD(1) * t39 - g(3)) * t104 + (qJD(1) * t29 - t214) * t107 + t163, t20 * qJ(3) + t43 * qJD(3) - t23 * pkin(2) - t29 * t39 - g(1) * (-pkin(2) * t183 + t55) - g(2) * (-pkin(2) * t184 + t53) - g(3) * t201 + t131 * qJD(1) * pkin(6), qJD(2) * t38 + t51 + (-g(3) + (-pkin(6) * qJD(2) - t24) * qJD(1)) * t104 + (-qJ(4) * qJDD(1) + t176 * qJD(1) - t214) * t107 + t163, -qJD(2) * t40 - 0.2e1 * t151 + ((-qJ(4) * qJD(2) + t24) * t107 + t140) * qJD(1) + t117, (-t189 + t219) * qJDD(1), -g(1) * t55 - g(2) * t53 - g(3) * t161 - t9 * qJ(3) - t158 * t28 - t19 * t24 - t209 * t8 - t21 * t40 + t214 * t219, t37 * t193 - t197, (-t10 - t205) * t106 + (-t11 - t204) * t103, (-t159 + t191) * qJD(1) + t128, (t160 - t192) * qJD(1) + t127, -t57 * t174, -t102 * t11 + t116 * t103 - t106 * t210 - t158 * t36 - t3 * t174, t102 * t10 + t103 * t210 + t116 * t106 - t158 * t37 + t4 * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t165, t50, -qJD(2) * t43 + t29 * t175 + t213 + t23, t50, t44, -t165, t188 - t151 + (-qJ(4) * t171 + t140) * qJD(1) + t117, 0, 0, 0, 0, 0, -t106 * t216 + t187 - t196, t103 * t216 + t186 - t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 + t165, 0.2e1 * t150 - t164, -t198 * t111, (-t107 * t28 + (t21 - t154) * t104) * qJD(1) + t123 + t152, 0, 0, 0, 0, 0, (-t160 - t192) * qJD(1) - t127, (-t159 - t191) * qJD(1) + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t36, -t36 ^ 2 + t37 ^ 2, t10 - t205, t11 - t204, t35, -g(1) * t32 + g(2) * t30 + t103 * t129 - t169 * t17 + t22 * t37 + t4 * t57 + t1, g(1) * t33 - g(2) * t31 + t103 * t147 + t106 * t129 - t22 * t36 + t3 * t57;];
tau_reg = t6;
