% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:32
% DurationCPUTime: 2.38s
% Computational Cost: add. (2445->294), mult. (3749->415), div. (0->0), fcn. (2513->14), ass. (0->196)
t123 = qJDD(1) + qJDD(2);
t126 = qJD(1) + qJD(2);
t131 = sin(pkin(9));
t137 = cos(qJ(5));
t138 = cos(qJ(4));
t195 = qJD(4) + qJD(5);
t167 = t195 * t138;
t133 = sin(qJ(5));
t211 = t133 * t138;
t134 = sin(qJ(4));
t225 = t123 * t134;
t199 = qJD(4) * t134;
t179 = t131 * t199;
t212 = t133 * t134;
t184 = t131 * t212;
t249 = -qJD(5) * t184 - t133 * t179;
t15 = t249 * t126 + (t123 * t211 + (t126 * t167 + t225) * t137) * t131;
t124 = t131 ^ 2;
t132 = cos(pkin(9));
t220 = t126 * t132;
t96 = -qJD(4) + t220;
t209 = qJD(4) + t96;
t135 = sin(qJ(2));
t236 = pkin(1) * qJD(1);
t190 = t135 * t236;
t88 = qJ(3) * t126 + t190;
t250 = (t124 * t126 + t209 * t132) * t88;
t139 = cos(qJ(2));
t203 = qJD(1) * t139;
t161 = -pkin(1) * t203 + qJD(3);
t214 = t132 * t138;
t91 = -pkin(3) * t132 - pkin(7) * t131 - pkin(2);
t230 = qJ(3) * t214 + t134 * t91;
t198 = qJD(4) * t138;
t213 = t132 * t139;
t247 = qJD(3) * t214 + t91 * t198 - (t134 * t135 + t138 * t213) * t236;
t221 = t126 * t131;
t130 = qJ(1) + qJ(2);
t119 = sin(t130);
t114 = g(1) * t119;
t202 = qJD(2) * t135;
t189 = pkin(1) * t202;
t243 = pkin(1) * t139;
t206 = -qJD(1) * t189 + qJDD(1) * t243;
t245 = t114 + t206;
t242 = pkin(2) * t123;
t121 = cos(t130);
t241 = g(2) * t121;
t240 = g(3) * t131;
t217 = t131 * t138;
t185 = t126 * t217;
t51 = t126 * t184 - t137 * t185;
t156 = t134 * t137 + t211;
t52 = t156 * t221;
t239 = t51 * t52;
t111 = pkin(1) * t135 + qJ(3);
t76 = t91 - t243;
t237 = t111 * t214 + t134 * t76;
t196 = qJDD(1) * t135;
t201 = qJD(2) * t139;
t54 = t123 * qJ(3) + t126 * qJD(3) + (qJD(1) * t201 + t196) * pkin(1);
t48 = t124 * t54;
t90 = -qJD(5) + t96;
t235 = t132 * t90;
t188 = t88 * t214;
t47 = t91 * t126 + t161;
t153 = -t134 * t47 - t188;
t218 = t131 * t134;
t193 = pkin(8) * t218;
t26 = -t126 * t193 - t153;
t234 = t137 * t26;
t216 = t132 * t123;
t95 = -qJDD(4) + t216;
t89 = -qJDD(5) + t95;
t233 = t89 * t132;
t232 = t95 * t132;
t229 = qJ(3) * t134;
t228 = t111 * t134;
t227 = t119 * t132;
t226 = t121 * t132;
t224 = t123 * t138;
t122 = t126 ^ 2;
t223 = t124 * t122;
t219 = t126 * t134;
t215 = t132 * t134;
t210 = t134 * t138;
t208 = t121 * pkin(2) + t119 * qJ(3);
t207 = g(1) * t121 + g(2) * t119;
t125 = t132 ^ 2;
t205 = t124 + t125;
t128 = t138 ^ 2;
t204 = t134 ^ 2 - t128;
t200 = qJD(3) * t134;
t197 = qJD(5) * t133;
t176 = qJDD(3) - t206;
t39 = t91 * t123 + t176;
t194 = t134 * t39 + t47 * t198 + t54 * t214;
t192 = pkin(8) * t217;
t104 = pkin(1) * t201 + qJD(3);
t191 = t104 * t214 + t134 * t189 + t76 * t198;
t183 = t123 * t217;
t182 = qJ(3) * t198;
t181 = t126 * t202;
t180 = t126 * t198;
t178 = t132 * t199;
t61 = t176 - t242;
t177 = -t61 - t241;
t175 = -pkin(2) * t119 + t121 * qJ(3);
t42 = t138 * t47;
t25 = -pkin(8) * t185 - t88 * t215 + t42;
t17 = -pkin(4) * t96 + t25;
t147 = -t88 * t178 + t194;
t150 = t180 + t225;
t9 = -t150 * t131 * pkin(8) + t147;
t174 = qJD(5) * t17 + t9;
t173 = t205 * t54;
t172 = t95 + t216;
t171 = t104 * t205;
t170 = t126 * t209;
t169 = t205 * t123;
t168 = t125 * t54 - t207 + t48;
t36 = t138 * t39;
t164 = -t54 * t215 + t36;
t80 = t138 * t91;
t40 = -t192 + t80 + (-pkin(4) - t229) * t132;
t163 = qJD(5) * t40 + (-qJ(3) * t215 - t192) * qJD(4) + t247;
t100 = pkin(8) * t179;
t46 = -t193 + t230;
t62 = (-t134 * t213 + t135 * t138) * t236;
t162 = qJD(4) * t230 + qJD(5) * t46 + t132 * t200 - t100 + t62;
t160 = -t133 * t17 - t234;
t73 = t138 * t76;
t30 = -t192 + t73 + (-pkin(4) - t228) * t132;
t37 = -t193 + t237;
t159 = -t133 * t37 + t137 * t30;
t158 = t133 * t30 + t137 * t37;
t157 = qJD(4) * (t96 + t220);
t155 = -t137 * t138 + t212;
t101 = t131 * pkin(4) * t198;
t154 = -t161 * t131 - t101;
t65 = t119 * t215 + t121 * t138;
t67 = t119 * t138 - t121 * t215;
t152 = -g(1) * t65 - g(2) * t67 + t147 * t132 + t138 * t48;
t66 = -t119 * t214 + t121 * t134;
t68 = t119 * t134 + t121 * t214;
t151 = t124 * t88 * t198 - g(1) * t66 - g(2) * t68 + t134 * t48;
t149 = t126 * t190 - t241;
t22 = t26 * t197;
t29 = (t150 * pkin(4) + t54) * t131;
t142 = t195 * t156;
t31 = t142 * t131;
t50 = (pkin(4) * t219 + t88) * t131;
t129 = qJ(4) + qJ(5);
t118 = sin(t129);
t120 = cos(t129);
t56 = t118 * t227 + t120 * t121;
t58 = -t118 * t226 + t119 * t120;
t71 = t155 * t131;
t8 = -pkin(8) * t183 - pkin(4) * t95 + (-t188 + (pkin(8) * t221 - t47) * t134) * qJD(4) + t164;
t148 = -g(1) * t56 - g(2) * t58 + (t133 * t8 + t174 * t137 - t22) * t132 - t29 * t71 - t50 * t31;
t146 = -t96 ^ 2 - t223;
t2 = t160 * qJD(5) - t133 * t9 + t137 * t8;
t32 = t137 * t131 * t167 + t249;
t57 = t118 * t121 - t120 * t227;
t59 = t118 * t119 + t120 * t226;
t70 = t156 * t131;
t145 = -g(1) * t57 - g(2) * t59 - t2 * t132 + t29 * t70 + t50 * t32;
t144 = t161 * t205;
t143 = t22 + t120 * t240 + g(1) * t59 + (t26 * t90 - t8) * t133 - g(2) * t57 + t50 * t52;
t141 = -g(1) * t58 + g(2) * t56 + t118 * t240 + t50 * t51 + t2;
t14 = t137 * t183 + (-t123 * t212 - t142 * t126) * t131;
t140 = cos(qJ(1));
t136 = sin(qJ(1));
t116 = -pkin(2) - t243;
t107 = pkin(4) * t218;
t103 = t138 * t189;
t99 = g(1) * t227;
t84 = qJ(3) * t131 + t107;
t83 = -pkin(2) * t126 + t161;
t74 = t111 * t131 + t107;
t60 = t104 * t131 + t101;
t55 = (t123 * t128 - 0.2e1 * t134 * t180) * t124;
t38 = 0.2e1 * (t204 * t126 * qJD(4) - t123 * t210) * t124;
t28 = (t172 * t134 + t138 * t157) * t131;
t27 = (t134 * t157 - t172 * t138) * t131;
t24 = -t237 * qJD(4) - t104 * t215 + t100 + t103;
t23 = (-t111 * t215 - t192) * qJD(4) + t191;
t16 = t51 ^ 2 - t52 ^ 2;
t13 = t51 * t90 - t15;
t12 = -t52 * t90 + t14;
t11 = t153 * qJD(4) + t164;
t6 = -t14 * t71 + t31 * t51;
t5 = t132 * t15 + t32 * t90 + t70 * t89;
t4 = -t132 * t14 + t31 * t90 + t71 * t89;
t3 = -t14 * t70 + t15 * t71 + t31 * t52 + t32 * t51;
t1 = [qJDD(1), g(1) * t136 - g(2) * t140, g(1) * t140 + g(2) * t136, t123, -t241 + (t123 * t139 - t181) * pkin(1) + t245, ((-qJDD(1) - t123) * t135 + (-qJD(1) - t126) * t201) * pkin(1) + t207, t99 + (-pkin(1) * t181 - t116 * t123 + t177) * t132, t111 * t169 + t126 * t171 + t168, t61 * t116 + t83 * t189 - g(1) * (-pkin(1) * t136 + t175) - g(2) * (pkin(1) * t140 + t208) + t111 * t173 + t88 * t171, t55, t38, t27, t28, t232, -(-t76 * t199 + t103) * t96 - t73 * t95 + (-(-t104 * t134 - t111 * t198) * t96 + t95 * t228 - t11) * t132 + (t104 * t219 + t150 * t111) * t124 + t151, (-t111 * t178 + t191) * t96 + t237 * t95 + ((t104 * t126 + t111 * t123) * t138 + (-t111 * t126 - t88) * t199) * t124 + t152, t6, t3, t4, t5, t233, -(-qJD(5) * t158 - t133 * t23 + t137 * t24) * t90 - t159 * t89 + t60 * t52 + t74 * t15 + t145, (qJD(5) * t159 + t133 * t24 + t137 * t23) * t90 + t158 * t89 - t60 * t51 + t74 * t14 + t148; 0, 0, 0, t123, t149 + t245, (-t196 + (-qJD(2) + t126) * t203) * pkin(1) + t207, t99 + (t149 - t61 + t242) * t132, qJ(3) * t169 + t144 * t126 + t168, -t61 * pkin(2) - g(1) * t175 - g(2) * t208 + qJ(3) * t173 + t144 * t88 - t83 * t190, t55, t38, t27, t28, t232, -t80 * t95 + (t91 * t199 + t62) * t96 + (-(-t182 - t200) * t96 + t95 * t229 - t11) * t132 + (qJ(3) * t225 + (t161 * t134 + t182) * t126) * t124 + t151, t230 * t95 + (-qJ(3) * t178 + t247) * t96 + (qJ(3) * t224 - t88 * t199 + (-qJ(3) * t199 + t161 * t138) * t126) * t124 + t152, t6, t3, t4, t5, t233, -(-t133 * t46 + t137 * t40) * t89 + t84 * t15 + (t133 * t163 + t137 * t162) * t90 - t154 * t52 + t145, (t133 * t40 + t137 * t46) * t89 + t84 * t14 + (-t133 * t162 + t137 * t163) * t90 + t154 * t51 + t148; 0, 0, 0, 0, 0, 0, -t216, -t205 * t122, -t205 * t88 * t126 - t114 - t177, 0, 0, 0, 0, 0, t146 * t134 - t138 * t95, t134 * t95 + t146 * t138, 0, 0, 0, 0, 0, t142 * t90 + t155 * t89 + (-t131 * t52 - t156 * t235) * t126, t51 * t221 + t156 * t89 + (t235 * t126 - t195 * t90) * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, t210 * t223, -t204 * t223, (-t134 * t170 + t224) * t131, (-t138 * t170 - t225) * t131, -t95, -g(1) * t67 + g(2) * t65 + t36 - t138 * t250 + (-t132 * t54 - t209 * t47 + t240) * t134, g(1) * t68 - g(2) * t66 + g(3) * t217 + t134 * t250 - t42 * t96 - t194, -t239, t16, t12, t13, -t89, (-t133 * t25 - t234) * t90 + (-t137 * t89 - t185 * t52 + t197 * t90) * pkin(4) + t141, (-t25 * t90 - t174) * t137 + (qJD(5) * t137 * t90 + t133 * t89 + t185 * t51) * pkin(4) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, t16, t12, t13, -t89, t160 * t90 + t141, (-t9 + (-qJD(5) - t90) * t17) * t137 + t143;];
tau_reg = t1;
