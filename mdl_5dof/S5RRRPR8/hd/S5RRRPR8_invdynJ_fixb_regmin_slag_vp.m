% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:29
% EndTime: 2019-12-31 21:20:36
% DurationCPUTime: 2.37s
% Computational Cost: add. (2817->335), mult. (6277->426), div. (0->0), fcn. (4393->10), ass. (0->188)
t143 = sin(qJ(3));
t144 = sin(qJ(2));
t147 = cos(qJ(2));
t240 = cos(qJ(3));
t96 = t143 * t147 + t240 * t144;
t87 = t96 * qJD(1);
t251 = qJD(5) + t87;
t137 = qJD(2) + qJD(3);
t142 = sin(qJ(5));
t146 = cos(qJ(5));
t202 = t240 * t147;
t188 = qJD(1) * t202;
t213 = qJD(1) * t144;
t201 = t143 * t213;
t85 = -t188 + t201;
t65 = t137 * t142 - t146 * t85;
t252 = t251 * t65;
t243 = pkin(7) + pkin(6);
t200 = qJD(3) * t240;
t103 = t243 * t147;
t99 = qJD(1) * t103;
t89 = t143 * t99;
t102 = t243 * t144;
t97 = qJD(1) * t102;
t58 = -t240 * t97 - t89;
t228 = pkin(2) * t200 + qJD(4) - t58;
t212 = qJD(3) * t143;
t92 = t240 * t99;
t57 = -t143 * t97 + t92;
t186 = pkin(2) * t212 - t57;
t233 = qJD(2) * pkin(2);
t93 = -t97 + t233;
t53 = -t240 * t93 + t89;
t217 = qJD(4) + t53;
t250 = -qJD(5) + t251;
t141 = qJ(2) + qJ(3);
t134 = sin(t141);
t135 = cos(t141);
t215 = t135 * pkin(3) + t134 * qJ(4);
t136 = qJDD(2) + qJDD(3);
t130 = t136 * qJ(4);
t249 = -t137 * qJD(4) - t130;
t128 = -t240 * pkin(2) - pkin(3);
t121 = -pkin(8) + t128;
t242 = t85 * pkin(4);
t221 = t143 * t144;
t178 = t137 * t221;
t195 = qJDD(1) * t240;
t207 = t147 * qJDD(1);
t189 = -t137 * t188 - t143 * t207 - t144 * t195;
t39 = qJD(1) * t178 + t189;
t35 = -qJDD(5) + t39;
t248 = (t242 + t186) * t251 - t121 * t35;
t145 = sin(qJ(1));
t148 = cos(qJ(1));
t183 = g(1) * t145 - g(2) * t148;
t203 = qJD(2) * t243;
t100 = t147 * t203;
t98 = t144 * t203;
t28 = t143 * t100 + t102 * t200 + t103 * t212 + t240 * t98;
t69 = -t143 * t102 + t240 * t103;
t247 = -t134 * t183 - t136 * t69 + t137 * t28;
t29 = t69 * qJD(3) + t240 * t100 - t143 * t98;
t68 = t240 * t102 + t143 * t103;
t246 = t135 * t183 - t136 * t68 - t137 * t29;
t245 = t87 ^ 2;
t244 = pkin(3) + pkin(8);
t241 = t87 * pkin(4);
t239 = pkin(2) * t147;
t238 = pkin(3) * t136;
t125 = g(3) * t134;
t126 = g(3) * t135;
t54 = t143 * t93 + t92;
t50 = -qJ(4) * t137 - t54;
t31 = -t50 - t242;
t95 = -t202 + t221;
t237 = t31 * t95;
t129 = pkin(1) + t239;
t174 = -qJ(4) * t96 - t129;
t41 = t244 * t95 + t174;
t236 = t41 * t35;
t235 = t251 * t85;
t234 = t87 * t85;
t231 = t137 * t54;
t210 = qJD(5) * t146;
t211 = qJD(5) * t142;
t208 = t144 * qJDD(1);
t177 = t143 * t208 - t147 * t195;
t60 = t137 * t96;
t40 = qJD(1) * t60 + t177;
t19 = t146 * t136 - t137 * t211 + t142 * t40 + t85 * t210;
t230 = t19 * t146;
t229 = t241 + t228;
t227 = t134 * t145;
t226 = t134 * t148;
t225 = t135 * t145;
t224 = t135 * t148;
t223 = t142 * t145;
t222 = t142 * t148;
t220 = t145 * t146;
t219 = t146 * t148;
t218 = t241 + t217;
t139 = t144 ^ 2;
t214 = -t147 ^ 2 + t139;
t209 = qJD(1) * qJD(2);
t132 = t144 * t233;
t205 = t95 * t210;
t25 = -t244 * t137 + t218;
t101 = t129 * qJD(1);
t163 = -qJ(4) * t87 - t101;
t27 = t244 * t85 + t163;
t13 = t142 * t25 + t146 * t27;
t198 = t147 * t209;
t63 = qJDD(2) * pkin(2) + t243 * (-t198 - t208);
t199 = t144 * t209;
t64 = t243 * (-t199 + t207);
t192 = -t143 * t63 - t93 * t200 + t99 * t212 - t240 * t64;
t15 = t192 + t249;
t8 = -pkin(4) * t40 - t15;
t204 = -t13 * t85 + t8 * t146;
t82 = pkin(2) * t199 - qJDD(1) * t129;
t155 = qJ(4) * t39 - qJD(4) * t87 + t82;
t2 = t244 * t40 + t155;
t197 = qJD(5) * t25 + t2;
t193 = t143 * t64 + t99 * t200 + t93 * t212 - t240 * t63;
t181 = qJDD(4) + t193;
t7 = -pkin(4) * t39 - t244 * t136 + t181;
t196 = -qJD(5) * t27 + t7;
t194 = t251 * t31;
t191 = t136 * t142 - t146 * t40;
t190 = t142 * t251;
t185 = -pkin(2) * t144 - pkin(3) * t134;
t51 = pkin(3) * t87 + qJ(4) * t85;
t184 = g(1) * t148 + g(2) * t145;
t182 = t251 * t60 - t35 * t95;
t180 = t244 * t35 - (t54 - t242) * t251;
t12 = -t142 * t27 + t146 * t25;
t179 = t12 * t85 + t8 * t142 + (t146 * t87 + t210) * t31;
t67 = t137 * t146 + t142 * t85;
t47 = pkin(2) * t213 + t51;
t173 = -t146 * t35 - t190 * t251;
t172 = t129 + t215;
t59 = -qJD(2) * t202 - t147 * t200 + t178;
t171 = qJ(4) * t59 - qJD(4) * t96 + t132;
t170 = -0.2e1 * pkin(1) * t209 - pkin(6) * qJDD(2);
t169 = -g(1) * t226 - g(2) * t227 + t126 + t193;
t168 = -g(1) * t224 - g(2) * t225 - t125 - t192;
t48 = t96 * pkin(4) + t68;
t165 = t31 * t60 + t48 * t35 + t8 * t95;
t164 = -t146 * t251 ^ 2 + t142 * t35;
t162 = t101 * t87 - t169;
t161 = -t101 * t85 - t168;
t160 = -t135 * t184 - t125;
t151 = qJD(2) ^ 2;
t159 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t151 + t183;
t152 = qJD(1) ^ 2;
t158 = pkin(1) * t152 - pkin(6) * qJDD(1) + t184;
t45 = pkin(3) * t85 + t163;
t157 = t45 * t87 + qJDD(4) + t169;
t156 = -t45 * t85 + t168 - t249;
t23 = -t189 + (-t201 + t85) * t137;
t80 = t87 * pkin(8);
t154 = (-qJD(5) * t121 + t47 + t80) * t251 + t160;
t153 = (qJD(5) * t244 + t51 + t80) * t251 + t160;
t123 = pkin(2) * t143 + qJ(4);
t105 = qJ(4) * t224;
t104 = qJ(4) * t225;
t79 = -t134 * t223 + t219;
t78 = t134 * t220 + t222;
t77 = t134 * t222 + t220;
t76 = t134 * t219 - t223;
t52 = pkin(3) * t95 + t174;
t49 = -t95 * pkin(4) + t69;
t46 = -pkin(3) * t137 + t217;
t42 = -t85 ^ 2 + t245;
t21 = pkin(3) * t60 + t171;
t20 = qJD(5) * t67 + t191;
t18 = -t59 * pkin(4) + t29;
t17 = -pkin(4) * t60 - t28;
t16 = t181 - t238;
t14 = t244 * t60 + t171;
t11 = pkin(3) * t40 + t155;
t10 = -t65 * t85 + t164;
t9 = t67 * t85 + t173;
t4 = t146 * t7;
t3 = -t190 * t67 + t230;
t1 = (-t251 * t67 - t20) * t146 + (-t19 + t252) * t142;
t5 = [qJDD(1), t183, t184, qJDD(1) * t139 + 0.2e1 * t144 * t198, 0.2e1 * t144 * t207 - 0.2e1 * t214 * t209, qJDD(2) * t144 + t147 * t151, qJDD(2) * t147 - t144 * t151, 0, t144 * t170 + t147 * t159, -t144 * t159 + t147 * t170, -t39 * t96 - t59 * t87, t39 * t95 - t40 * t96 + t59 * t85 - t60 * t87, t136 * t96 - t137 * t59, -t136 * t95 - t137 * t60, 0, -t101 * t60 - t129 * t40 + t132 * t85 + t82 * t95 + t246, t101 * t59 + t129 * t39 + t87 * t132 + t82 * t96 + t247, t15 * t95 + t16 * t96 + t28 * t85 + t29 * t87 - t39 * t68 - t40 * t69 - t46 * t59 + t50 * t60 - t184, -t11 * t95 - t21 * t85 - t40 * t52 - t45 * t60 - t246, -t11 * t96 - t21 * t87 + t39 * t52 + t45 * t59 - t247, t11 * t52 - t15 * t69 + t16 * t68 + t45 * t21 + t50 * t28 + t46 * t29 + (-g(1) * t243 - g(2) * t172) * t148 + (g(1) * t172 - g(2) * t243) * t145, t67 * t205 + (t19 * t95 + t60 * t67) * t142, (-t142 * t65 + t146 * t67) * t60 + (-t142 * t20 + t230 + (-t142 * t67 - t146 * t65) * qJD(5)) * t95, t142 * t182 + t19 * t96 + t205 * t251 - t59 * t67, -t211 * t251 * t95 + t146 * t182 - t20 * t96 + t59 * t65, -t251 * t59 - t35 * t96, -g(1) * t79 - g(2) * t77 - t12 * t59 + t17 * t65 + t49 * t20 + t4 * t96 + (-t14 * t251 - t2 * t96 + t236) * t142 + (t18 * t251 - t165) * t146 + ((-t142 * t48 - t146 * t41) * t251 - t13 * t96 + t142 * t237) * qJD(5), g(1) * t78 - g(2) * t76 + t13 * t59 + t17 * t67 + t49 * t19 + (-(qJD(5) * t48 + t14) * t251 + t236 - t197 * t96 + qJD(5) * t237) * t146 + (-(-qJD(5) * t41 + t18) * t251 - t196 * t96 + t165) * t142; 0, 0, 0, -t144 * t152 * t147, t214 * t152, t208, t207, qJDD(2), -g(3) * t147 + t144 * t158, g(3) * t144 + t147 * t158, t234, t42, t23, -t177, t136, t57 * t137 + (t240 * t136 - t137 * t212 - t85 * t213) * pkin(2) + t162, t58 * t137 + (-t136 * t143 - t137 * t200 - t87 * t213) * pkin(2) + t161, -t123 * t40 - t128 * t39 + (t186 - t50) * t87 + (t46 - t228) * t85, t47 * t85 + t186 * t137 + (-pkin(3) + t128) * t136 + t157, t123 * t136 + t228 * t137 + t47 * t87 + t156, -t15 * t123 + t16 * t128 - t45 * t47 - g(1) * (t148 * t185 + t105) - g(2) * (t145 * t185 + t104) - g(3) * (t215 + t239) - t228 * t50 + t186 * t46, t3, t1, t9, t10, t235, t123 * t20 + t154 * t142 + t146 * t248 + t229 * t65 + t179, t123 * t19 + t229 * t67 + t154 * t146 + (-t194 - t248) * t142 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t42, t23, -t177, t136, t162 + t231, -t137 * t53 + t161, pkin(3) * t39 - qJ(4) * t40 + (-t50 - t54) * t87 + (t46 - t217) * t85, t51 * t85 + t157 - t231 - 0.2e1 * t238, t217 * t137 + t51 * t87 + t130 + t156, -t15 * qJ(4) - t16 * pkin(3) - t45 * t51 - t46 * t54 - g(1) * (-pkin(3) * t226 + t105) - g(2) * (-pkin(3) * t227 + t104) - g(3) * t215 - t217 * t50, t3, t1, t9, t10, t235, qJ(4) * t20 + t142 * t153 + t146 * t180 + t218 * t65 + t179, qJ(4) * t19 + t218 * t67 + (-t194 - t180) * t142 + t153 * t146 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t136 - t234, -t137 ^ 2 - t245, t137 * t50 + t157 - t238, 0, 0, 0, 0, 0, -t137 * t65 + t173, -t137 * t67 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t65, -t65 ^ 2 + t67 ^ 2, t19 + t252, t250 * t67 - t191, -t35, -g(1) * t76 - g(2) * t78 + t126 * t146 + t13 * t250 - t142 * t2 - t31 * t67 + t4, g(1) * t77 - g(2) * t79 + t12 * t251 + t31 * t65 - t197 * t146 + (-t196 - t126) * t142;];
tau_reg = t5;
