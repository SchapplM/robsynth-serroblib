% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:22
% EndTime: 2019-12-31 18:04:26
% DurationCPUTime: 2.22s
% Computational Cost: add. (2578->318), mult. (6091->417), div. (0->0), fcn. (4584->10), ass. (0->171)
t157 = sin(qJ(5));
t160 = cos(qJ(5));
t155 = cos(pkin(8));
t161 = cos(qJ(4));
t154 = sin(pkin(8));
t158 = sin(qJ(4));
t214 = t154 * t158;
t95 = t155 * t161 + t214;
t85 = t95 * qJD(1);
t207 = qJD(1) * t155;
t191 = t158 * t207;
t208 = qJD(1) * t154;
t193 = t161 * t208;
t87 = -t191 + t193;
t235 = t157 * t87 + t160 * t85;
t243 = t235 ^ 2;
t151 = qJD(4) + qJD(5);
t242 = t151 * t235;
t179 = -t157 * t85 + t160 * t87;
t223 = t235 * t179;
t136 = t154 * qJ(3);
t233 = -t155 * pkin(2) - pkin(1) - t136;
t241 = -t155 * pkin(3) + t233;
t237 = t179 ^ 2;
t240 = t237 - t243;
t200 = qJD(5) * t160;
t201 = qJD(5) * t157;
t134 = t154 * qJDD(1);
t197 = t155 * qJDD(1);
t180 = -t161 * t134 + t158 * t197;
t83 = t95 * qJD(4);
t42 = qJD(1) * t83 + t180;
t172 = qJD(4) * t191 - qJDD(1) * t95;
t202 = qJD(4) * t161;
t192 = t154 * t202;
t43 = qJD(1) * t192 - t172;
t8 = t157 * t43 + t160 * t42 + t85 * t200 + t201 * t87;
t239 = -t8 + t242;
t221 = -pkin(6) + qJ(2);
t104 = t221 * t155;
t97 = qJD(1) * t104;
t218 = t158 * t97;
t106 = qJ(2) * t208 + qJD(3);
t91 = -pkin(6) * t208 + t106;
t44 = t161 * t91 - t218;
t25 = -pkin(7) * t87 + t44;
t24 = qJD(4) * pkin(4) + t25;
t45 = t158 * t91 + t161 * t97;
t26 = -t85 * pkin(7) + t45;
t199 = qJD(1) * qJD(2);
t89 = qJ(2) * t134 + t154 * t199 + qJDD(3);
t64 = -pkin(6) * t134 + t89;
t66 = (t221 * qJDD(1) + t199) * t155;
t184 = -t158 * t66 + t161 * t64;
t17 = -t45 * qJD(4) + t184;
t6 = qJDD(4) * pkin(4) + t42 * pkin(7) + t17;
t196 = -t158 * t64 - t161 * t66 - t91 * t202;
t203 = qJD(4) * t158;
t16 = -t203 * t97 - t196;
t7 = -pkin(7) * t43 + t16;
t1 = (qJD(5) * t24 + t7) * t160 + t157 * t6 - t26 * t201;
t79 = -qJD(1) * pkin(1) - pkin(2) * t207 - qJ(3) * t208 + qJD(2);
t58 = pkin(3) * t207 - t79;
t31 = pkin(4) * t85 + t58;
t159 = sin(qJ(1));
t152 = qJ(4) + qJ(5);
t137 = sin(t152);
t138 = cos(t152);
t178 = t137 * t154 + t138 * t155;
t61 = t178 * t159;
t162 = cos(qJ(1));
t63 = t178 * t162;
t215 = t138 * t154;
t81 = t137 * t155 - t215;
t238 = g(1) * t63 + g(2) * t61 - g(3) * t81 + t31 * t235 - t1;
t236 = t151 * t179;
t145 = g(2) * t159;
t209 = g(1) * t162 + t145;
t96 = t154 * t161 - t155 * t158;
t75 = t96 * t162;
t234 = -g(1) * t75 + g(3) * t95;
t198 = qJD(1) * qJD(3);
t121 = t154 * t198;
t232 = -pkin(2) * t197 - qJ(3) * t134 - t121;
t217 = t160 * t26;
t11 = t157 * t24 + t217;
t2 = -qJD(5) * t11 - t157 * t7 + t160 * t6;
t213 = t155 * t159;
t60 = t137 * t213 - t159 * t215;
t62 = t81 * t162;
t231 = g(1) * t62 + g(2) * t60 + g(3) * t178 - t31 * t179 + t2;
t9 = qJD(5) * t179 - t157 * t42 + t160 * t43;
t230 = -t9 + t236;
t229 = t85 ^ 2;
t228 = t87 ^ 2;
t227 = pkin(4) * t43;
t224 = g(1) * t159;
t146 = g(2) * t162;
t222 = t87 * t85;
t103 = t221 * t154;
t53 = t158 * t103 + t161 * t104;
t220 = t157 * t26;
t153 = qJDD(1) * pkin(1);
t212 = t155 * t162;
t150 = t155 ^ 2;
t194 = 0.2e1 * t199;
t116 = t150 * t194;
t135 = t150 * qJDD(1);
t166 = qJ(2) ^ 2;
t211 = qJ(2) * t116 + t166 * t135;
t210 = t162 * pkin(1) + t159 * qJ(2);
t206 = qJD(2) * t158;
t205 = qJD(2) * t161;
t204 = qJD(3) * t154;
t132 = qJDD(2) - t153;
t195 = pkin(4) * t214;
t189 = t154 * t197;
t142 = t162 * qJ(2);
t188 = -t159 * pkin(1) + t142;
t187 = -t146 + t224;
t52 = t161 * t103 - t104 * t158;
t149 = t154 ^ 2;
t165 = qJD(1) ^ 2;
t101 = (-t149 - t150) * t165;
t183 = 0.2e1 * qJ(2) * t135 + t116 - t209;
t182 = pkin(2) * t212 + t162 * t136 + t210;
t29 = -pkin(7) * t96 + t52;
t30 = -pkin(7) * t95 + t53;
t14 = -t157 * t30 + t160 * t29;
t15 = t157 * t29 + t160 * t30;
t47 = -t157 * t95 + t160 * t96;
t99 = t157 * t161 + t160 * t158;
t98 = -t157 * t158 + t160 * t161;
t56 = t132 + t232;
t176 = -t187 + t132;
t175 = -t132 + t153 - t146;
t174 = -qJDD(1) * t233 - t146 - t56;
t127 = pkin(3) * t197;
t48 = t127 - t56;
t27 = t103 * t202 - t104 * t203 + t154 * t206 + t155 * t205;
t173 = (qJ(2) * qJDD(1) + t199) * t149;
t171 = t176 + t232;
t170 = -t127 + t171;
t28 = -qJD(4) * t53 + t154 * t205 - t155 * t206;
t164 = qJD(4) ^ 2;
t163 = -pkin(7) - pkin(6);
t148 = qJDD(4) + qJDD(5);
t143 = g(3) * t155;
t133 = t149 * qJDD(1);
t131 = pkin(4) * t161 + pkin(3);
t118 = g(1) * t213;
t84 = -t155 * t203 + t192;
t76 = t95 * t162;
t74 = t95 * t159;
t73 = t96 * t159;
t57 = pkin(4) * t84 + t204;
t51 = t151 * t99;
t50 = t151 * t98;
t49 = pkin(4) * t95 - t241;
t46 = t157 * t96 + t160 * t95;
t22 = t83 * pkin(7) + t28;
t21 = -pkin(7) * t84 + t27;
t20 = t48 + t227;
t19 = qJD(5) * t47 - t157 * t83 + t160 * t84;
t18 = t157 * t84 + t160 * t83 + t200 * t95 + t201 * t96;
t13 = t160 * t25 - t220;
t12 = -t157 * t25 - t217;
t10 = t160 * t24 - t220;
t4 = -qJD(5) * t15 - t157 * t21 + t160 * t22;
t3 = qJD(5) * t14 + t157 * t22 + t160 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), t187, t209, 0, 0, t133, 0.2e1 * t189, 0, t135, 0, 0, t155 * t175 + t118, (-t175 - t224) * t154, 0.2e1 * t173 + t183, -t132 * pkin(1) - g(1) * t188 - g(2) * t210 + (qJ(2) * t194 + qJDD(1) * t166) * t149 + t211, t133, 0, -0.2e1 * t189, 0, 0, t135, t118 + (t174 + t121) * t155, t154 * t89 + t173 + t183, t149 * t198 + (t174 + t224) * t154, t56 * t233 - g(1) * (-pkin(2) * t213 + t188) - g(2) * t182 + (t89 * qJ(2) + qJ(3) * t224 + t106 * qJD(2) - t79 * qJD(3)) * t154 + t211, -t42 * t96 - t83 * t87, t42 * t95 - t43 * t96 + t83 * t85 - t84 * t87, -qJD(4) * t83 + qJDD(4) * t96, t43 * t95 + t84 * t85, -qJD(4) * t84 - qJDD(4) * t95, 0, g(1) * t74 - g(2) * t76 + qJD(4) * t28 + qJDD(4) * t52 + t204 * t85 - t241 * t43 + t48 * t95 + t58 * t84, g(1) * t73 - g(2) * t75 - qJD(4) * t27 - qJDD(4) * t53 + t204 * t87 + t241 * t42 + t48 * t96 - t58 * t83, -t16 * t95 - t17 * t96 - t27 * t85 - t28 * t87 + t42 * t52 - t43 * t53 + t44 * t83 - t45 * t84 + t209, t16 * t53 + t45 * t27 + t17 * t52 + t44 * t28 - t48 * t241 + t58 * t204 - g(1) * (-t162 * pkin(6) + t142) - g(2) * (pkin(3) * t212 + t182) + (g(2) * pkin(6) - g(1) * t241) * t159, -t179 * t18 - t47 * t8, -t179 * t19 + t18 * t235 + t46 * t8 - t47 * t9, t148 * t47 - t151 * t18, t19 * t235 + t46 * t9, -t148 * t46 - t151 * t19, 0, g(1) * t61 - g(2) * t63 + t14 * t148 + t151 * t4 + t19 * t31 + t20 * t46 + t235 * t57 + t49 * t9, -g(1) * t60 + g(2) * t62 - t148 * t15 - t151 * t3 + t179 * t57 - t18 * t31 + t20 * t47 - t49 * t8, -t1 * t46 + t10 * t18 - t11 * t19 + t14 * t8 - t15 * t9 - t179 * t4 - t2 * t47 - t235 * t3 + t209, t1 * t15 + t11 * t3 + t2 * t14 + t10 * t4 + t20 * t49 + t31 * t57 - g(1) * (t162 * t163 + t142) - g(2) * (t131 * t212 + t162 * t195 + t182) + (-g(1) * (-t131 * t155 - t195 + t233) - g(2) * t163) * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t134, t101, qJ(2) * t101 + t176, 0, 0, 0, 0, 0, 0, -t197, t101, -t134, -qJ(2) * t165 * t150 - t106 * t208 + t171, 0, 0, 0, 0, 0, 0, (-t87 - t193) * qJD(4) + t172, 0.2e1 * t85 * qJD(4) + t180, t228 + t229, -t44 * t87 - t45 * t85 + t170, 0, 0, 0, 0, 0, 0, -t9 - t236, t8 + t242, t237 + t243, -t10 * t179 - t11 * t235 + t170 - t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154 * t165 * t155, t134, -t149 * t165, t143 + (qJD(1) * t79 - t209) * t154 + t89, 0, 0, 0, 0, 0, 0, qJDD(4) * t161 - t164 * t158 - t208 * t85, -qJDD(4) * t158 - t164 * t161 - t208 * t87, -t158 * t43 + t161 * t42 + (t158 * t87 - t161 * t85) * qJD(4), t16 * t158 + t17 * t161 + t143 + (-t158 * t44 + t161 * t45) * qJD(4) + (-qJD(1) * t58 - t209) * t154, 0, 0, 0, 0, 0, 0, t148 * t98 - t151 * t51 - t208 * t235, -t148 * t99 - t151 * t50 - t179 * t208, t179 * t51 - t235 * t50 + t8 * t98 - t9 * t99, t1 * t99 - t10 * t51 + t11 * t50 + t2 * t98 + t143 + (-qJD(1) * t31 - t209) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t228 - t229, -t180, -t222, (t87 - t193) * qJD(4) + t172, qJDD(4), -g(2) * t73 - t58 * t87 + t184 + t234, g(1) * t76 + g(2) * t74 + g(3) * t96 + t58 * t85 + (t44 + t218) * qJD(4) + t196, 0, 0, t223, t240, t239, -t223, t230, t148, -t12 * t151 + (t148 * t160 - t151 * t201 - t235 * t87) * pkin(4) + t231, t13 * t151 + (-t148 * t157 - t151 * t200 - t179 * t87) * pkin(4) + t238, -t10 * t235 + t11 * t179 + t12 * t179 + t13 * t235 + (-t157 * t9 + t160 * t8 + (t157 * t179 - t160 * t235) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t157 + t2 * t160 - t31 * t87 - t96 * t145 + (-t10 * t157 + t11 * t160) * qJD(5) + t234) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, t240, t239, -t223, t230, t148, t11 * t151 + t231, t10 * t151 + t238, 0, 0;];
tau_reg = t5;
