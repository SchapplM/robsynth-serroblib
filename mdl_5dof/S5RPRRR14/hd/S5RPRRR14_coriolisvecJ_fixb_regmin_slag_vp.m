% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR14_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:43
% EndTime: 2019-12-31 19:19:55
% DurationCPUTime: 3.68s
% Computational Cost: add. (6366->357), mult. (21429->546), div. (0->0), fcn. (18500->12), ass. (0->181)
t140 = sin(pkin(5));
t139 = sin(pkin(6));
t231 = cos(pkin(5));
t198 = t231 * t139;
t244 = cos(qJ(3));
t172 = t244 * t198;
t142 = cos(pkin(6));
t141 = cos(pkin(11));
t206 = t244 * t141;
t192 = t142 * t206;
t255 = t140 * t192 + t172;
t138 = sin(pkin(11));
t145 = sin(qJ(3));
t225 = t142 * t145;
t164 = t138 * t225 - t206;
t158 = t140 * t164;
t110 = qJD(1) * t158;
t202 = qJD(3) * t244;
t254 = t139 * t202 + t110;
t222 = qJD(1) * t140;
t205 = t138 * t222;
t253 = t255 * qJD(1) - t145 * t205;
t90 = qJD(4) - t253;
t147 = cos(qJ(4));
t196 = qJD(1) * t231;
t204 = t141 * t222;
t124 = t139 * t204;
t213 = qJD(3) - t124;
t160 = -t142 * t196 - t213;
t111 = t147 * t160;
t144 = sin(qJ(4));
t207 = t244 * t138;
t102 = (t141 * t225 + t207) * t140 + t145 * t198;
t96 = qJD(1) * t102;
t68 = t144 * t96 + t111;
t67 = qJD(5) + t68;
t227 = t139 * t145;
t117 = -t147 * t142 + t144 * t227;
t191 = t139 * t205;
t252 = t117 * qJD(4) + t144 * t191 - t254 * t147;
t118 = t144 * t142 + t147 * t227;
t251 = t118 * qJD(4) + t254 * t144 + t147 * t191;
t210 = pkin(1) * t231;
t133 = t141 * t210;
t229 = t138 * t140;
t151 = t231 * pkin(2) + (-pkin(8) * t142 - qJ(2)) * t229;
t103 = t133 + t151;
t230 = t138 * t139;
t112 = (-pkin(2) * t141 - pkin(8) * t230 - pkin(1)) * t140;
t208 = t142 * t244;
t209 = t139 * t244;
t226 = t140 * t141;
t156 = (t142 * t226 + t198) * pkin(8);
t224 = qJ(2) * t226 + t138 * t210;
t99 = t156 + t224;
t250 = t103 * t208 + t112 * t209 - t145 * t99;
t106 = qJD(1) * t112 + qJD(2);
t188 = pkin(1) * t196;
t115 = qJ(2) * t204 + t138 * t188;
t86 = qJD(1) * t156 + t115;
t130 = t141 * t188;
t91 = qJD(1) * t151 + t130;
t249 = t106 * t209 - t145 * t86 + t91 * t208;
t248 = t140 ^ 2 * (t138 ^ 2 + t141 ^ 2);
t161 = t141 * t145 + t142 * t207;
t157 = t140 * t161;
t109 = qJD(1) * t157;
t220 = qJD(3) * t145;
t247 = -t139 * t220 + t109;
t95 = t102 * qJD(3);
t66 = t142 * t106 - t139 * t91;
t31 = -pkin(3) * t253 - t96 * pkin(9) + t66;
t79 = t244 * t86;
t46 = t106 * t227 + t91 * t225 + t79;
t33 = -pkin(9) * t160 + t46;
t15 = t144 * t31 + t147 * t33;
t153 = qJD(2) * t158;
t28 = -qJD(1) * t153 + qJD(3) * t249;
t221 = qJD(2) * t140;
t190 = t221 * t230;
t84 = t253 * qJD(3);
t85 = qJD(1) * t95;
t57 = t85 * pkin(3) - t84 * pkin(9) + qJD(1) * t190;
t199 = t144 * t28 - t147 * t57;
t6 = -t85 * pkin(4) + t15 * qJD(4) + t199;
t70 = -t144 * t160 + t147 * t96;
t246 = (t70 * pkin(4) + t67 * pkin(10)) * t67 + t6;
t245 = (t106 * t139 + t142 * t91) * t145 + t79;
t49 = qJD(4) * t70 + t144 * t84;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t218 = qJD(4) * t144;
t48 = -qJD(4) * t111 + t147 * t84 - t96 * t218;
t52 = t143 * t90 + t146 * t70;
t20 = t52 * qJD(5) + t143 * t48 - t146 * t85;
t152 = qJD(2) * t157;
t29 = qJD(1) * t152 + qJD(3) * t245;
t10 = t49 * pkin(4) - t48 * pkin(10) + t29;
t12 = t90 * pkin(10) + t15;
t32 = pkin(3) * t160 - t249;
t18 = t68 * pkin(4) - t70 * pkin(10) + t32;
t182 = t143 * t12 - t146 * t18;
t216 = qJD(4) * t147;
t168 = t144 * t57 + t147 * t28 + t31 * t216 - t33 * t218;
t5 = t85 * pkin(10) + t168;
t1 = -t182 * qJD(5) + t143 * t10 + t146 * t5;
t50 = t143 * t70 - t146 * t90;
t243 = t50 * t67;
t242 = t52 * t67;
t241 = t68 * t90;
t240 = t70 * t90;
t228 = t138 * t145;
t101 = t140 * t228 - t255;
t71 = -t139 * t103 + t142 * t112;
t40 = t101 * pkin(3) - t102 * pkin(9) + t71;
t197 = t231 * t142;
t116 = t139 * t226 - t197;
t150 = t244 * t99 + (t103 * t142 + t112 * t139) * t145;
t44 = -t116 * pkin(9) + t150;
t177 = t144 * t40 + t147 * t44;
t65 = t96 * pkin(3) - pkin(9) * t253;
t239 = t144 * t65 + t147 * t249;
t238 = pkin(9) * qJD(4);
t237 = t143 * t49;
t235 = t146 * t49;
t234 = t147 * t253;
t214 = qJD(5) * t146;
t215 = qJD(5) * t143;
t19 = t143 * t85 + t146 * t48 + t90 * t214 - t70 * t215;
t233 = t19 * t143;
t232 = -t46 + t90 * (pkin(4) * t144 - pkin(10) * t147);
t219 = qJD(4) * t143;
t217 = qJD(4) * t146;
t212 = t67 * t219;
t211 = t67 * t217;
t195 = t146 * t67;
t194 = t147 * t90;
t126 = -t147 * pkin(4) - t144 * pkin(10) - pkin(3);
t193 = t96 * pkin(10) - qJD(5) * t126 + t239;
t148 = qJD(1) ^ 2;
t186 = t140 * t148 * t231;
t4 = t146 * t12 + t143 * t18;
t17 = t101 * pkin(10) + t177;
t43 = t116 * pkin(3) - t250;
t72 = t102 * t144 + t116 * t147;
t73 = t102 * t147 - t116 * t144;
t23 = t72 * pkin(4) - t73 * pkin(10) + t43;
t181 = t143 * t23 + t146 * t17;
t180 = -t143 * t17 + t146 * t23;
t14 = -t144 * t33 + t147 * t31;
t35 = t250 * qJD(3) - t153;
t94 = (t172 + (t192 - t228) * t140) * qJD(3);
t61 = t95 * pkin(3) - t94 * pkin(9) + t190;
t179 = -t144 * t35 + t147 * t61;
t178 = -t144 * t44 + t147 * t40;
t176 = t101 * t146 - t73 * t143;
t59 = t101 * t143 + t73 * t146;
t173 = (-qJ(2) * t205 + t130) * t138 - t115 * t141;
t171 = -t67 * t214 - t237;
t170 = -t67 * t215 + t235;
t169 = -pkin(9) * t85 + t90 * t32;
t167 = t144 * t61 + t147 * t35 + t40 * t216 - t44 * t218;
t166 = -0.2e1 * t196 * t221;
t165 = -t143 * t118 - t146 * t209;
t162 = -t146 * t118 + t143 * t209;
t11 = -t90 * pkin(4) - t14;
t159 = -pkin(10) * t49 + (t11 + t14) * t67;
t2 = -t4 * qJD(5) + t146 * t10 - t143 * t5;
t36 = qJD(3) * t150 + t152;
t63 = t143 * t96 + t146 * t234;
t62 = t143 * t234 - t146 * t96;
t56 = -t72 * qJD(4) + t94 * t147;
t55 = t73 * qJD(4) + t94 * t144;
t25 = t176 * qJD(5) + t95 * t143 + t56 * t146;
t24 = t59 * qJD(5) + t56 * t143 - t95 * t146;
t21 = -t96 * pkin(4) + t144 * t249 - t147 * t65;
t16 = -t101 * pkin(4) - t178;
t13 = t55 * pkin(4) - t56 * pkin(10) + t36;
t8 = -t95 * pkin(4) + t177 * qJD(4) - t179;
t7 = t95 * pkin(10) + t167;
t3 = [0, 0, 0, t138 * t166, t141 * t166, 0.2e1 * qJD(2) * qJD(1) * t248, ((t141 * t224 + (qJ(2) * t229 - t133) * t138) * qJD(1) - t173) * t221, t84 * t102 + t96 * t94, -t84 * t101 - t102 * t85 + t253 * t94 - t96 * t95, -t84 * t116 - t160 * t94, t85 * t116 + t160 * t95, 0, t36 * t160 + t29 * t116 + t71 * t85 + t66 * t95 + (qJD(1) * t101 - t253) * t190, t28 * t116 + t35 * t160 + 0.2e1 * t96 * t190 + t66 * t94 + t71 * t84, t48 * t73 + t70 * t56, -t48 * t72 - t73 * t49 - t70 * t55 - t56 * t68, t48 * t101 + t56 * t90 + t70 * t95 + t73 * t85, -t49 * t101 - t55 * t90 - t68 * t95 - t72 * t85, t85 * t101 + t90 * t95, t179 * t90 + t178 * t85 - t199 * t101 + t14 * t95 + t36 * t68 + t43 * t49 + t29 * t72 + t32 * t55 + (-t101 * t15 - t177 * t90) * qJD(4), -t168 * t101 - t15 * t95 - t167 * t90 - t177 * t85 + t29 * t73 + t32 * t56 + t36 * t70 + t43 * t48, t19 * t59 + t52 * t25, t176 * t19 - t59 * t20 - t52 * t24 - t25 * t50, t19 * t72 + t25 * t67 + t59 * t49 + t52 * t55, t176 * t49 - t20 * t72 - t24 * t67 - t50 * t55, t49 * t72 + t67 * t55, (-qJD(5) * t181 + t146 * t13 - t143 * t7) * t67 + t180 * t49 + t2 * t72 - t182 * t55 + t8 * t50 + t16 * t20 - t6 * t176 + t11 * t24, -(qJD(5) * t180 + t143 * t13 + t146 * t7) * t67 - t181 * t49 - t1 * t72 - t4 * t55 + t8 * t52 + t16 * t19 + t6 * t59 + t11 * t25; 0, 0, 0, t138 * t186, t141 * t186, -t148 * t248, t173 * t222, 0, 0, 0, 0, 0, t142 * t85 - t109 * t160 + (t160 * t220 + t205 * t253) * t139, t142 * t84 + t110 * t160 + (t160 * t202 - t205 * t96) * t139, 0, 0, 0, 0, 0, -t109 * t68 - t117 * t85 - t251 * t90 + (t68 * t220 - t244 * t49) * t139, -t109 * t70 - t118 * t85 + t252 * t90 + (t70 * t220 - t244 * t48) * t139, 0, 0, 0, 0, 0, t117 * t20 + t165 * t49 + (qJD(5) * t162 + t252 * t143 - t247 * t146) * t67 + t251 * t50, t117 * t19 + t162 * t49 + (-qJD(5) * t165 + t247 * t143 + t252 * t146) * t67 + t251 * t52; 0, 0, 0, 0, 0, 0, 0, -t96 * t253, -t253 ^ 2 + t96 ^ 2, t160 * t253 + t84, t96 * t213 + (t96 * t197 - t95) * qJD(1), 0, -t46 * t124 - t66 * t96 + (-t161 * t221 + t197 * t46) * qJD(1) + (-t245 + t46) * qJD(3), -t249 * t124 - t66 * t253 + (t164 * t221 + t197 * t249) * qJD(1), t48 * t144 + t194 * t70, (t48 - t241) * t147 + (-t49 - t240) * t144, t144 * t85 + t194 * t90 - t70 * t96, -t90 ^ 2 * t144 + t147 * t85 + t68 * t96, -t90 * t96, -pkin(3) * t49 - t14 * t96 - t46 * t68 + (-t29 + (-t65 - t238) * t90) * t147 + (t249 * t90 + t169) * t144, -pkin(3) * t48 + t239 * t90 + t15 * t96 - t46 * t70 + (t90 * t238 + t29) * t144 + t169 * t147, t19 * t146 * t144 + (-t144 * t215 + t146 * t216 - t63) * t52, t63 * t50 + t52 * t62 + (-t143 * t52 - t146 * t50) * t216 + (-t233 - t146 * t20 + (t143 * t50 - t146 * t52) * qJD(5)) * t144, -t63 * t67 + (-t19 + t211) * t147 + (t52 * t90 + t170) * t144, t62 * t67 + (t20 - t212) * t147 + (-t50 * t90 + t171) * t144, t144 * t67 * t90 - t49 * t147, t126 * t235 - t11 * t62 - t21 * t50 + (t143 * t193 + t146 * t232) * t67 + (t11 * t219 - t2 + (qJD(4) * t50 + t171) * pkin(9)) * t147 + (t11 * t214 + t6 * t143 - t90 * t182 + (t20 + t212) * pkin(9)) * t144, -t126 * t237 - t11 * t63 - t21 * t52 + (-t143 * t232 + t146 * t193) * t67 + (t11 * t217 + t1 + (qJD(4) * t52 - t170) * pkin(9)) * t147 + (-t11 * t215 + t6 * t146 - t90 * t4 + (t19 + t211) * pkin(9)) * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t68, -t68 ^ 2 + t70 ^ 2, t48 + t241, t240 - t49, t85, -t32 * t70 - t199 + (-qJD(4) + t90) * t15, t14 * t90 + t32 * t68 - t168, t195 * t52 + t233, (t19 - t243) * t146 + (-t20 - t242) * t143, t195 * t67 - t52 * t70 + t237, -t143 * t67 ^ 2 + t50 * t70 + t235, -t67 * t70, -pkin(4) * t20 + t159 * t143 - t246 * t146 - t15 * t50 + t182 * t70, -pkin(4) * t19 + t246 * t143 + t159 * t146 - t15 * t52 + t4 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t50, -t50 ^ 2 + t52 ^ 2, t19 + t243, -t20 + t242, t49, -t11 * t52 + t4 * t67 + t2, t11 * t50 - t182 * t67 - t1;];
tauc_reg = t3;
