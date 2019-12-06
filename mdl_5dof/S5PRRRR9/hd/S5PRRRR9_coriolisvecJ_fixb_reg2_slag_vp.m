% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:11
% EndTime: 2019-12-05 17:21:23
% DurationCPUTime: 3.89s
% Computational Cost: add. (4392->386), mult. (11189->570), div. (0->0), fcn. (8187->10), ass. (0->189)
t144 = sin(qJ(3));
t148 = cos(qJ(3));
t165 = pkin(3) * t144 - pkin(8) * t148;
t112 = t165 * qJD(3);
t116 = -t148 * pkin(3) - t144 * pkin(8) - pkin(2);
t143 = sin(qJ(4));
t145 = sin(qJ(2));
t147 = cos(qJ(4));
t195 = t147 * qJD(3);
t199 = qJD(4) * t147;
t201 = qJD(4) * t143;
t140 = sin(pkin(5));
t207 = qJD(1) * t140;
t149 = cos(qJ(2));
t211 = t148 * t149;
t240 = t116 * t199 + t143 * t112 + (-t144 * t195 - t148 * t201) * pkin(7) - (t143 * t145 + t147 * t211) * t207;
t196 = t143 * qJD(3);
t262 = (-t143 * t211 + t145 * t147) * t207 - t144 * pkin(7) * t196 - t147 * t112;
t212 = t147 * t148;
t130 = pkin(7) * t212;
t160 = pkin(4) * t144 - pkin(9) * t212;
t261 = -t160 * qJD(3) - (-t130 + (pkin(9) * t144 - t116) * t143) * qJD(4) + t262;
t180 = t148 * t196;
t183 = t144 * t199;
t157 = t180 + t183;
t260 = -t157 * pkin(9) + t240;
t204 = qJD(2) * t144;
t104 = t143 * t204 - t195;
t106 = t147 * t204 + t196;
t142 = sin(qJ(5));
t146 = cos(qJ(5));
t162 = t142 * t104 - t146 * t106;
t52 = t146 * t104 + t142 * t106;
t243 = t52 * t162;
t248 = -pkin(9) - pkin(8);
t187 = qJD(4) * t248;
t194 = t148 * qJD(2);
t109 = t165 * qJD(2);
t182 = t145 * t207;
t114 = qJD(2) * pkin(7) + t182;
t141 = cos(pkin(5));
t206 = qJD(1) * t141;
t81 = -t144 * t114 + t148 * t206;
t46 = t143 * t109 + t147 * t81;
t259 = -t46 + (pkin(9) * t194 + t187) * t143;
t45 = t147 * t109 - t143 * t81;
t258 = t160 * qJD(2) - t147 * t187 + t45;
t193 = qJD(2) * qJD(3);
t257 = qJD(3) * qJD(4) + t148 * t193;
t256 = t162 ^ 2 - t52 ^ 2;
t129 = -qJD(4) + t194;
t179 = t144 * t206;
t82 = t148 * t114 + t179;
t73 = qJD(3) * pkin(8) + t82;
t181 = t149 * t207;
t84 = t116 * qJD(2) - t181;
t38 = -t143 * t73 + t147 * t84;
t24 = -t106 * pkin(9) + t38;
t20 = -t129 * pkin(4) + t24;
t39 = t143 * t84 + t147 * t73;
t25 = -t104 * pkin(9) + t39;
t234 = t146 * t25;
t10 = t142 * t20 + t234;
t205 = qJD(2) * t140;
t185 = t149 * t205;
t156 = qJD(1) * (qJD(3) * t141 + t185);
t203 = qJD(3) * t144;
t49 = -t114 * t203 + t148 * t156;
t80 = (t112 + t182) * qJD(2);
t14 = t143 * t80 + t147 * t49 + t84 * t199 - t73 * t201;
t200 = qJD(4) * t144;
t174 = qJD(2) * t200;
t188 = t257 * t143 + t147 * t174;
t11 = -t188 * pkin(9) + t14;
t132 = t144 * t193;
t15 = -qJD(4) * t39 - t143 * t49 + t147 * t80;
t77 = t143 * t174 - t257 * t147;
t8 = pkin(4) * t132 + t77 * pkin(9) + t15;
t2 = -t10 * qJD(5) - t142 * t11 + t146 * t8;
t72 = -qJD(3) * pkin(3) - t81;
t48 = t104 * pkin(4) + t72;
t255 = t48 * t162 + t2;
t125 = -qJD(5) + t129;
t197 = qJD(5) * t146;
t198 = qJD(5) * t142;
t18 = t104 * t197 + t106 * t198 + t142 * t188 + t146 * t77;
t254 = -t125 * t52 - t18;
t1 = (qJD(5) * t20 + t11) * t146 + t142 * t8 - t25 * t198;
t253 = t48 * t52 - t1;
t153 = t162 * qJD(5) + t142 * t77 - t146 * t188;
t252 = t125 * t162 + t153;
t251 = t38 * t129 + t14;
t250 = t39 * t129 - t15;
t86 = t143 * t116 + t130;
t178 = t148 * t195;
t184 = t143 * t200;
t249 = t178 - t184;
t191 = qJD(4) + qJD(5);
t103 = t147 * t116;
t213 = t144 * t147;
t56 = -pkin(9) * t213 + t103 + (-pkin(7) * t143 - pkin(4)) * t148;
t215 = t143 * t144;
t64 = -pkin(9) * t215 + t86;
t28 = t142 * t56 + t146 * t64;
t247 = t28 * qJD(5) + t260 * t142 + t261 * t146;
t27 = -t142 * t64 + t146 * t56;
t246 = -t27 * qJD(5) + t261 * t142 - t260 * t146;
t245 = pkin(4) * t143;
t202 = qJD(3) * t148;
t50 = t114 * t202 + t144 * t156;
t219 = t140 * t145;
t94 = -t141 * t148 + t144 * t219;
t244 = t50 * t94;
t120 = t248 * t143;
t121 = t248 * t147;
t76 = t142 * t120 - t146 * t121;
t242 = t76 * qJD(5) + t259 * t142 + t258 * t146;
t75 = t146 * t120 + t142 * t121;
t241 = -t75 * qJD(5) + t258 * t142 - t259 * t146;
t239 = -t86 * qJD(4) - t262;
t216 = t142 * t143;
t107 = -t146 * t147 + t216;
t238 = -t107 * t194 - t146 * t199 - t147 * t197 + t191 * t216;
t108 = t142 * t147 + t146 * t143;
t59 = t191 * t108;
t237 = -t108 * t194 + t59;
t236 = qJD(2) * pkin(2);
t235 = t142 * t25;
t233 = t147 * t72;
t230 = t50 * t143;
t229 = t50 * t144;
t228 = t50 * t147;
t227 = t72 * t143;
t226 = t77 * t143;
t224 = t104 * t129;
t223 = t106 * t104;
t222 = t106 * t129;
t221 = t129 * t143;
t220 = t129 * t147;
t218 = t140 * t149;
t151 = qJD(2) ^ 2;
t217 = t140 * t151;
t214 = t143 * t148;
t150 = qJD(3) ^ 2;
t210 = t150 * t144;
t209 = t150 * t148;
t138 = t144 ^ 2;
t139 = t148 ^ 2;
t208 = t138 - t139;
t190 = t145 * t217;
t189 = t144 * t151 * t148;
t186 = t145 * t205;
t171 = t144 * t185;
t170 = t148 * t185;
t169 = t144 * t181;
t168 = t148 * t132;
t167 = pkin(4) * t201 - t179 - (qJD(2) * t245 + t114) * t148;
t166 = t188 * t147;
t115 = -t181 - t236;
t164 = -t115 - t181;
t95 = t141 * t144 + t148 * t219;
t159 = t143 * t218 - t95 * t147;
t62 = -t95 * t143 - t147 * t218;
t29 = t142 * t159 + t146 * t62;
t30 = t142 * t62 - t146 * t159;
t163 = -t143 * t39 - t147 * t38;
t161 = qJD(2) * t138 - t129 * t148;
t158 = qJD(2) * t164;
t154 = qJD(3) * (-t164 - t236);
t152 = t229 + t49 * t148 + (-t144 * t82 - t148 * t81) * qJD(3);
t136 = -t147 * pkin(4) - pkin(3);
t113 = (pkin(7) + t245) * t144;
t90 = t107 * t144;
t89 = t108 * t144;
t85 = -pkin(7) * t214 + t103;
t83 = t157 * pkin(4) + pkin(7) * t202;
t61 = t95 * qJD(3) + t171;
t60 = -t94 * qJD(3) + t170;
t34 = -t198 * t215 + (t191 * t213 + t180) * t146 + t249 * t142;
t33 = t142 * t180 + t59 * t144 - t146 * t178;
t32 = t188 * pkin(4) + t50;
t23 = t62 * qJD(4) + t143 * t186 + t60 * t147;
t22 = t159 * qJD(4) - t60 * t143 + t147 * t186;
t13 = t146 * t24 - t235;
t12 = -t142 * t24 - t234;
t9 = t146 * t20 - t235;
t4 = -t30 * qJD(5) - t142 * t23 + t146 * t22;
t3 = t29 * qJD(5) + t142 * t22 + t146 * t23;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t149 * t217, 0, 0, 0, 0, 0, 0, 0, 0, -t148 * t190 + (-t61 - t171) * qJD(3), t144 * t190 + (-t60 - t170) * qJD(3), (t144 * t61 + t148 * t60 + (-t144 * t95 + t148 * t94) * qJD(3)) * qJD(2), t49 * t95 + t244 + t82 * t60 - t81 * t61 + (t115 - t181) * t186, 0, 0, 0, 0, 0, 0, t61 * t104 - t22 * t129 + t132 * t62 + t188 * t94, t61 * t106 + t23 * t129 + t132 * t159 - t94 * t77, -t23 * t104 - t22 * t106 + t159 * t188 + t62 * t77, -t14 * t159 + t15 * t62 + t38 * t22 + t39 * t23 + t72 * t61 + t244, 0, 0, 0, 0, 0, 0, -t4 * t125 + t132 * t29 - t153 * t94 + t61 * t52, t3 * t125 - t132 * t30 - t162 * t61 - t94 * t18, t153 * t30 + t162 * t4 + t29 * t18 - t3 * t52, t1 * t30 + t10 * t3 + t2 * t29 + t32 * t94 + t9 * t4 + t48 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t168, -0.2e1 * t208 * t193, t209, -0.2e1 * t168, -t210, 0, -pkin(7) * t209 + t144 * t154, pkin(7) * t210 + t148 * t154, (-t138 - t139) * qJD(2) * t181 + t152, ((t144 * t81 - t148 * t82) * t149 + (-t115 - t236) * t145) * t207 + t152 * pkin(7), t249 * t106 - t77 * t213, (-t147 * t104 - t106 * t143) * t202 + (-t166 + t226 + (t104 * t143 - t106 * t147) * qJD(4)) * t144, t129 * t184 + t77 * t148 + (t106 * t144 + t161 * t147) * qJD(3), t157 * t104 + t188 * t215, t129 * t183 + t188 * t148 + (-t104 * t144 - t161 * t143) * qJD(3), (-t129 - t194) * t203, -t239 * t129 + (-t15 + (pkin(7) * t104 + t227) * qJD(3)) * t148 + (pkin(7) * t188 + t230 + t72 * t199 - t104 * t181 + (t85 * qJD(2) + t38) * qJD(3)) * t144, t240 * t129 + (t14 + (pkin(7) * t106 + t233) * qJD(3)) * t148 + (-t106 * t181 - t72 * t201 - pkin(7) * t77 + t228 + (-qJD(2) * t86 - t39) * qJD(3)) * t144, -t86 * t188 + t85 * t77 - t239 * t106 - t240 * t104 + t163 * t202 + (-t14 * t143 - t15 * t147 + (t143 * t38 - t147 * t39) * qJD(4)) * t144, -t72 * t169 + t14 * t86 + t15 * t85 + t240 * t39 + t239 * t38 + (t202 * t72 + t229) * pkin(7), t162 * t33 + t18 * t90, -t153 * t90 + t162 * t34 + t18 * t89 + t33 * t52, t33 * t125 + t18 * t148 + (-qJD(2) * t90 - t162) * t203, -t153 * t89 + t52 * t34, t34 * t125 - t153 * t148 + (-qJD(2) * t89 - t52) * t203, (-t125 - t194) * t203, -t113 * t153 - t2 * t148 + t32 * t89 + t48 * t34 + t83 * t52 + t247 * t125 + (-t52 * t181 + (qJD(2) * t27 + t9) * qJD(3)) * t144, t1 * t148 - t113 * t18 - t32 * t90 - t48 * t33 - t83 * t162 - t246 * t125 + (t162 * t181 + (-qJD(2) * t28 - t10) * qJD(3)) * t144, -t1 * t89 - t10 * t34 + t153 * t28 - t162 * t247 + t27 * t18 + t2 * t90 + t246 * t52 + t9 * t33, t1 * t28 + t32 * t113 + t2 * t27 - t247 * t9 + (t83 - t169) * t48 - t246 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t208 * t151, 0, t189, 0, 0, t144 * t158, t148 * t158, 0, 0, -t106 * t220 - t226, (-t77 + t224) * t147 + (-t188 + t222) * t143, -t129 * t199 + (t129 * t212 + (-t106 + t196) * t144) * qJD(2), -t104 * t221 - t166, t129 * t201 + (-t129 * t214 + (t104 + t195) * t144) * qJD(2), t129 * t204, -pkin(3) * t188 - t228 + t45 * t129 - t82 * t104 + (pkin(8) * t220 + t227) * qJD(4) + (-t38 * t144 + (-pkin(8) * t203 - t148 * t72) * t143) * qJD(2), pkin(3) * t77 - t82 * t106 - t46 * t129 + t230 + (-pkin(8) * t221 + t233) * qJD(4) + (-t72 * t212 + (-pkin(8) * t195 + t39) * t144) * qJD(2), t46 * t104 + t45 * t106 + ((t106 * qJD(4) - t188) * pkin(8) + t251) * t147 + ((t104 * qJD(4) - t77) * pkin(8) + t250) * t143, -t50 * pkin(3) - t38 * t45 - t39 * t46 - t72 * t82 + (qJD(4) * t163 + t14 * t147 - t15 * t143) * pkin(8), -t18 * t108 + t162 * t238, t18 * t107 + t108 * t153 + t162 * t237 + t238 * t52, t238 * t125 + (qJD(3) * t108 + t162) * t204, -t107 * t153 + t237 * t52, t237 * t125 + (-qJD(3) * t107 + t52) * t204, t125 * t204, t32 * t107 - t136 * t153 + t167 * t52 + t237 * t48 + t242 * t125 + (qJD(3) * t75 - t9) * t204, t32 * t108 - t136 * t18 - t167 * t162 - t238 * t48 - t241 * t125 + (-qJD(3) * t76 + t10) * t204, -t1 * t107 - t10 * t237 - t2 * t108 + t153 * t76 - t162 * t242 + t75 * t18 + t238 * t9 + t241 * t52, t1 * t76 - t10 * t241 + t32 * t136 + t167 * t48 + t2 * t75 - t242 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t104 ^ 2 + t106 ^ 2, -t77 - t224, -t223, -t188 - t222, t132, -t72 * t106 - t250, t72 * t104 - t251, 0, 0, -t243, t256, t254, t243, t252, t132, t12 * t125 + (-t106 * t52 + t125 * t198 + t132 * t146) * pkin(4) + t255, -t13 * t125 + (t106 * t162 + t125 * t197 - t132 * t142) * pkin(4) + t253, -t10 * t162 - t12 * t162 + t13 * t52 - t9 * t52 + (t142 * t153 + t146 * t18 + (-t142 * t162 - t146 * t52) * qJD(5)) * pkin(4), -t10 * t13 - t9 * t12 + (t1 * t142 - t106 * t48 + t146 * t2 + (t10 * t146 - t142 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t256, t254, t243, t252, t132, -t10 * t125 + t255, -t9 * t125 + t253, 0, 0;];
tauc_reg = t5;
