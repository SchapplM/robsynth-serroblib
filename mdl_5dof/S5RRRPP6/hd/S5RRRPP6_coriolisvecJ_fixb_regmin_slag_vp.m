% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:07
% EndTime: 2021-01-15 22:37:17
% DurationCPUTime: 2.40s
% Computational Cost: add. (3573->347), mult. (8997->475), div. (0->0), fcn. (5798->6), ass. (0->170)
t147 = sin(pkin(8));
t148 = cos(pkin(8));
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t101 = t147 * t151 + t148 * t149;
t152 = cos(qJ(2));
t191 = t152 * qJD(1);
t90 = t101 * qJD(3);
t227 = t101 * t191 - t90;
t180 = t149 * t191;
t211 = t148 * t151;
t196 = qJD(3) * t151;
t197 = qJD(3) * t149;
t91 = -t147 * t197 + t148 * t196;
t226 = -t147 * t180 + t191 * t211 - t91;
t131 = -qJD(3) + t191;
t150 = sin(qJ(2));
t200 = qJD(1) * t150;
t181 = t149 * t200;
t193 = t151 * qJD(2);
t105 = t181 - t193;
t194 = t149 * qJD(2);
t107 = t151 * t200 + t194;
t61 = t148 * t105 + t147 * t107;
t239 = t131 * t61;
t182 = t152 * t194;
t185 = t150 * t196;
t238 = t182 + t185;
t164 = -t147 * t105 + t148 * t107;
t237 = t164 ^ 2;
t190 = qJD(1) * qJD(2);
t236 = -0.2e1 * t190;
t233 = qJ(4) + pkin(7);
t174 = qJD(3) * t233;
t157 = -t149 * qJD(4) - t151 * t174;
t207 = t151 * t152;
t160 = pkin(3) * t150 - qJ(4) * t207;
t166 = pkin(2) * t150 - pkin(7) * t152;
t108 = t166 * qJD(1);
t220 = pkin(6) * t181 + t151 * t108;
t55 = t160 * qJD(1) + t220;
t208 = t150 * t151;
t209 = t149 * t152;
t92 = t149 * t108;
t65 = t92 + (-pkin(6) * t208 - qJ(4) * t209) * qJD(1);
t192 = t151 * qJD(4);
t87 = -t149 * t174 + t192;
t230 = (-t157 + t55) * t148 + (-t65 + t87) * t147;
t119 = -qJD(2) * pkin(2) + pkin(6) * t200;
t74 = t105 * pkin(3) + qJD(4) + t119;
t23 = t61 * pkin(4) - qJ(5) * t164 + t74;
t235 = t23 * t164;
t141 = pkin(6) * t191;
t167 = -t141 + (-t180 + t197) * pkin(3);
t189 = qJD(2) * qJD(3);
t76 = t238 * qJD(1) + t149 * t189;
t144 = t150 * pkin(6);
t120 = qJD(2) * pkin(7) + t141;
t116 = -t152 * pkin(2) - t150 * pkin(7) - pkin(1);
t96 = t116 * qJD(1);
t223 = t149 * t96;
t69 = t151 * t120 + t223;
t49 = -t105 * qJ(4) + t69;
t46 = t148 * t49;
t68 = -t149 * t120 + t151 * t96;
t48 = -t107 * qJ(4) + t68;
t21 = t147 * t48 + t46;
t234 = t21 * t164;
t175 = t150 * t190;
t170 = pkin(6) * t175;
t109 = t166 * qJD(2);
t97 = qJD(1) * t109;
t221 = -t149 * t170 - t151 * t97;
t156 = -t69 * qJD(3) - t221;
t177 = t152 * t190;
t75 = -qJD(3) * t181 + (t177 + t189) * t151;
t15 = pkin(3) * t175 - t75 * qJ(4) - t107 * qJD(4) + t156;
t162 = -t120 * t197 + t149 * t97 + t96 * t196;
t155 = -t151 * t170 + t162;
t20 = -t76 * qJ(4) - t105 * qJD(4) + t155;
t3 = -t147 * t20 + t148 * t15;
t4 = t147 * t15 + t148 * t20;
t30 = t147 * t55 + t148 * t65;
t25 = qJ(5) * t200 + t30;
t53 = t147 * t157 + t148 * t87;
t232 = t25 - t53;
t231 = pkin(4) * t200 + t230;
t229 = -t30 + t53;
t134 = pkin(6) * t207;
t219 = t151 * t109 + t194 * t144;
t31 = -t150 * t192 + t160 * qJD(2) + (-t134 + (qJ(4) * t150 - t116) * t149) * qJD(3) + t219;
t225 = t149 * t109 + t116 * t196;
t36 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t208 + (-qJD(4) * t150 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t152) * t149 + t225;
t9 = t147 * t31 + t148 * t36;
t228 = t227 * pkin(4) - t226 * qJ(5) + t101 * qJD(5) - t167;
t44 = -t131 * pkin(3) + t48;
t17 = t147 * t44 + t46;
t103 = t151 * t116;
t66 = -qJ(4) * t208 + t103 + (-pkin(6) * t149 - pkin(3)) * t152;
t202 = t149 * t116 + t134;
t210 = t149 * t150;
t70 = -qJ(4) * t210 + t202;
t38 = t147 * t66 + t148 * t70;
t224 = t147 * t49;
t222 = t75 * t149;
t118 = t233 * t151;
t178 = t233 * t149;
t72 = t147 * t118 + t148 * t178;
t218 = qJD(2) * t72;
t73 = t148 * t118 - t147 * t178;
t217 = qJD(2) * t73;
t216 = t105 * t131;
t215 = t107 * t131;
t214 = t119 * t149;
t213 = t119 * t151;
t212 = t131 * t151;
t154 = qJD(1) ^ 2;
t206 = t152 * t154;
t153 = qJD(2) ^ 2;
t205 = t153 * t150;
t204 = t153 * t152;
t22 = t148 * t48 - t224;
t203 = qJD(5) - t22;
t110 = pkin(3) * t210 + t144;
t145 = t150 ^ 2;
t201 = -t152 ^ 2 + t145;
t199 = qJD(2) * t150;
t198 = qJD(2) * t152;
t195 = t119 * qJD(3);
t188 = qJ(5) * t175 + t4;
t77 = t238 * pkin(3) + pkin(6) * t198;
t139 = -t151 * pkin(3) - pkin(2);
t186 = t150 * t197;
t179 = t152 * t193;
t59 = t76 * pkin(3) + pkin(6) * t177;
t42 = t147 * t75 + t148 * t76;
t173 = pkin(1) * t236;
t172 = t105 + t193;
t171 = -t107 + t194;
t43 = -t147 * t76 + t148 * t75;
t169 = -t73 * t42 + t72 * t43 - t53 * t61;
t168 = -t21 * t131 + t3;
t165 = -t61 ^ 2 - t237;
t8 = -t147 * t36 + t148 * t31;
t16 = t148 * t44 - t224;
t37 = -t147 * t70 + t148 * t66;
t163 = qJD(1) * t145 - t131 * t152;
t161 = -t131 * t164 + t42;
t2 = -pkin(4) * t175 - t3;
t159 = -t43 - t239;
t5 = t42 * pkin(4) - t43 * qJ(5) - qJD(5) * t164 + t59;
t137 = -t148 * pkin(3) - pkin(4);
t135 = t147 * pkin(3) + qJ(5);
t100 = t147 * t149 - t211;
t83 = -t147 * t210 + t148 * t208;
t82 = t101 * t150;
t57 = t100 * pkin(4) - t101 * qJ(5) + t139;
t51 = t147 * t182 - t148 * t179 + t150 * t90;
t50 = -t147 * t179 - t148 * t182 - t150 * t91;
t45 = t82 * pkin(4) - t83 * qJ(5) + t110;
t35 = t152 * pkin(4) - t37;
t34 = -t152 * qJ(5) + t38;
t24 = t107 * pkin(3) + pkin(4) * t164 + qJ(5) * t61;
t12 = -t131 * qJ(5) + t17;
t11 = t131 * pkin(4) + qJD(5) - t16;
t10 = -t50 * pkin(4) + t51 * qJ(5) - t83 * qJD(5) + t77;
t7 = -pkin(4) * t199 - t8;
t6 = qJ(5) * t199 - t152 * qJD(5) + t9;
t1 = -t131 * qJD(5) + t188;
t13 = [0, 0, 0, 0.2e1 * t152 * t175, t201 * t236, t204, -t205, 0, -pkin(6) * t204 + t150 * t173, pkin(6) * t205 + t152 * t173, t75 * t208 + (t179 - t186) * t107, (-t105 * t151 - t107 * t149) * t198 + (-t222 - t151 * t76 + (t105 * t149 - t107 * t151) * qJD(3)) * t150, t131 * t186 - t75 * t152 + (t107 * t150 + t163 * t151) * qJD(2), t131 * t185 + t76 * t152 + (-t105 * t150 - t163 * t149) * qJD(2), (-t131 - t191) * t199, -(-t116 * t197 + t219) * t131 + (t151 * t195 + pkin(6) * t76 + (qJD(1) * t103 + t68) * qJD(2)) * t150 + ((pkin(6) * t105 + t214) * qJD(2) + (t223 + (pkin(6) * t131 + t120) * t151) * qJD(3) + t221) * t152, (-t152 * pkin(6) * t197 + t225) * t131 + t162 * t152 + (pkin(6) * t75 - t149 * t195) * t150 + ((pkin(6) * t107 + t213) * t152 + (-pkin(6) * t212 - t202 * qJD(1) - t69) * t150) * qJD(2), t110 * t42 - t8 * t131 - t3 * t152 - t74 * t50 + t59 * t82 + t77 * t61 + (qJD(1) * t37 + t16) * t199, t110 * t43 + t9 * t131 + t4 * t152 - t74 * t51 + t59 * t83 + t77 * t164 + (-qJD(1) * t38 - t17) * t199, t16 * t51 - t164 * t8 + t17 * t50 - t3 * t83 - t37 * t43 - t38 * t42 - t4 * t82 - t9 * t61, t59 * t110 + t16 * t8 + t17 * t9 + t3 * t37 + t4 * t38 + t74 * t77, t10 * t61 + t7 * t131 + t2 * t152 - t23 * t50 + t45 * t42 + t5 * t82 + (-qJD(1) * t35 - t11) * t199, -t1 * t82 - t11 * t51 + t12 * t50 + t164 * t7 + t2 * t83 - t34 * t42 + t35 * t43 - t6 * t61, -t1 * t152 - t10 * t164 - t6 * t131 + t23 * t51 - t45 * t43 - t5 * t83 + (qJD(1) * t34 + t12) * t199, t1 * t34 + t23 * t10 + t11 * t7 + t12 * t6 + t2 * t35 + t5 * t45; 0, 0, 0, -t150 * t206, t201 * t154, 0, 0, 0, t154 * pkin(1) * t150, pkin(1) * t206, -t107 * t212 + t222, (t75 + t216) * t151 + (-t76 + t215) * t149, -t131 * t196 + (t131 * t207 + t171 * t150) * qJD(1), t131 * t197 + (-t131 * t209 + t172 * t150) * qJD(1), t131 * t200, -pkin(2) * t76 + t220 * t131 + (pkin(7) * t212 + t214) * qJD(3) + ((-pkin(7) * t194 - t68) * t150 + (-t172 * pkin(6) - t214) * t152) * qJD(1), -pkin(2) * t75 - t92 * t131 + (-t149 * pkin(7) * t131 + t213) * qJD(3) + (-t119 * t207 + (-pkin(7) * t193 + t69) * t150 + (t131 * t208 + t171 * t152) * pkin(6)) * qJD(1), t59 * t100 + t139 * t42 - t227 * t74 + t167 * t61 + t230 * t131 + (-t16 - t218) * t200, t59 * t101 + t139 * t43 - t226 * t74 + t167 * t164 + t229 * t131 + (t17 - t217) * t200, -t4 * t100 - t3 * t101 + t226 * t16 + t164 * t230 + t227 * t17 + t30 * t61 + t169, t59 * t139 - t230 * t16 + t167 * t74 + t229 * t17 - t3 * t72 + t4 * t73, t5 * t100 + t57 * t42 - t228 * t61 - t227 * t23 + t231 * t131 + (t11 - t218) * t200, -t1 * t100 + t2 * t101 - t226 * t11 + t227 * t12 + t164 * t231 + t25 * t61 + t169, -t5 * t101 - t57 * t43 + t228 * t164 + t226 * t23 + t232 * t131 + (-t12 + t217) * t200, t1 * t73 + t231 * t11 - t232 * t12 + t2 * t72 - t228 * t23 + t5 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t105, -t105 ^ 2 + t107 ^ 2, t75 - t216, -t215 - t76, t175, -t119 * t107 - t69 * t131 + t156, t119 * t105 - t68 * t131 - t155, -t74 * t164 + (-t107 * t61 + t148 * t175) * pkin(3) + t168, -t22 * t131 + t74 * t61 + (-t107 * t164 - t147 * t175) * pkin(3) - t4, t17 * t164 - t234 + (-t147 * t42 - t148 * t43) * pkin(3) + (-t16 + t22) * t61, t16 * t21 - t17 * t22 + (-t107 * t74 + t147 * t4 + t148 * t3) * pkin(3), -t235 - t24 * t61 + (pkin(4) - t137) * t175 + t168, t12 * t164 - t135 * t42 + t137 * t43 - t234 + (t11 - t203) * t61, t135 * t175 - t23 * t61 + t24 * t164 + (-0.2e1 * qJD(5) + t22) * t131 + t188, t1 * t135 - t11 * t21 + t12 * t203 + t2 * t137 - t23 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t159, t165, t16 * t164 + t17 * t61 + t59, t161, t165, t159, -t11 * t164 + t12 * t61 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164 * t61 - t175, t43 - t239, -t131 ^ 2 - t237, t12 * t131 + t2 + t235;];
tauc_reg = t13;
