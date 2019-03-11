% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:47
% EndTime: 2019-03-09 21:50:04
% DurationCPUTime: 5.78s
% Computational Cost: add. (5578->466), mult. (14890->782), div. (0->0), fcn. (13619->8), ass. (0->201)
t128 = sin(qJ(3));
t117 = qJD(3) * t128;
t127 = sin(qJ(4));
t130 = cos(qJ(4));
t131 = cos(qJ(3));
t208 = qJD(4) * t131;
t185 = t130 * t208;
t254 = -t127 * t117 + t185;
t116 = qJD(4) * t130;
t212 = qJD(3) * t131;
t253 = t128 * t116 + t127 * t212;
t229 = qJ(5) * t127;
t240 = pkin(4) + qJ(6);
t159 = -t240 * t130 - t229;
t210 = qJD(4) * t127;
t125 = sin(pkin(6));
t132 = cos(qJ(2));
t215 = qJD(2) * t132;
t183 = t125 * t215;
t129 = sin(qJ(2));
t226 = t125 * t129;
t230 = cos(pkin(6));
t78 = t230 * t128 + t131 * t226;
t48 = t78 * qJD(3) + t128 * t183;
t77 = t128 * t226 - t230 * t131;
t149 = -t130 * t48 + t77 * t210;
t252 = pkin(10) * t149;
t150 = t77 * t116 + t127 * t48;
t251 = pkin(10) * t150;
t47 = t48 * qJ(5);
t72 = t77 * qJD(5);
t250 = t47 + t72;
t225 = t125 * t132;
t49 = -qJD(3) * t77 + t131 * t183;
t249 = -qJD(4) * t225 + t49;
t115 = qJD(5) * t130;
t165 = -pkin(4) * t130 - t229;
t248 = t165 * qJD(4) + t115;
t123 = t130 ^ 2;
t219 = t127 ^ 2 - t123;
t177 = t219 * qJD(4);
t187 = t127 * t208;
t213 = qJD(3) * t130;
t142 = t128 * t213 + t187;
t241 = pkin(10) * t131;
t168 = pkin(3) * t128 - t241;
t89 = t168 * qJD(3);
t242 = pkin(10) * t128;
t169 = -pkin(3) * t131 - t242;
t94 = -pkin(2) + t169;
t238 = t94 * t116 + t127 * t89;
t43 = t142 * pkin(9) - t238;
t113 = pkin(9) * t212;
t209 = qJD(4) * t128;
t188 = t127 * t209;
t173 = t253 * pkin(4) + qJ(5) * t188 + t113;
t195 = qJ(5) * t212;
t42 = (-qJD(5) * t128 - t195) * t130 + t173;
t90 = -pkin(3) + t165;
t247 = qJD(4) * (t128 * t90 + t241) - t42;
t246 = 0.2e1 * t125;
t245 = pkin(5) + pkin(10);
t244 = pkin(4) * t48;
t243 = pkin(9) * t125;
t196 = pkin(1) * t230;
t68 = -t230 * pkin(2) + pkin(8) * t226 - t132 * t196;
t38 = t77 * pkin(3) - t78 * pkin(10) + t68;
t141 = pkin(8) * t225 + t129 * t196;
t69 = t230 * pkin(9) + t141;
t70 = (-pkin(2) * t132 - pkin(9) * t129 - pkin(1)) * t125;
t239 = t128 * t70 + t131 * t69;
t40 = -pkin(10) * t225 + t239;
t15 = t127 * t38 + t130 * t40;
t216 = qJD(2) * t129;
t184 = t125 * t216;
t28 = t127 * t184 + t249 * t130 - t78 * t210;
t237 = t127 * t28;
t228 = qJ(5) * t130;
t160 = qJ(6) * t127 - t228;
t206 = qJD(6) * t127;
t30 = t160 * t212 + (t206 + (qJ(6) * qJD(4) - qJD(5)) * t130) * t128 + t173;
t236 = t127 * t30;
t50 = t127 * t78 + t130 * t225;
t235 = t127 * t50;
t27 = t78 * t116 + t249 * t127 - t130 * t184;
t234 = t130 * t27;
t233 = t130 * t28;
t232 = t130 * t30;
t51 = -t127 * t225 + t130 * t78;
t231 = t130 * t51;
t221 = t130 * t131;
t108 = pkin(9) * t221;
t66 = t127 * t94 + t108;
t227 = qJD(6) * t77;
t224 = t127 * t128;
t223 = t127 * t131;
t222 = t128 * t130;
t220 = pkin(4) * t224 + t128 * pkin(9);
t122 = t128 ^ 2;
t218 = -t131 ^ 2 + t122;
t217 = qJ(5) * qJD(5);
t214 = qJD(3) * t127;
t211 = qJD(3) * t132;
t207 = qJD(5) * t127;
t205 = qJD(6) * t131;
t204 = t131 * qJD(5);
t203 = -0.2e1 * pkin(2) * qJD(3);
t202 = -0.2e1 * pkin(3) * qJD(4);
t201 = pkin(5) * t221;
t107 = pkin(9) * t223;
t200 = 0.2e1 * t240;
t199 = pkin(4) * t117;
t198 = pkin(10) * t210;
t83 = -pkin(3) + t159;
t197 = t83 * t116;
t120 = t125 ^ 2;
t194 = t120 * t215;
t191 = qJ(5) * t117;
t190 = t130 * t212;
t182 = t127 * t116;
t181 = t128 * t212;
t14 = -t127 * t40 + t130 * t38;
t180 = -t128 * t69 + t131 * t70;
t73 = (pkin(2) * t129 - pkin(9) * t132) * t125 * qJD(2);
t178 = qJD(2) * t230;
t170 = t132 * t178;
t74 = -pkin(1) * t170 + pkin(8) * t184;
t21 = t69 * t117 - t128 * t73 + t131 * t74 - t70 * t212;
t19 = pkin(10) * t184 - t21;
t75 = t141 * qJD(2);
t26 = t48 * pkin(3) - t49 * pkin(10) + t75;
t179 = t40 * t116 + t127 * t19 - t130 * t26 + t38 * t210;
t65 = t130 * t94 - t107;
t176 = t218 * qJD(3);
t175 = 0.2e1 * t181;
t174 = t254 * pkin(9) - t130 * t89 + t94 * t210;
t172 = t130 * t181;
t171 = t129 * t194;
t12 = -qJ(5) * t77 - t15;
t39 = pkin(3) * t225 - t180;
t54 = qJ(5) * t131 - t66;
t10 = -pkin(5) * t50 - t12;
t9 = pkin(5) * t51 - t240 * t77 - t14;
t166 = -t10 * t127 + t130 * t9;
t13 = -pkin(4) * t77 - t14;
t164 = t12 * t127 + t13 * t130;
t163 = -t127 * t51 - t130 * t50;
t119 = t131 * pkin(4);
t55 = t119 - t65;
t162 = t127 * t54 + t130 * t55;
t161 = -qJ(5) * t27 - qJD(5) * t50;
t22 = -t70 * t117 + t128 * t74 + t131 * t73 - t69 * t212;
t158 = -qJD(6) * t130 - t207;
t157 = -pkin(5) * t28 - t179;
t145 = -qJ(5) * t51 + t39;
t11 = t240 * t50 + t145;
t20 = -pkin(3) * t184 - t22;
t134 = -qJ(5) * t28 - qJD(5) * t51 + t20;
t3 = qJD(6) * t50 + t240 * t27 + t134;
t156 = t11 * t116 + t127 * t3;
t155 = t11 * t210 - t130 * t3;
t16 = pkin(4) * t50 + t145;
t8 = pkin(4) * t27 + t134;
t154 = -t16 * t116 - t127 * t8;
t153 = t130 * t8 - t16 * t210;
t152 = t39 * t116 + t127 * t20;
t151 = -t130 * t20 + t39 * t210;
t111 = pkin(4) * t210;
t52 = t160 * qJD(4) + t111 + t158;
t148 = -t130 * t52 + t83 * t210;
t53 = t160 * t128 + t220;
t147 = -qJD(4) * t53 - t83 * t212;
t6 = -t38 * t116 - t127 * t26 - t130 * t19 + t40 * t210;
t144 = t128 * t211 + t131 * t216;
t143 = t128 * t216 - t131 * t211;
t4 = t6 - t250;
t140 = -pkin(5) * t188 + t174;
t139 = -pkin(5) * t27 - t6;
t5 = t179 - t244;
t138 = t164 * qJD(4) + t127 * t5 - t130 * t4;
t71 = -qJ(5) * t222 + t220;
t76 = -qJ(5) * t116 + t111 - t207;
t137 = -qJD(4) * t71 - t128 * t76 + (-t131 * t90 + t242) * qJD(3);
t35 = -t191 + t43 + t204;
t41 = t174 - t199;
t136 = t162 * qJD(4) + t127 * t41 - t130 * t35;
t135 = 0.2e1 * t191 - 0.2e1 * t204 - t43;
t133 = 0.2e1 * qJD(5);
t112 = pkin(10) * t116;
t96 = t245 * t130;
t95 = t245 * t127;
t88 = pkin(5) * t116 + t112;
t87 = t245 * t210;
t79 = -t188 + t190;
t46 = -pkin(5) * t224 - t54;
t45 = qJ(6) * t131 + t107 + t119 + (pkin(5) * t128 - t94) * t130;
t29 = -t204 + (-pkin(5) * t222 - t107) * qJD(4) + (-pkin(5) * t223 + (-pkin(9) * t130 + qJ(5)) * t128) * qJD(3) + t238;
t25 = t205 + (-t240 * t128 + t201) * qJD(3) + t140;
t2 = t139 + t250;
t1 = -t240 * t48 - t157 - t227;
t7 = [0, 0, 0, 0.2e1 * t171, 0.2e1 * (-t129 ^ 2 + t132 ^ 2) * t120 * qJD(2), t170 * t246, -0.2e1 * t178 * t226, 0, -0.2e1 * t120 * pkin(1) * t216 - 0.2e1 * t75 * t230, -0.2e1 * pkin(1) * t194 + 0.2e1 * t74 * t230, 0.2e1 * t78 * t49, -0.2e1 * t48 * t78 - 0.2e1 * t49 * t77 (-t132 * t49 + t216 * t78) * t246 (t132 * t48 - t216 * t77) * t246, -0.2e1 * t171, 0.2e1 * t68 * t48 + 0.2e1 * t75 * t77 + 0.2e1 * (-t22 * t132 + t180 * t216) * t125, 0.2e1 * t68 * t49 + 0.2e1 * t75 * t78 + 0.2e1 * (-t21 * t132 - t216 * t239) * t125, 0.2e1 * t51 * t28, -0.2e1 * t27 * t51 - 0.2e1 * t28 * t50, 0.2e1 * t28 * t77 + 0.2e1 * t48 * t51, -0.2e1 * t27 * t77 - 0.2e1 * t48 * t50, 0.2e1 * t77 * t48, 0.2e1 * t14 * t48 - 0.2e1 * t179 * t77 + 0.2e1 * t20 * t50 + 0.2e1 * t27 * t39, -0.2e1 * t15 * t48 + 0.2e1 * t20 * t51 + 0.2e1 * t28 * t39 + 0.2e1 * t6 * t77, 0.2e1 * t12 * t27 + 0.2e1 * t13 * t28 + 0.2e1 * t4 * t50 + 0.2e1 * t5 * t51, 0.2e1 * t13 * t48 - 0.2e1 * t16 * t27 + 0.2e1 * t5 * t77 - 0.2e1 * t50 * t8, -0.2e1 * t12 * t48 - 0.2e1 * t16 * t28 - 0.2e1 * t4 * t77 - 0.2e1 * t51 * t8, 0.2e1 * t12 * t4 + 0.2e1 * t13 * t5 + 0.2e1 * t16 * t8, 0.2e1 * t1 * t51 - 0.2e1 * t10 * t27 - 0.2e1 * t2 * t50 + 0.2e1 * t28 * t9, 0.2e1 * t10 * t48 - 0.2e1 * t11 * t28 + 0.2e1 * t2 * t77 - 0.2e1 * t3 * t51, -0.2e1 * t1 * t77 + 0.2e1 * t11 * t27 + 0.2e1 * t3 * t50 - 0.2e1 * t48 * t9, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, t183, -t184, 0, -t75, t74, t128 * t49 + t78 * t212, -t128 * t48 + t131 * t49 + (-t128 * t78 - t131 * t77) * qJD(3), t143 * t125, t144 * t125, 0, -pkin(2) * t48 + t117 * t68 - t131 * t75 - t143 * t243, -pkin(2) * t49 + t128 * t75 - t144 * t243 + t212 * t68, t51 * t190 + (-t210 * t51 + t233) * t128, t163 * t212 + (-t237 - t234 + (-t231 + t235) * qJD(4)) * t128 (t213 * t77 - t28) * t131 + (qJD(3) * t51 - t149) * t128 (-t214 * t77 + t27) * t131 + (-qJD(3) * t50 - t150) * t128, t117 * t77 - t131 * t48, -t174 * t77 + t48 * t65 + (t179 + (pkin(9) * t50 + t127 * t39) * qJD(3)) * t131 + (pkin(9) * t27 + qJD(3) * t14 + t152) * t128, t43 * t77 - t48 * t66 + (-t6 + (pkin(9) * t51 + t130 * t39) * qJD(3)) * t131 + (pkin(9) * t28 - qJD(3) * t15 - t151) * t128, t27 * t54 + t28 * t55 + t35 * t50 + t41 * t51 + t164 * t212 + (t127 * t4 + t130 * t5 + (t12 * t130 - t127 * t13) * qJD(4)) * t128, -t27 * t71 + t41 * t77 - t42 * t50 + t48 * t55 + (-t16 * t214 - t5) * t131 + (qJD(3) * t13 + t154) * t128, -t28 * t71 - t35 * t77 - t42 * t51 - t48 * t54 + (-t16 * t213 + t4) * t131 + (-qJD(3) * t12 - t153) * t128, t12 * t35 + t13 * t41 + t16 * t42 + t4 * t54 + t5 * t55 + t71 * t8, t25 * t51 - t27 * t46 + t28 * t45 - t29 * t50 + t166 * t212 + (t1 * t130 - t127 * t2 + (-t10 * t130 - t127 * t9) * qJD(4)) * t128, -t28 * t53 + t29 * t77 - t30 * t51 + t46 * t48 + (-t11 * t213 - t2) * t131 + (qJD(3) * t10 + t155) * t128, -t25 * t77 + t27 * t53 + t30 * t50 - t45 * t48 + (t11 * t214 + t1) * t131 + (-qJD(3) * t9 + t156) * t128, t1 * t45 + t10 * t29 + t11 * t30 + t2 * t46 + t25 * t9 + t3 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, -0.2e1 * t176, 0, 0, 0, t128 * t203, t131 * t203, -0.2e1 * t122 * t182 + 0.2e1 * t123 * t181, 0.2e1 * t122 * t177 - 0.4e1 * t127 * t172, 0.2e1 * t128 * t187 + 0.2e1 * t213 * t218, -0.2e1 * t127 * t176 + 0.2e1 * t128 * t185, -0.2e1 * t181, 0.2e1 * t65 * t117 + 0.2e1 * t131 * t174 + 0.2e1 * (t116 * t122 + t127 * t175) * pkin(9), -0.2e1 * t66 * t117 - 0.2e1 * t131 * t43 + 0.2e1 * (-t122 * t210 + 0.2e1 * t172) * pkin(9), 0.2e1 * t162 * t212 + 0.2e1 * (t127 * t35 + t130 * t41 + (-t127 * t55 + t130 * t54) * qJD(4)) * t128, 0.2e1 * (-t214 * t71 - t41) * t131 + 0.2e1 * (qJD(3) * t55 - t116 * t71 - t127 * t42) * t128, 0.2e1 * (-t213 * t71 + t35) * t131 + 0.2e1 * (-qJD(3) * t54 - t130 * t42 + t210 * t71) * t128, 0.2e1 * t35 * t54 + 0.2e1 * t41 * t55 + 0.2e1 * t42 * t71, 0.2e1 * (-t127 * t46 + t130 * t45) * t212 + 0.2e1 * (-t127 * t29 + t130 * t25 + (-t127 * t45 - t130 * t46) * qJD(4)) * t128, 0.2e1 * (-t213 * t53 - t29) * t131 + 0.2e1 * (qJD(3) * t46 + t210 * t53 - t232) * t128, 0.2e1 * (t214 * t53 + t25) * t131 + 0.2e1 * (-qJD(3) * t45 + t116 * t53 + t236) * t128, 0.2e1 * t25 * t45 + 0.2e1 * t29 * t46 + 0.2e1 * t30 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, t184, t22, t21, t116 * t51 + t237, qJD(4) * t163 - t127 * t27 + t233, t150, -t149, 0, -pkin(3) * t27 + t151 - t251, -pkin(3) * t28 + t152 + t252 (t237 - t234 + (t231 + t235) * qJD(4)) * pkin(10) + t138, -t27 * t90 - t50 * t76 + t153 + t251, -t28 * t90 - t51 * t76 + t154 - t252, pkin(10) * t138 + t16 * t76 + t8 * t90, qJD(4) * t166 + t1 * t127 + t130 * t2 - t27 * t96 + t28 * t95 + t50 * t87 + t51 * t88, -t28 * t83 + t48 * t96 - t51 * t52 - t77 * t87 - t156, t27 * t83 - t48 * t95 + t50 * t52 - t77 * t88 + t155, t1 * t95 - t10 * t87 + t11 * t52 + t2 * t96 + t3 * t83 + t88 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, -t117, 0, -t113, pkin(9) * t117, t127 * t190 - t128 * t177, -0.4e1 * t128 * t182 - t212 * t219, -t254, t142, 0 (pkin(10) * t221 + (-pkin(3) * t130 + pkin(9) * t127) * t128) * qJD(4) + (t127 * t169 - t108) * qJD(3) (pkin(9) * t222 + t127 * t168) * qJD(4) + (t130 * t169 + t107) * qJD(3), t136, t137 * t127 - t247 * t130, t247 * t127 + t137 * t130, pkin(10) * t136 + t42 * t90 + t71 * t76 (t95 * t212 + t128 * t88 + t29 + (-t128 * t96 + t45) * qJD(4)) * t130 + (-t96 * t212 + t128 * t87 + t25 + (-t128 * t95 - t46) * qJD(4)) * t127, -t236 + t131 * t87 + t147 * t130 + (qJD(3) * t96 + t148) * t128, -t232 + t131 * t88 + (-qJD(3) * t95 + t197) * t128 + (t128 * t52 - t147) * t127, t25 * t95 + t29 * t96 + t30 * t83 + t45 * t88 - t46 * t87 + t52 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t182, -0.2e1 * t177, 0, 0, 0, t127 * t202, t130 * t202, 0, 0.2e1 * t130 * t76 - 0.2e1 * t210 * t90, -0.2e1 * t116 * t90 - 0.2e1 * t127 * t76, 0.2e1 * t90 * t76, 0.2e1 * t127 * t88 - 0.2e1 * t130 * t87 + 0.2e1 * (-t127 * t96 + t130 * t95) * qJD(4), -0.2e1 * t127 * t52 - 0.2e1 * t197, 0.2e1 * t148, 0.2e1 * t52 * t83 - 0.2e1 * t87 * t96 + 0.2e1 * t88 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t48, -t179, t6, -pkin(4) * t28 + t161, t179 - 0.2e1 * t244, -t4 + t250, -pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t12, -qJD(6) * t51 - t240 * t28 + t161, t139 + 0.2e1 * t47 + 0.2e1 * t72, t200 * t48 + t157 + 0.2e1 * t227, qJ(5) * t2 + qJD(5) * t10 - qJD(6) * t9 - t1 * t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t253, t117, -t174, t43 (-pkin(4) * t212 - qJ(5) * t209) * t130 + (-t195 + (pkin(4) * qJD(4) - qJD(5)) * t128) * t127, t174 - 0.2e1 * t199, t135, -pkin(4) * t41 - qJ(5) * t35 - qJD(5) * t54, t159 * t212 + ((t127 * t240 - t228) * qJD(4) + t158) * t128, -pkin(5) * t253 + t135, -0.2e1 * t205 + (t128 * t200 - t201) * qJD(3) - t140, qJ(5) * t29 + qJD(5) * t46 - qJD(6) * t45 - t240 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, -t210, 0, -t112, t198, t248, t112, -t198, t248 * pkin(10), qJD(4) * t159 + t115 - t206, -t87, -t88, -qJ(5) * t87 + qJD(5) * t96 - t95 * qJD(6) - t240 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0.2e1 * t217, 0, t133, 0.2e1 * qJD(6), 0.2e1 * qJD(6) * t240 + 0.2e1 * t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t48, 0, t5, t28, 0, -t48, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t117, 0, t41, t79, 0, -t117, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, t112, t116, 0, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t48, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t117, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, 0, 0, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
