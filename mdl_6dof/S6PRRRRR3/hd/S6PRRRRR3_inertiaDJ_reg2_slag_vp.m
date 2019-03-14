% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:33
% EndTime: 2019-03-09 00:51:50
% DurationCPUTime: 6.76s
% Computational Cost: add. (8255->460), mult. (20864->812), div. (0->0), fcn. (20021->12), ass. (0->202)
t265 = sin(qJ(5));
t210 = t265 * qJD(5);
t211 = t265 * qJD(4);
t275 = t211 + t210;
t267 = cos(qJ(5));
t212 = t267 * qJD(5);
t274 = t267 * qJD(4) + t212;
t139 = cos(qJ(4));
t227 = t267 * t139;
t268 = -pkin(10) - pkin(9);
t136 = sin(qJ(4));
t199 = t268 * t265;
t95 = t136 * t199;
t75 = -t268 * t227 + t95;
t273 = t75 * qJD(5);
t200 = t268 * t267;
t226 = t265 * t139;
t74 = t136 * t200 + t268 * t226;
t137 = sin(qJ(3));
t244 = qJD(4) * t139;
t140 = cos(qJ(3));
t247 = qJD(3) * t140;
t272 = t136 * t247 + t137 * t244;
t94 = t267 * t136 + t226;
t215 = qJD(3) * t265;
t192 = t140 * t215;
t216 = qJD(3) * t267;
t193 = t140 * t216;
t198 = t137 * t227;
t202 = t136 * t193 + t139 * t192 + (qJD(4) + qJD(5)) * t198;
t271 = t136 * t275;
t155 = -t137 * t271 + t202;
t129 = t136 ^ 2;
t131 = t139 ^ 2;
t252 = t129 - t131;
t208 = qJD(4) * t252;
t263 = pkin(3) * t140;
t270 = t268 * t137 - pkin(2) - t263;
t115 = t136 * t211;
t269 = t136 * t210 - t274 * t139 + t115;
t201 = t274 * t136 + t275 * t139;
t266 = cos(qJ(6));
t264 = pkin(3) * t137;
t262 = t136 * pkin(8);
t128 = t137 * pkin(8);
t261 = t139 * pkin(3);
t133 = sin(pkin(6));
t141 = cos(qJ(2));
t249 = qJD(2) * t141;
t221 = t133 * t249;
t134 = cos(pkin(6));
t138 = sin(qJ(2));
t258 = t133 * t138;
t88 = t134 * t137 + t140 * t258;
t70 = t88 * qJD(3) + t137 * t221;
t87 = -t134 * t140 + t137 * t258;
t50 = t87 * t70;
t170 = t270 * t139;
t229 = -pkin(4) - t262;
t156 = t229 * t140 + t170;
t59 = t265 * t156;
t256 = t136 * t137;
t255 = t139 * t140;
t122 = pkin(8) * t255;
t191 = -pkin(9) * t137 - t263;
t181 = pkin(2) - t191;
t81 = -t136 * t181 + t122;
t73 = -pkin(10) * t256 + t81;
t38 = t267 * t73 + t59;
t260 = t137 * t70;
t259 = qJD(3) * t87;
t257 = t133 * t141;
t254 = t94 * t137;
t97 = pkin(4) * t256 + t128;
t130 = t137 ^ 2;
t251 = -t140 ^ 2 + t130;
t250 = qJD(2) * t138;
t127 = qJD(3) * t137;
t248 = qJD(3) * t139;
t246 = qJD(3) * t141;
t245 = qJD(4) * t136;
t243 = qJD(4) * t140;
t135 = sin(qJ(6));
t242 = qJD(6) * t135;
t241 = -0.2e1 * pkin(2) * qJD(3);
t240 = -0.2e1 * pkin(3) * qJD(4);
t239 = t140 * t262;
t238 = t139 * t128;
t237 = t267 * pkin(4);
t236 = t265 * pkin(4);
t235 = pkin(5) * t127;
t234 = pkin(4) * t245;
t233 = pkin(5) * t242;
t126 = pkin(8) * t247;
t232 = t268 * t140;
t231 = t87 * t245;
t77 = t272 * pkin(4) + t126;
t125 = -t139 * pkin(4) - pkin(3);
t124 = t265 * t136;
t225 = t136 * t243;
t223 = t139 * t243;
t222 = t133 * t250;
t219 = t136 * t244;
t218 = t137 * t247;
t217 = t139 * t247;
t214 = qJD(6) * t266;
t209 = -0.2e1 * t233;
t207 = t251 * qJD(3);
t206 = 0.2e1 * t218;
t205 = pkin(4) * t212;
t204 = pkin(4) * t210;
t203 = pkin(5) * t214;
t197 = t137 * t217;
t196 = t130 * t219;
t195 = t237 + pkin(5);
t194 = t266 * t265;
t60 = t267 * t156;
t37 = -t265 * t73 + t60;
t190 = -pkin(9) * t140 + t264;
t189 = t266 * t254;
t188 = t124 - t227;
t178 = t136 * t257 - t88 * t139;
t71 = -t88 * t136 - t139 * t257;
t187 = t136 * t178 - t139 * t71;
t164 = t139 * t181 + t239;
t186 = -t136 * t81 + t139 * t164;
t183 = qJD(4) * t200;
t182 = qJD(4) * t199;
t147 = (-t270 * t136 - t122) * qJD(4) + (t139 * t232 + (-t229 + t261) * t137) * qJD(3);
t148 = (t170 - t239) * qJD(4) + (-t238 + (t232 + t264) * t136) * qJD(3);
t144 = t267 * t147 - t265 * t148;
t143 = -qJD(5) * t59 - t73 * t212 + t144;
t45 = t136 * t192 + t201 * t137 - t139 * t193;
t142 = t45 * pkin(11) + t143 + t235;
t14 = -qJD(5) * t60 - t265 * t147 - t267 * t148 + t73 * t210;
t150 = -t155 * pkin(11) - t14;
t83 = -t137 * t124 + t198;
t168 = -t140 * pkin(5) - t83 * pkin(11) + t37;
t163 = t266 * t168;
t29 = -t254 * pkin(11) + t38;
t1 = -qJD(6) * t163 - t135 * t142 - t266 * t150 + t29 * t242;
t174 = t178 * t267 - t265 * t71;
t40 = t178 * t265 + t267 * t71;
t21 = t135 * t40 - t174 * t266;
t180 = t136 * t70 + t87 * t244;
t179 = -t139 * t70 + t231;
t177 = t136 * t190;
t176 = t135 * t188;
t175 = t266 * t195;
t173 = t266 * t188;
t53 = -t135 * t254 + t266 * t83;
t46 = -t74 * qJD(5) - t136 * t183 - t139 * t182;
t165 = t135 * t168;
t69 = t140 * t221 - t259;
t33 = t178 * qJD(4) - t69 * t136 + t139 * t222;
t34 = t71 * qJD(4) + t136 * t222 + t69 * t139;
t160 = t187 * qJD(4) - t136 * t33 + t139 * t34;
t48 = t164 * qJD(4) + (-t177 + t238) * qJD(3);
t49 = -t81 * qJD(4) + (pkin(8) * t256 + t139 * t190) * qJD(3);
t159 = t186 * qJD(4) - t136 * t49 - t139 * t48;
t158 = t260 + t140 * t69 + (-t137 * t88 + t140 * t87) * qJD(3);
t57 = -t266 * t205 - qJD(6) * t175 + (t265 * qJD(6) + t210) * t135 * pkin(4);
t157 = -t94 * pkin(11) + t74;
t17 = t266 * t29 + t165;
t154 = t135 * t157;
t153 = -t201 * pkin(11) - t46;
t152 = t266 * t157;
t151 = (qJD(5) + qJD(6)) * (-t267 * t135 - t194) * pkin(4);
t149 = -t135 * t150 + t266 * t142;
t146 = pkin(11) * t269 - t136 * t182 + t139 * t183 - t273;
t145 = -qJD(6) * t165 - t29 * t214 + t149;
t114 = -0.2e1 * t218;
t86 = pkin(4) * t194 + t135 * t195;
t85 = -t135 * t236 + t175;
t82 = t188 * pkin(5) + t125;
t76 = -t136 * t217 + t137 * t208;
t67 = t254 * pkin(5) + t97;
t66 = t266 * t94 - t176;
t65 = t135 * t94 + t173;
t58 = t151 - t233;
t55 = t201 * pkin(5) + t234;
t54 = -t188 * pkin(11) + t75;
t52 = t135 * t83 + t189;
t47 = -t273 + (t139 * t200 - t95) * qJD(4);
t30 = pkin(5) * t155 + t77;
t28 = t266 * t54 + t154;
t27 = -t135 * t54 + t152;
t23 = -qJD(6) * t176 - t135 * t269 + t266 * t201 + t94 * t214;
t22 = qJD(6) * t173 + t135 * t201 + t94 * t242 + t266 * t269;
t20 = t135 * t174 + t266 * t40;
t19 = t53 * qJD(6) - t135 * t45 + t266 * t155;
t18 = qJD(6) * t189 + t135 * t155 + t83 * t242 + t266 * t45;
t16 = -t135 * t29 + t163;
t15 = -t38 * qJD(5) + t144;
t13 = t174 * qJD(5) - t265 * t34 + t267 * t33;
t12 = -t178 * t210 - t71 * t212 - t265 * t33 - t267 * t34;
t8 = -qJD(6) * t154 - t135 * t153 + t266 * t146 - t54 * t214;
t7 = -qJD(6) * t152 - t135 * t146 - t266 * t153 + t54 * t242;
t4 = -t21 * qJD(6) + t135 * t12 + t266 * t13;
t3 = t266 * t12 - t135 * t13 - t174 * t242 - t40 * t214;
t2 = -qJD(6) * t17 + t149;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t133 ^ 2 * t138 * t249 + 0.2e1 * t88 * t69 + 0.2e1 * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t178 * t34 + 0.2e1 * t33 * t71 + 0.2e1 * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t174 + 0.2e1 * t13 * t40 + 0.2e1 * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t20 * t4 - 0.2e1 * t21 * t3 + 0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, -t221, 0, 0, 0, 0, 0, 0, 0, 0 (-t137 * t246 - t140 * t250) * t133 (t137 * t250 - t140 * t246) * t133, t158, -pkin(2) * t222 + pkin(8) * t158, 0, 0, 0, 0, 0, 0 (t136 * t259 - t33) * t140 + (qJD(3) * t71 + t180) * t137 (t87 * t248 + t34) * t140 + (qJD(3) * t178 - t179) * t137, t187 * t247 + (-t136 * t34 - t139 * t33 + (t136 * t71 + t139 * t178) * qJD(4)) * t137, -t33 * t164 + t34 * t81 + t48 * t178 + t49 * t71 + (t87 * t247 + t260) * pkin(8), 0, 0, 0, 0, 0, 0, -t13 * t140 + t70 * t254 + t87 * t202 + (t40 * qJD(3) - t271 * t87) * t137, -t12 * t140 + t127 * t174 - t45 * t87 + t70 * t83, t12 * t254 - t13 * t83 + t155 * t174 + t40 * t45, -t12 * t38 + t13 * t37 + t14 * t174 + t15 * t40 + t70 * t97 + t77 * t87, 0, 0, 0, 0, 0, 0, t127 * t20 - t140 * t4 + t19 * t87 + t52 * t70, -t127 * t21 - t140 * t3 - t18 * t87 + t53 * t70, t18 * t20 - t19 * t21 + t3 * t52 - t4 * t53, -t1 * t21 + t16 * t4 - t17 * t3 + t2 * t20 + t30 * t87 + t67 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, -0.2e1 * t207, 0, t114, 0, 0, t137 * t241, t140 * t241, 0, 0, 0.2e1 * t131 * t218 - 0.2e1 * t196, 0.2e1 * t130 * t208 - 0.4e1 * t136 * t197, 0.2e1 * t137 * t225 + 0.2e1 * t251 * t248, 0.2e1 * t129 * t218 + 0.2e1 * t196, -0.2e1 * t136 * t207 + 0.2e1 * t137 * t223, t114, -0.2e1 * t164 * t127 - 0.2e1 * t140 * t49 + 0.2e1 * (t130 * t244 + t136 * t206) * pkin(8), -0.2e1 * t81 * t127 - 0.2e1 * t140 * t48 + 0.2e1 * (-t130 * t245 + 0.2e1 * t197) * pkin(8), 0.2e1 * t186 * t247 + 0.2e1 * (t136 * t48 - t139 * t49 + (-t136 * t164 - t139 * t81) * qJD(4)) * t137, 0.2e1 * pkin(8) ^ 2 * t218 - 0.2e1 * t164 * t49 - 0.2e1 * t48 * t81, -0.2e1 * t83 * t45, -0.2e1 * t83 * t155 + 0.2e1 * t45 * t254, 0.2e1 * t83 * t127 + 0.2e1 * t140 * t45, 0.2e1 * t254 * t155, 0.2e1 * t202 * t140 + 0.2e1 * (-t254 * qJD(3) - t140 * t271) * t137, t114, -0.2e1 * t15 * t140 + 0.2e1 * t77 * t254 + 0.2e1 * t97 * t202 + 0.2e1 * (t37 * qJD(3) - t271 * t97) * t137, -0.2e1 * t127 * t38 - 0.2e1 * t14 * t140 - 0.2e1 * t45 * t97 + 0.2e1 * t77 * t83, 0.2e1 * t14 * t254 - 0.2e1 * t15 * t83 - 0.2e1 * t155 * t38 + 0.2e1 * t37 * t45, -0.2e1 * t14 * t38 + 0.2e1 * t15 * t37 + 0.2e1 * t77 * t97, -0.2e1 * t53 * t18, 0.2e1 * t18 * t52 - 0.2e1 * t19 * t53, 0.2e1 * t127 * t53 + 0.2e1 * t140 * t18, 0.2e1 * t52 * t19, -0.2e1 * t127 * t52 + 0.2e1 * t140 * t19, t114, 0.2e1 * t127 * t16 - 0.2e1 * t140 * t2 + 0.2e1 * t19 * t67 + 0.2e1 * t30 * t52, -0.2e1 * t1 * t140 - 0.2e1 * t127 * t17 - 0.2e1 * t18 * t67 + 0.2e1 * t30 * t53, 0.2e1 * t1 * t52 + 0.2e1 * t16 * t18 - 0.2e1 * t17 * t19 - 0.2e1 * t2 * t53, -0.2e1 * t1 * t17 + 0.2e1 * t16 * t2 + 0.2e1 * t30 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0, 0, 0, 0, 0, 0, 0, t179, t180, t160, -pkin(3) * t70 + pkin(9) * t160, 0, 0, 0, 0, 0, 0, t188 * t70 + t201 * t87, -t269 * t87 + t70 * t94, t12 * t188 - t13 * t94 + t174 * t201 + t269 * t40, pkin(4) * t231 - t12 * t75 + t125 * t70 + t13 * t74 + t174 * t46 + t40 * t47, 0, 0, 0, 0, 0, 0, t23 * t87 + t65 * t70, -t22 * t87 + t66 * t70, t20 * t22 - t21 * t23 + t3 * t65 - t4 * t66, t20 * t8 - t21 * t7 + t27 * t4 - t28 * t3 + t55 * t87 + t70 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, 0, -t127, 0, -t126, pkin(8) * t127, 0, 0, -t76, -0.4e1 * t137 * t219 - t252 * t247, t136 * t127 - t223, t76, t137 * t248 + t225, 0 (pkin(9) * t255 + (-t261 + t262) * t137) * qJD(4) + (t136 * t191 - t122) * qJD(3) (t177 + t238) * qJD(4) + (t139 * t191 + t239) * qJD(3), t159, -pkin(3) * t126 + pkin(9) * t159, -t269 * t83 - t45 * t94, -t94 * t155 + t45 * t188 - t83 * t201 + t254 * t269, t94 * t127 + t140 * t269, t155 * t188 + t254 * t201, -t188 * t127 + t201 * t140, 0, t125 * t155 + t74 * t127 - t47 * t140 + t77 * t188 + t97 * t201 + t254 * t234, -t125 * t45 - t127 * t75 - t46 * t140 + t234 * t83 - t269 * t97 + t77 * t94, t14 * t188 - t15 * t94 - t155 * t75 - t201 * t38 + t254 * t46 + t269 * t37 + t74 * t45 - t47 * t83, t125 * t77 - t14 * t75 + t15 * t74 + t234 * t97 + t37 * t47 - t38 * t46, -t18 * t66 - t22 * t53, t18 * t65 - t19 * t66 + t22 * t52 - t23 * t53, t127 * t66 + t140 * t22, t19 * t65 + t23 * t52, -t127 * t65 + t140 * t23, 0, t127 * t27 - t140 * t8 + t19 * t82 + t23 * t67 + t30 * t65 + t52 * t55, -t127 * t28 - t140 * t7 - t18 * t82 - t22 * t67 + t30 * t66 + t53 * t55, t1 * t65 + t16 * t22 - t17 * t23 + t18 * t27 - t19 * t28 - t2 * t66 + t52 * t7 - t53 * t8, -t1 * t28 + t16 * t8 - t17 * t7 + t2 * t27 + t30 * t82 + t55 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t219, -0.2e1 * t208, 0, -0.2e1 * t219, 0, 0, t136 * t240, t139 * t240, 0, 0, -0.2e1 * t94 * t269, 0.2e1 * t188 * t269 - 0.2e1 * t94 * t201, 0, 0.2e1 * t188 * t201, 0, 0, 0.2e1 * t125 * t201 + 0.2e1 * t188 * t234, -0.2e1 * t125 * t269 + 0.2e1 * t234 * t94, 0.2e1 * t188 * t46 - 0.2e1 * t201 * t75 + 0.2e1 * t269 * t74 - 0.2e1 * t47 * t94, 0.2e1 * t125 * t234 - 0.2e1 * t46 * t75 + 0.2e1 * t47 * t74, -0.2e1 * t66 * t22, 0.2e1 * t22 * t65 - 0.2e1 * t23 * t66, 0, 0.2e1 * t65 * t23, 0, 0, 0.2e1 * t23 * t82 + 0.2e1 * t55 * t65, -0.2e1 * t22 * t82 + 0.2e1 * t55 * t66, 0.2e1 * t22 * t27 - 0.2e1 * t23 * t28 + 0.2e1 * t65 * t7 - 0.2e1 * t66 * t8, 0.2e1 * t27 * t8 - 0.2e1 * t28 * t7 + 0.2e1 * t55 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, 0 (-t265 * t12 + t267 * t13 + (-t174 * t267 - t265 * t40) * qJD(5)) * pkin(4), 0, 0, 0, 0, 0, 0, t4, t3, 0, t20 * t58 - t21 * t57 - t3 * t86 + t4 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t245 + t217, 0, -t272, t127, t49, t48, 0, 0, 0, 0, -t45, 0, -t155, t127, pkin(4) * t137 * t216 + t140 * t204 + t143 (-t137 * t215 + t140 * t212) * pkin(4) + t14 (t265 * (t115 * t137 - t202) + t45 * t267 + (-t267 * t254 + (t256 * t265 + t83) * t265) * qJD(5)) * pkin(4) (-t265 * t14 + t267 * t15 + (-t265 * t37 + t267 * t38) * qJD(5)) * pkin(4), 0, 0, -t18, 0, -t19, t127, t127 * t85 - t58 * t140 + t145, -t127 * t86 - t140 * t57 + t1, t18 * t85 - t19 * t86 + t52 * t57 - t53 * t58, -t1 * t86 + t16 * t58 - t17 * t57 + t2 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, 0, -t245, 0, -pkin(9) * t244, pkin(9) * t245, 0, 0, 0, 0, -t269, 0, -t201, 0, t47, t46, -t188 * t205 - t201 * t236 + t204 * t94 + t237 * t269 (-t265 * t46 + t267 * t47 + (-t265 * t74 + t267 * t75) * qJD(5)) * pkin(4), 0, 0, -t22, 0, -t23, 0, t8, t7, t22 * t85 - t23 * t86 + t57 * t65 - t58 * t66, t27 * t58 - t28 * t57 - t7 * t86 + t8 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t204, -0.2e1 * t205, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t58, 0.2e1 * t57, 0, -0.2e1 * t57 * t86 + 0.2e1 * t58 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0 (t266 * t4 - t135 * t3 + (-t135 * t20 + t266 * t21) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t155, t127, t15, t14, 0, 0, 0, 0, -t18, 0, -t19, t127, t140 * t233 + t266 * t235 + t145 (-t127 * t135 + t140 * t214) * pkin(5) + t1 (t266 * t18 - t135 * t19 + (t135 * t53 - t266 * t52) * qJD(6)) * pkin(5) (t266 * t2 - t1 * t135 + (-t135 * t16 + t266 * t17) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, 0, -t201, 0, t47, t46, 0, 0, 0, 0, -t22, 0, -t23, 0, t8, t7 (t266 * t22 - t135 * t23 + (t135 * t66 - t266 * t65) * qJD(6)) * pkin(5) (t266 * t8 - t135 * t7 + (-t135 * t27 + t266 * t28) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, -t205, 0, 0, 0, 0, 0, 0, 0, 0, t151 + t209, -t203 + t57, 0 (t266 * t58 - t135 * t57 + (-t135 * t85 + t266 * t86) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, -0.2e1 * t203, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, t127, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t23, 0, t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, -t203, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;