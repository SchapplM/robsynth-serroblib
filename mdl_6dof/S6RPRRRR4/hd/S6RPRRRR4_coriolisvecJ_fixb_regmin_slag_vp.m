% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:54
% EndTime: 2019-03-09 07:06:04
% DurationCPUTime: 4.31s
% Computational Cost: add. (10176->340), mult. (27580->452), div. (0->0), fcn. (23052->10), ass. (0->202)
t160 = cos(qJ(6));
t218 = qJD(6) * t160;
t157 = sin(qJ(5));
t161 = cos(qJ(5));
t155 = cos(pkin(11));
t163 = cos(qJ(3));
t227 = t163 * t155;
t154 = sin(pkin(11));
t159 = sin(qJ(3));
t228 = t159 * t154;
t181 = -t227 + t228;
t130 = t181 * qJD(1);
t136 = t163 * t154 + t159 * t155;
t131 = t136 * qJD(1);
t158 = sin(qJ(4));
t162 = cos(qJ(4));
t185 = -t158 * t130 + t162 * t131;
t186 = -t162 * t130 - t158 * t131;
t70 = t157 * t185 - t161 * t186;
t294 = t160 * t70;
t298 = t218 + t294;
t223 = qJD(4) * t158;
t217 = t163 * qJD(3);
t143 = t155 * qJD(1) * t217;
t212 = qJD(1) * t228;
t124 = -qJD(3) * t212 + t143;
t176 = t136 * qJD(2);
t175 = qJD(1) * t176;
t250 = pkin(7) + qJ(2);
t140 = t250 * t154;
t137 = qJD(1) * t140;
t141 = t250 * t155;
t138 = qJD(1) * t141;
t183 = t159 * t137 - t163 * t138;
t82 = -t124 * pkin(8) + qJD(3) * t183 - t175;
t98 = -t130 * pkin(8) - t183;
t203 = t158 * t82 - t98 * t223;
t133 = t136 * qJD(3);
t125 = qJD(1) * t133;
t268 = -t163 * t137 - t159 * t138;
t81 = -t125 * pkin(8) - qJD(2) * t130 + t268 * qJD(3);
t97 = -t131 * pkin(8) + t268;
t96 = qJD(3) * pkin(3) + t97;
t261 = (qJD(4) * t96 + t81) * t162 + t203;
t67 = qJD(4) * t185 + t158 * t124 + t162 * t125;
t17 = -t67 * pkin(9) + t261;
t95 = t162 * t98;
t189 = -t158 * t96 - t95;
t204 = -t158 * t81 + t162 * t82;
t170 = t189 * qJD(4) + t204;
t222 = qJD(4) * t162;
t66 = t162 * t124 - t158 * t125 - t130 * t222 - t131 * t223;
t18 = -t66 * pkin(9) + t170;
t221 = qJD(5) * t157;
t276 = pkin(9) * t186;
t49 = -t189 + t276;
t206 = t157 * t18 - t49 * t221;
t153 = qJD(3) + qJD(4);
t93 = t158 * t98;
t202 = t162 * t96 - t93;
t277 = pkin(9) * t185;
t48 = t202 - t277;
t46 = t153 * pkin(4) + t48;
t2 = (qJD(5) * t46 + t17) * t161 + t206;
t146 = -t155 * pkin(2) - pkin(1);
t139 = t146 * qJD(1) + qJD(2);
t115 = t130 * pkin(3) + t139;
t83 = -pkin(4) * t186 + t115;
t286 = t83 * t70;
t297 = t286 - t2;
t150 = qJD(5) + t153;
t237 = t70 * t150;
t220 = qJD(5) * t161;
t34 = -t157 * t67 + t161 * t66 - t185 * t221 + t186 * t220;
t280 = t34 + t237;
t156 = sin(qJ(6));
t219 = qJD(6) * t156;
t260 = t157 * t186 + t161 * t185;
t14 = t150 * t218 + t160 * t34 - t219 * t260;
t59 = t156 * t150 + t160 * t260;
t15 = qJD(6) * t59 + t156 * t34;
t57 = -t160 * t150 + t156 * t260;
t296 = t14 * t160 - t156 * t15 - t298 * t57;
t12 = t14 * t156;
t291 = t298 * t59 + t12;
t35 = qJD(5) * t260 + t157 * t66 + t161 * t67;
t31 = t156 * t35;
t293 = -qJD(6) - t70;
t61 = t293 * t218;
t248 = t31 - t61;
t253 = t59 * t260;
t290 = -t293 * t294 + t248 - t253;
t243 = t157 * t49;
t23 = t161 * t46 - t243;
t21 = -t150 * pkin(5) - t23;
t287 = t21 * t70;
t274 = t83 * t260;
t207 = t157 * t17 - t161 * t18;
t239 = t161 * t49;
t24 = t157 * t46 + t239;
t3 = t24 * qJD(5) + t207;
t295 = -t274 - t3;
t284 = t260 * t70;
t238 = t260 * t150;
t262 = -t35 + t238;
t292 = qJD(6) + t293;
t281 = t260 ^ 2 - t70 ^ 2;
t45 = pkin(5) * t260 + t70 * pkin(10);
t252 = t260 * t57;
t283 = t156 * t293;
t33 = t160 * t35;
t289 = -t283 * t293 + t252 + t33;
t288 = t283 * t59 + t296;
t275 = t293 * t260;
t233 = t153 * t185;
t282 = -t67 + t233;
t22 = t150 * pkin(10) + t24;
t38 = t70 * pkin(5) - pkin(10) * t260 + t83;
t9 = t156 * t38 + t160 * t22;
t265 = t3 * t156 + t21 * t218 + t260 * t9;
t190 = t156 * t22 - t160 * t38;
t266 = t190 * t260 + t21 * t219;
t257 = t185 * pkin(4);
t273 = -t219 * t293 - t33;
t232 = t153 * t186;
t272 = t66 - t232;
t271 = t115 * t185;
t270 = t115 * t186;
t269 = t185 * t186;
t267 = t185 ^ 2 - t186 ^ 2;
t104 = -t136 * pkin(8) - t163 * t140 - t159 * t141;
t182 = t159 * t140 - t163 * t141;
t105 = -pkin(8) * t181 - t182;
t114 = t162 * t136 - t158 * t181;
t52 = -t114 * pkin(9) + t162 * t104 - t158 * t105;
t184 = -t158 * t136 - t162 * t181;
t188 = -t158 * t104 - t162 * t105;
t53 = pkin(9) * t184 - t188;
t37 = t157 * t52 + t161 * t53;
t132 = t181 * qJD(3);
t77 = qJD(4) * t184 - t162 * t132 - t158 * t133;
t78 = qJD(4) * t114 - t158 * t132 + t162 * t133;
t79 = t157 * t114 - t161 * t184;
t39 = -qJD(5) * t79 - t157 * t78 + t161 * t77;
t172 = -t140 * t217 + qJD(2) * t227 + (-qJD(2) * t154 - qJD(3) * t141) * t159;
t87 = -t133 * pkin(8) + t172;
t165 = qJD(3) * t182 - t176;
t88 = t132 * pkin(8) + t165;
t178 = t104 * t222 - t105 * t223 + t158 * t88 + t162 * t87;
t29 = -t78 * pkin(9) + t178;
t169 = qJD(4) * t188 - t158 * t87 + t162 * t88;
t30 = -t77 * pkin(9) + t169;
t36 = t157 * t53 - t161 * t52;
t4 = -t36 * qJD(5) + t157 * t30 + t161 * t29;
t80 = t161 * t114 + t157 * t184;
t119 = pkin(3) * t181 + t146;
t91 = -pkin(4) * t184 + t119;
t43 = t79 * pkin(5) - t80 * pkin(10) + t91;
t259 = (qJD(6) * t43 + t4) * t293 - (qJD(6) * t38 + t2) * t79 + t21 * t39 + t3 * t80 - t37 * t35;
t256 = t131 * pkin(3);
t255 = t21 * t80;
t254 = t43 * t35;
t251 = t80 * t35;
t247 = t162 * t97 - t93;
t245 = t156 * t59;
t149 = t162 * pkin(3) + pkin(4);
t231 = t157 * t158;
t201 = -t158 * t97 - t95;
t50 = t201 - t276;
t51 = t247 - t277;
t236 = t157 * t50 + t161 * t51 - t149 * t220 - (-t158 * t221 + (t161 * t162 - t231) * qJD(4)) * pkin(3);
t230 = t158 * t161;
t235 = -t157 * t51 + t161 * t50 + t149 * t221 + (t158 * t220 + (t157 * t162 + t230) * qJD(4)) * pkin(3);
t224 = t154 ^ 2 + t155 ^ 2;
t216 = qJD(1) * qJD(2);
t213 = t80 * t219;
t211 = -pkin(3) * t153 - t96;
t210 = -pkin(4) * t150 - t46;
t196 = t224 * qJD(1) ^ 2;
t127 = pkin(3) * t230 + t157 * t149 + pkin(10);
t89 = t256 + t257;
t195 = qJD(6) * t127 + t45 + t89;
t147 = t157 * pkin(4) + pkin(10);
t194 = qJD(6) * t147 + t257 + t45;
t25 = t157 * t48 + t239;
t192 = pkin(4) * t221 - t25;
t54 = t125 * pkin(3) + t67 * pkin(4);
t60 = t133 * pkin(3) + t78 * pkin(4);
t191 = -t293 * t39 + t251;
t180 = t283 * t70 - t273;
t179 = 0.2e1 * t224 * t216;
t174 = -t127 * t35 - t236 * t293 + t287;
t26 = t161 * t48 - t243;
t167 = -t147 * t35 + t287 - (-pkin(4) * t220 + t26) * t293;
t148 = -t161 * pkin(4) - pkin(5);
t126 = pkin(3) * t231 - t161 * t149 - pkin(5);
t40 = qJD(5) * t80 + t157 * t77 + t161 * t78;
t10 = t40 * pkin(5) - t39 * pkin(10) + t60;
t7 = t35 * pkin(5) - t34 * pkin(10) + t54;
t6 = t160 * t7;
t5 = t37 * qJD(5) + t157 * t29 - t161 * t30;
t1 = [0, 0, 0, 0, 0, t179, qJ(2) * t179, t124 * t136 - t131 * t132, -t124 * t181 - t136 * t125 + t132 * t130 - t131 * t133, -t132 * qJD(3), -t133 * qJD(3), 0, qJD(3) * t165 + t146 * t125 + t139 * t133, -qJD(3) * t172 + t146 * t124 - t139 * t132, t66 * t114 + t185 * t77, -t114 * t67 + t184 * t66 - t185 * t78 + t186 * t77, t153 * t77, -t153 * t78, 0, t119 * t67 + t115 * t78 + t169 * t153 + (-t125 * t184 - t133 * t186) * pkin(3), t119 * t66 + t115 * t77 - t178 * t153 + (t125 * t114 + t133 * t185) * pkin(3), t260 * t39 + t34 * t80, -t260 * t40 - t34 * t79 - t39 * t70 - t251, t39 * t150, -t40 * t150, 0, -t5 * t150 + t91 * t35 + t83 * t40 + t54 * t79 + t60 * t70, -t4 * t150 + t260 * t60 + t91 * t34 + t83 * t39 + t54 * t80, -t59 * t213 + (t14 * t80 + t39 * t59) * t160 (-t160 * t57 - t245) * t39 + (-t12 - t15 * t160 + (t156 * t57 - t160 * t59) * qJD(6)) * t80, t14 * t79 + t160 * t191 + t213 * t293 + t59 * t40, -t15 * t79 - t156 * t191 - t57 * t40 + t80 * t61, -t293 * t40 + t35 * t79, t36 * t15 - t190 * t40 + t5 * t57 + t6 * t79 + (-t10 * t293 + t254 + (-t22 * t79 + t293 * t37 + t255) * qJD(6)) * t160 + t259 * t156, t36 * t14 - t9 * t40 + t5 * t59 + ((-qJD(6) * t37 + t10) * t293 - t254 - (-qJD(6) * t22 + t7) * t79 - qJD(6) * t255) * t156 + t259 * t160; 0, 0, 0, 0, 0, -t196, -qJ(2) * t196, 0, 0, 0, 0, 0, 0.2e1 * t131 * qJD(3), t143 + (-t130 - t212) * qJD(3), 0, 0, 0, 0, 0, t67 + t233, t66 + t232, 0, 0, 0, 0, 0, t35 + t238, t34 - t237, 0, 0, 0, 0, 0, t180 - t252, -t160 * t293 ^ 2 - t253 - t31; 0, 0, 0, 0, 0, 0, 0, t131 * t130, -t130 ^ 2 + t131 ^ 2, t143 + (t130 - t212) * qJD(3), 0, 0, -t139 * t131 - t175, t139 * t130 + t181 * t216, -t269, t267, t272, t282, 0, t186 * t256 - t271 - t201 * t153 + (t211 * t158 - t95) * qJD(4) + t204, -t185 * t256 - t270 + t247 * t153 + (t211 * qJD(4) - t81) * t162 - t203, t284, t281, t280, t262, 0, -t235 * t150 - t89 * t70 + t295, t236 * t150 - t260 * t89 + t297, t291, t245 * t293 + t296, t290, t180 + t252, t275, t126 * t15 + t235 * t57 + (t195 * t293 - t3) * t160 + t174 * t156 + t266, t126 * t14 + t174 * t160 - t195 * t283 + t235 * t59 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t269, t267, t272, t282, 0, -t189 * t153 + t170 - t271, t202 * t153 - t261 - t270, t284, t281, t280, t262, 0, -t70 * t257 + t25 * t150 - t274 + (t210 * t157 - t239) * qJD(5) - t207, -t260 * t257 + t26 * t150 + t286 + (t210 * qJD(5) - t17) * t161 - t206, t291, t288, t290, t289, t275, t148 * t15 + t192 * t57 + (t194 * t293 - t3) * t160 + t167 * t156 + t266, t148 * t14 + t160 * t167 + t192 * t59 - t194 * t283 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t281, t280, t262, 0, t24 * t150 + t295, t23 * t150 + t297, t291, t288, t290, t289, t275, -pkin(5) * t15 - t3 * t160 + (-t156 * t23 + t160 * t45) * t293 - t24 * t57 + t156 * t287 - t248 * pkin(10) + t266, -pkin(5) * t14 - (t156 * t45 + t160 * t23) * t293 - t24 * t59 + t21 * t294 + t273 * pkin(10) + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t57, -t57 ^ 2 + t59 ^ 2, -t293 * t57 + t14, -t293 * t59 - t15, t35, -t156 * t2 - t21 * t59 - t292 * t9 + t6, -t156 * t7 - t160 * t2 + t292 * t190 + t21 * t57;];
tauc_reg  = t1;
