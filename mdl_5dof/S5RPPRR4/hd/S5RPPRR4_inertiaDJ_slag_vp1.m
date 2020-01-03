% Calculate time derivative of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:27
% EndTime: 2020-01-03 11:30:48
% DurationCPUTime: 5.83s
% Computational Cost: add. (11955->505), mult. (12820->750), div. (0->0), fcn. (12113->10), ass. (0->224)
t194 = sin(pkin(8));
t197 = -pkin(6) - qJ(3);
t189 = -pkin(7) + t197;
t240 = t189 - t197;
t282 = t240 * t194;
t191 = pkin(9) + qJ(4);
t181 = sin(t191);
t196 = cos(pkin(8));
t182 = cos(t191);
t267 = Icges(5,4) * t182;
t141 = -Icges(5,6) * t196 + (-Icges(5,2) * t181 + t267) * t194;
t268 = Icges(5,4) * t181;
t142 = -Icges(5,5) * t196 + (Icges(5,1) * t182 - t268) * t194;
t281 = -t141 * t182 - t142 * t181;
t195 = cos(pkin(9));
t177 = t195 * pkin(3) + pkin(2);
t163 = pkin(4) * t182 + t177;
t261 = t163 * t196;
t280 = (rSges(6,3) - t189) * t194 + t261;
t193 = sin(pkin(9));
t201 = (rSges(4,3) + qJ(3)) * t194 + (rSges(4,1) * t195 - rSges(4,2) * t193 + pkin(2)) * t196;
t279 = 2 * m(5);
t278 = 2 * m(6);
t277 = m(5) / 0.2e1;
t276 = m(6) / 0.2e1;
t275 = pkin(3) * t193;
t183 = qJ(5) + t191;
t176 = cos(t183);
t199 = cos(qJ(1));
t175 = sin(t183);
t198 = sin(qJ(1));
t252 = t198 * t175;
t150 = -t176 * t199 - t196 * t252;
t251 = t198 * t176;
t151 = -t175 * t199 + t196 * t251;
t256 = t194 * t198;
t102 = t151 * rSges(6,1) + t150 * rSges(6,2) + rSges(6,3) * t256;
t254 = t196 * t199;
t152 = t175 * t254 - t251;
t209 = t176 * t254 + t252;
t213 = rSges(6,1) * t209 - t152 * rSges(6,2);
t255 = t194 * t199;
t103 = -rSges(6,3) * t255 - t213;
t56 = t102 * t255 + t103 * t256;
t274 = pkin(4) * qJD(4);
t192 = qJD(4) + qJD(5);
t257 = t192 * t194;
t266 = Icges(6,4) * t175;
t131 = (-Icges(6,2) * t176 - t266) * t257;
t138 = -Icges(6,5) * t196 + (Icges(6,1) * t176 - t266) * t194;
t130 = (-Icges(6,5) * t175 - Icges(6,6) * t176) * t257;
t265 = Icges(6,4) * t176;
t132 = (-Icges(6,1) * t175 - t265) * t257;
t224 = t194 * t176 * t132 - t196 * t130;
t137 = -Icges(6,6) * t196 + (-Icges(6,2) * t175 + t265) * t194;
t232 = t192 * t176 * t137;
t42 = (-t232 + (-t138 * t192 - t131) * t175) * t194 + t224;
t273 = t42 * t196;
t161 = t198 * t261;
t165 = pkin(4) * t181 + t275;
t259 = t177 * t196;
t206 = -t259 - t282;
t88 = t161 + (-t165 + t275) * t199 + t206 * t198;
t270 = -t102 - t88;
t139 = -t196 * rSges(6,3) + (rSges(6,1) * t176 - rSges(6,2) * t175) * t194;
t238 = qJD(1) * t198;
t230 = t194 * t238;
t106 = t150 * qJD(1) + t209 * t192;
t107 = t151 * qJD(1) + t152 * t192;
t214 = -t107 * rSges(6,1) - t106 * rSges(6,2);
t66 = rSges(6,3) * t230 - t214;
t269 = t139 * t230 + t196 * t66;
t133 = (-rSges(6,1) * t175 - rSges(6,2) * t176) * t257;
t264 = t133 * t198;
t260 = t165 * t199;
t235 = qJD(4) * t194;
t147 = (-Icges(5,2) * t182 - t268) * t235;
t258 = t181 * t147;
t253 = t198 * t165;
t250 = t198 * t181;
t249 = t198 * t182;
t246 = t163 - t177;
t124 = t246 * t194 + t240 * t196;
t248 = -t124 - t139;
t237 = qJD(1) * t199;
t247 = t165 * t238 + t237 * t261;
t174 = t198 * t275;
t245 = t177 * t254 + t174;
t244 = t197 * t230 + t237 * t275;
t243 = pkin(1) * t237 + qJ(2) * t238;
t242 = qJ(2) * t237 + qJD(2) * t198;
t241 = t199 * pkin(1) + t198 * qJ(2);
t239 = qJD(1) * t194;
t236 = qJD(3) * t194;
t229 = t194 * t237;
t108 = -t152 * qJD(1) - t151 * t192;
t109 = t209 * qJD(1) + t150 * t192;
t67 = t109 * rSges(6,1) + t108 * rSges(6,2) + rSges(6,3) * t229;
t234 = t103 * t229 + t67 * t255 + t66 * t256;
t233 = t181 * t274;
t158 = -t181 * t199 + t196 * t249;
t159 = t181 * t254 - t249;
t122 = -t159 * qJD(1) - t158 * qJD(4);
t157 = -t182 * t199 - t196 * t250;
t203 = t157 * qJD(4);
t208 = t182 * t254 + t250;
t123 = t208 * qJD(1) + t203;
t79 = t123 * rSges(5,1) + t122 * rSges(5,2) + rSges(5,3) * t229;
t116 = t158 * rSges(5,1) + t157 * rSges(5,2) + rSges(5,3) * t256;
t231 = t199 * t236 + t242;
t225 = t248 * t199;
t146 = (-Icges(5,5) * t181 - Icges(5,6) * t182) * t235;
t148 = (-Icges(5,1) * t181 - t267) * t235;
t223 = t194 * t182 * t148 - t196 * t146;
t222 = t198 * t233;
t219 = -qJD(2) * t199 + t243;
t218 = rSges(3,1) * t196 - rSges(3,2) * t194;
t217 = rSges(4,1) * t193 + rSges(4,2) * t195;
t120 = t157 * qJD(1) + t208 * qJD(4);
t204 = t159 * qJD(4);
t121 = t158 * qJD(1) + t204;
t216 = -rSges(5,1) * t121 - rSges(5,2) * t120;
t215 = rSges(5,1) * t208 - t159 * rSges(5,2);
t212 = -t194 * t197 + t259;
t170 = t198 * t236;
t210 = t170 + t219;
t101 = -Icges(6,1) * t209 + Icges(6,4) * t152 - Icges(6,5) * t255;
t60 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t230;
t62 = Icges(6,4) * t107 + Icges(6,2) * t106 + Icges(6,6) * t230;
t64 = Icges(6,1) * t107 + Icges(6,4) * t106 + Icges(6,5) * t230;
t99 = -Icges(6,4) * t209 + Icges(6,2) * t152 - Icges(6,6) * t255;
t10 = -t196 * t60 + ((-t192 * t99 + t64) * t176 + (-t101 * t192 - t62) * t175) * t194;
t136 = -Icges(6,3) * t196 + (Icges(6,5) * t176 - Icges(6,6) * t175) * t194;
t18 = t106 * t137 + t107 * t138 + t152 * t131 - t209 * t132 + (-t130 * t199 + t136 * t238) * t194;
t19 = t108 * t137 + t109 * t138 + t150 * t131 + t151 * t132 + (t130 * t198 + t136 * t237) * t194;
t100 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t256;
t96 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t256;
t98 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t256;
t40 = -t196 * t96 + (t100 * t176 - t175 * t98) * t194;
t97 = -Icges(6,5) * t209 + Icges(6,6) * t152 - Icges(6,3) * t255;
t41 = -t196 * t97 + (t101 * t176 - t175 * t99) * t194;
t52 = t136 * t256 + t137 * t150 + t138 * t151;
t53 = -t136 * t255 + t152 * t137 - t138 * t209;
t61 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t229;
t63 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t229;
t65 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t229;
t9 = -t196 * t61 + ((-t192 * t98 + t65) * t176 + (-t100 * t192 - t63) * t175) * t194;
t207 = (t19 + t9) * t256 / 0.2e1 - (t10 + t18) * t255 / 0.2e1 + ((t41 + t53) * t198 + (t40 + t52) * t199) * t239 / 0.2e1;
t24 = t100 * t151 + t150 * t98 + t96 * t256;
t25 = t101 * t151 + t150 * t99 + t97 * t256;
t26 = -t100 * t209 + t152 * t98 - t96 * t255;
t27 = -t101 * t209 + t152 * t99 - t97 * t255;
t205 = -t196 * (-t273 + (-t10 * t199 + t198 * t9 + (t198 * t41 + t199 * t40) * qJD(1)) * t194) - (-t18 * t196 + ((t107 * t100 + t106 * t98 + t152 * t63 - t209 * t65) * t198 + t26 * t237 - (t107 * t101 + t106 * t99 + t152 * t62 - t209 * t64) * t199 + t27 * t238 + ((-t199 * t61 + t96 * t238) * t198 - (-t199 * t60 + t97 * t238) * t199) * t194) * t194) * t255 + (-t19 * t196 + ((t109 * t100 + t108 * t98 + t150 * t63 + t151 * t65) * t198 + t24 * t237 - (t109 * t101 + t108 * t99 + t150 * t62 + t151 * t64) * t199 + t25 * t238 + ((t198 * t61 + t96 * t237) * t198 - (t198 * t60 + t97 * t237) * t199) * t194) * t194) * t256 + (-t53 * t196 + (t198 * t26 - t199 * t27) * t194) * t230 + (-t52 * t196 + (t198 * t24 - t199 * t25) * t194) * t229;
t202 = t198 * rSges(3,3) + t218 * t199;
t200 = t217 * t198 + t201 * t199;
t190 = t194 ^ 2;
t187 = t198 * pkin(1);
t149 = (-rSges(5,1) * t181 - rSges(5,2) * t182) * t235;
t143 = -t196 * rSges(5,3) + (rSges(5,1) * t182 - rSges(5,2) * t181) * t194;
t140 = -Icges(5,3) * t196 + (Icges(5,5) * t182 - Icges(5,6) * t181) * t194;
t135 = t202 + t241;
t134 = t187 + (-rSges(3,3) - qJ(2)) * t199 + t218 * t198;
t126 = t202 * qJD(1) + t219;
t125 = (t199 * rSges(3,3) + (-pkin(1) - t218) * t198) * qJD(1) + t242;
t117 = -rSges(5,3) * t255 - t215;
t115 = -Icges(5,1) * t208 + Icges(5,4) * t159 - Icges(5,5) * t255;
t114 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t256;
t113 = -Icges(5,4) * t208 + Icges(5,2) * t159 - Icges(5,6) * t255;
t112 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t256;
t111 = -Icges(5,5) * t208 + Icges(5,6) * t159 - Icges(5,3) * t255;
t110 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t256;
t95 = t196 * t103;
t91 = t200 + t241;
t90 = t187 + (-qJ(2) - t217) * t199 + t201 * t198;
t89 = -t253 + (-t261 + t282) * t199 + t245;
t87 = t196 * t117 - t143 * t255;
t86 = -t116 * t196 - t143 * t256;
t85 = (rSges(5,3) - t197) * t255 + t215 + t241 + t245;
t84 = t187 + (-qJ(2) - t275) * t199 + t212 * t198 + t116;
t83 = t200 * qJD(1) + t210;
t82 = (t217 * t199 + (-pkin(1) - t201) * t198) * qJD(1) + t231;
t81 = -t139 * t255 + t95;
t80 = -t196 * t102 - t139 * t256;
t78 = rSges(5,3) * t230 - t216;
t77 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t229;
t76 = Icges(5,1) * t121 + Icges(5,4) * t120 + Icges(5,5) * t230;
t75 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t229;
t74 = Icges(5,4) * t121 + Icges(5,2) * t120 + Icges(5,6) * t230;
t73 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t229;
t72 = Icges(5,5) * t121 + Icges(5,6) * t120 + Icges(5,3) * t230;
t71 = pkin(4) * t203 + (t206 * t199 - t174) * qJD(1) + t247;
t70 = pkin(4) * t204 + (-t260 + (-t189 * t194 + t246 * t196) * t198) * qJD(1) + t244;
t69 = t199 * t280 + t213 + t241 + t253;
t68 = -t189 * t256 + t161 + t187 + (-qJ(2) - t165) * t199 + t102;
t59 = -t140 * t255 + t159 * t141 - t142 * t208;
t58 = t140 * t256 + t141 * t157 + t142 * t158;
t49 = (t212 * t199 + t174) * qJD(1) + t210 + t79;
t48 = (-rSges(5,3) * t194 - pkin(1) - t259) * t238 + t216 + t231 + t244;
t47 = (t281 * qJD(4) - t258) * t194 + t223;
t46 = -t196 * t79 + (-t143 * t237 - t149 * t198) * t194;
t45 = t196 * t78 + (t143 * t238 - t149 * t199) * t194;
t44 = -t196 * t111 + (-t113 * t181 + t115 * t182) * t194;
t43 = -t196 * t110 + (-t112 * t181 + t114 * t182) * t194;
t39 = -t196 * t67 + (-t139 * t237 - t264) * t194;
t38 = -t133 * t255 + t269;
t35 = -t111 * t255 + t159 * t113 - t115 * t208;
t34 = -t110 * t255 + t159 * t112 - t114 * t208;
t33 = t111 * t256 + t113 * t157 + t115 * t158;
t32 = t110 * t256 + t112 * t157 + t114 * t158;
t31 = -t196 * t222 + t170 + (-t182 * t274 - t189 * t239 - qJD(2)) * t199 + t67 + t243 + t247;
t30 = -t159 * t274 + (t260 + (-pkin(1) - t280) * t198) * qJD(1) + t214 + t231;
t29 = t194 * t225 + t196 * t89 + t95;
t28 = t270 * t196 + t248 * t256;
t23 = (t198 * t89 + t199 * t88) * t194 + t56;
t22 = t122 * t141 + t123 * t142 + t157 * t147 + t158 * t148 + (t140 * t237 + t146 * t198) * t194;
t21 = t120 * t141 + t121 * t142 + t159 * t147 - t208 * t148 + (t140 * t238 - t146 * t199) * t194;
t20 = (t198 * t78 + t199 * t79 + (-t116 * t198 + t117 * t199) * qJD(1)) * t194;
t15 = t190 * t222 + (-t67 - t71) * t196 + (qJD(1) * t225 - t264) * t194;
t14 = t190 * t199 * t233 + t196 * t70 + (t124 * t238 - t133 * t199) * t194 + t269;
t13 = -t102 * t230 + t234;
t12 = -t196 * t72 + (-t181 * t74 + t182 * t76 + (-t113 * t182 - t115 * t181) * qJD(4)) * t194;
t11 = -t196 * t73 + (-t181 * t75 + t182 * t77 + (-t112 * t182 - t114 * t181) * qJD(4)) * t194;
t4 = (t198 * t70 + t199 * t71 + (t270 * t198 + t199 * t89) * qJD(1)) * t194 + t234;
t1 = [0.2e1 * m(3) * (t125 * t135 + t126 * t134) + 0.2e1 * m(4) * (t82 * t91 + t83 * t90) + (t48 * t85 + t49 * t84) * t279 - t175 * t138 * t257 + (t30 * t69 + t31 * t68) * t278 + t223 + t224 + t281 * t235 + (-t175 * t131 - t232 - t258) * t194; m(3) * (-t199 * t125 - t198 * t126 + (-t134 * t199 + t135 * t198) * qJD(1)) + m(4) * (-t198 * t83 - t199 * t82 + (t198 * t91 - t199 * t90) * qJD(1)) + m(5) * (-t198 * t49 - t199 * t48 + (t198 * t85 - t199 * t84) * qJD(1)) + m(6) * (-t198 * t31 - t199 * t30 + (t198 * t69 - t199 * t68) * qJD(1)); 0; 0.2e1 * (m(4) * (t198 * t82 - t199 * t83 + t91 * t237 + t90 * t238) / 0.2e1 + (t198 * t48 - t199 * t49 + t85 * t237 + t84 * t238) * t277 + (t198 * t30 - t199 * t31 + t69 * t237 + t68 * t238) * t276) * t194; 0; 0; (-t47 - t42) * t196 + m(5) * (t45 * t85 + t46 * t84 + t48 * t87 + t49 * t86) + m(6) * (t14 * t69 + t15 * t68 + t28 * t31 + t29 * t30) + ((-t12 / 0.2e1 - t21 / 0.2e1) * t199 + (t11 / 0.2e1 + t22 / 0.2e1) * t198 + ((t43 / 0.2e1 + t58 / 0.2e1) * t199 + (t44 / 0.2e1 + t59 / 0.2e1) * t198) * qJD(1)) * t194 + t207; m(5) * (-t46 * t198 - t45 * t199 + (t198 * t87 - t199 * t86) * qJD(1)) + m(6) * (-t14 * t199 - t15 * t198 + (t198 * t29 - t199 * t28) * qJD(1)); 0.2e1 * (-m(5) * t20 / 0.2e1 - m(6) * t4 / 0.2e1) * t196 + 0.2e1 * ((t198 * t45 - t199 * t46 + t87 * t237 + t86 * t238) * t277 + (t14 * t198 - t15 * t199 + t29 * t237 + t28 * t238) * t276) * t194; (t87 * t45 + t86 * t46 + (t116 * t199 + t117 * t198) * t20 * t194) * t279 - t196 * (-t47 * t196 + (t11 * t198 - t12 * t199 + (t198 * t44 + t199 * t43) * qJD(1)) * t194) + (-t58 * t196 + (t198 * t32 - t199 * t33) * t194) * t229 + (-t22 * t196 + ((t122 * t112 + t123 * t114 + t157 * t75 + t158 * t77) * t198 + t32 * t237 - (t122 * t113 + t123 * t115 + t157 * t74 + t158 * t76) * t199 + t33 * t238 + ((t110 * t237 + t198 * t73) * t198 - (t111 * t237 + t198 * t72) * t199) * t194) * t194) * t256 + (-t59 * t196 + (t198 * t34 - t199 * t35) * t194) * t230 - (-t21 * t196 + ((t120 * t112 + t121 * t114 + t159 * t75 - t208 * t77) * t198 + t34 * t237 - (t120 * t113 + t121 * t115 + t159 * t74 - t208 * t76) * t199 + t35 * t238 + ((t110 * t238 - t199 * t73) * t198 - (t111 * t238 - t199 * t72) * t199) * t194) * t194) * t255 + (t14 * t29 + t15 * t28 + t23 * t4) * t278 + t205; -t273 + m(6) * (t30 * t81 + t31 * t80 + t38 * t69 + t39 * t68) + t207; m(6) * (-t39 * t198 - t38 * t199 + (t198 * t81 - t199 * t80) * qJD(1)); m(6) * (-t13 * t196 + (t198 * t38 - t199 * t39 + (t198 * t80 + t199 * t81) * qJD(1)) * t194); m(6) * (t13 * t23 + t14 * t81 + t15 * t80 + t28 * t39 + t29 * t38 + t4 * t56) + t205; (t13 * t56 + t38 * t81 + t39 * t80) * t278 + t205;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
