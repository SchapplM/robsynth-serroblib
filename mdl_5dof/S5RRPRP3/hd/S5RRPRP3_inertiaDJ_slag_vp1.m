% Calculate time derivative of joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:51:05
% DurationCPUTime: 5.14s
% Computational Cost: add. (6008->306), mult. (4872->425), div. (0->0), fcn. (3598->8), ass. (0->182)
t163 = pkin(8) + qJ(4);
t159 = cos(t163);
t158 = sin(t163);
t248 = Icges(6,5) * t158;
t250 = Icges(5,4) * t158;
t301 = t248 - t250 + (-Icges(5,2) - Icges(6,3)) * t159;
t247 = Icges(6,5) * t159;
t249 = Icges(5,4) * t159;
t300 = -t247 + t249 + (Icges(5,1) + Icges(6,1)) * t158;
t167 = cos(pkin(8));
t299 = -rSges(4,1) * t167 - pkin(2) + rSges(4,2) * sin(pkin(8));
t193 = Icges(6,3) * t158 + t247;
t196 = -Icges(5,2) * t158 + t249;
t197 = Icges(6,1) * t159 + t248;
t198 = Icges(5,1) * t159 - t250;
t298 = (t197 + t198) * t158 - (t193 - t196) * t159;
t252 = rSges(6,3) + qJ(5);
t291 = rSges(6,1) + pkin(4);
t282 = t158 * t252 + t159 * t291;
t164 = qJD(1) + qJD(2);
t165 = qJ(1) + qJ(2);
t160 = sin(t165);
t161 = cos(t165);
t289 = rSges(4,3) + qJ(3);
t66 = t289 * t160 - t299 * t161;
t295 = t164 * t66;
t294 = -t301 * t158 - t300 * t159;
t239 = t291 * t158 - t252 * t159;
t109 = Icges(5,5) * t158 + Icges(5,6) * t159;
t110 = Icges(6,4) * t158 - Icges(6,6) * t159;
t288 = t109 + t110;
t287 = t299 * t160;
t194 = Icges(5,5) * t159 - Icges(5,6) * t158;
t195 = Icges(6,4) * t159 + Icges(6,6) * t158;
t281 = t294 * t164 + (t194 + t195) * qJD(4);
t187 = t196 * t161;
t74 = Icges(5,6) * t160 + t187;
t189 = t198 * t161;
t78 = Icges(5,5) * t160 + t189;
t200 = t158 * t74 - t159 * t78;
t280 = t160 * t200;
t184 = t193 * t161;
t68 = Icges(6,6) * t160 + t184;
t188 = t197 * t161;
t76 = Icges(6,4) * t160 + t188;
t202 = t158 * t68 + t159 * t76;
t279 = t160 * t202;
t73 = -Icges(5,6) * t161 + t160 * t196;
t77 = -Icges(5,5) * t161 + t160 * t198;
t201 = t158 * t73 - t159 * t77;
t278 = t161 * t201;
t67 = -Icges(6,6) * t161 + t160 * t193;
t75 = -Icges(6,4) * t161 + t160 * t197;
t203 = t158 * t67 + t159 * t75;
t277 = t161 * t203;
t274 = 2 * m(3);
t273 = 2 * m(4);
t272 = 2 * m(5);
t271 = 2 * m(6);
t259 = rSges(5,2) * t158;
t261 = rSges(5,1) * t159;
t105 = (-t259 + t261) * qJD(4);
t267 = m(5) * t105;
t116 = rSges(5,1) * t158 + rSges(5,2) * t159;
t266 = m(5) * t116;
t169 = sin(qJ(1));
t265 = pkin(1) * t169;
t152 = t161 * rSges(6,2);
t264 = t282 * t160 - t152;
t149 = t160 * rSges(6,2);
t245 = t159 * t161;
t246 = t158 * t161;
t263 = t291 * t245 + t252 * t246 + t149;
t258 = pkin(1) * qJD(1);
t147 = t160 * rSges(5,3);
t64 = t239 * t161;
t257 = t164 * t64;
t251 = -t282 * qJD(4) + qJD(5) * t159;
t244 = t160 * t164;
t243 = t161 * t164;
t168 = -pkin(7) - qJ(3);
t242 = t161 * t168;
t241 = t164 * t168;
t131 = t160 * t259;
t240 = rSges(5,3) * t243 + t164 * t131;
t232 = qJD(4) * t160;
t224 = t158 * t232;
t238 = t291 * t224;
t144 = qJD(3) * t161;
t237 = t160 * t241 + t144;
t236 = t161 * rSges(5,3) + t131;
t235 = t160 ^ 2 + t161 ^ 2;
t234 = qJD(4) * t158;
t233 = qJD(4) * t159;
t230 = qJD(5) * t158;
t228 = rSges(5,2) * t246;
t227 = t169 * t258;
t170 = cos(qJ(1));
t226 = t170 * t258;
t223 = t159 * t232;
t225 = -rSges(5,1) * t224 - rSges(5,2) * t223 - t164 * t228;
t143 = qJD(3) * t160;
t154 = pkin(3) * t167 + pkin(2);
t183 = -t154 - t282;
t173 = t160 * t183 - t242;
t182 = rSges(6,2) * t243 + (-t239 * qJD(4) + t230) * t161;
t12 = t164 * t173 + t143 + t182;
t10 = t12 - t227;
t39 = t152 + t173;
t37 = t39 - t265;
t218 = t164 * t37 - t10;
t13 = (-t233 * t252 - t230) * t160 + (t161 * t183 - t149) * t164 + t237 + t238;
t11 = t13 - t226;
t162 = t170 * pkin(1);
t207 = t161 * t154 - t160 * t168;
t40 = t207 + t263;
t38 = t162 + t40;
t217 = t164 * t38 + t11;
t216 = t164 * t39 - t12;
t215 = t164 * t40 + t13;
t190 = t116 * qJD(4);
t208 = -t154 - t261;
t199 = t208 * t160;
t25 = t143 + t164 * t199 + (-t190 - t241) * t161 + t240;
t23 = t25 - t227;
t59 = t199 + t236 - t242;
t57 = t59 - t265;
t214 = t164 * t57 - t23;
t26 = (t161 * t208 - t147) * t164 - t225 + t237;
t24 = t26 - t226;
t82 = rSges(5,1) * t245 + t147 - t228;
t60 = t207 + t82;
t58 = t162 + t60;
t213 = t164 * t58 + t24;
t212 = t164 * t59 - t25;
t211 = t164 * t60 + t26;
t27 = t161 * t251 + t239 * t244;
t63 = t239 * t160;
t210 = -t164 * t63 + t27;
t28 = t160 * t251 - t257;
t209 = -t28 - t257;
t118 = t161 * rSges(3,1) - rSges(3,2) * t160;
t95 = -rSges(3,1) * t243 + rSges(3,2) * t244;
t117 = -rSges(3,1) * t160 - rSges(3,2) * t161;
t94 = t117 * t164;
t186 = t195 * t161;
t185 = t194 * t161;
t65 = t289 * t161 + t287;
t179 = Icges(6,2) * t164 - qJD(4) * t110;
t176 = Icges(5,3) * t164 - qJD(4) * t109;
t55 = t287 * t164 + t289 * t243 + t143;
t175 = t298 * qJD(4) + t300 * t233 + t301 * t234;
t174 = (t281 * t160 + (-t200 + t202) * qJD(4) - t298 * t244) * t160 / 0.2e1 - (-t281 * t161 + (-t201 + t203) * qJD(4) + ((-t184 + t187) * t159 + (t188 + t189) * t158) * t164) * t161 / 0.2e1 + (-t288 * t161 - t294 * t160 + (-t67 + t73) * t159 + (t75 + t77) * t158) * t244 / 0.2e1 + (-t294 * t161 + t288 * t160 + (-t68 + t74) * t159 + (t76 + t78) * t158) * t243 / 0.2e1;
t56 = t144 - t295;
t97 = t118 + t162;
t96 = t117 - t265;
t84 = t95 - t226;
t83 = t94 - t227;
t80 = t160 * t261 - t236;
t72 = Icges(6,2) * t160 + t186;
t71 = -Icges(6,2) * t161 + t160 * t195;
t70 = Icges(5,3) * t160 + t185;
t69 = -Icges(5,3) * t161 + t160 * t194;
t62 = t162 + t66;
t61 = t65 - t265;
t48 = t160 * t179 + t164 * t186;
t47 = t161 * t179 - t195 * t244;
t46 = t160 * t176 + t164 * t185;
t45 = t161 * t176 - t194 * t244;
t42 = t56 - t226;
t41 = t55 - t227;
t22 = t160 * t70 - t161 * t200;
t21 = t160 * t69 - t278;
t20 = t160 * t72 + t161 * t202;
t19 = t160 * t71 + t277;
t18 = -t161 * t70 - t280;
t17 = -t160 * t201 - t161 * t69;
t16 = -t161 * t72 + t279;
t15 = t160 * t203 - t161 * t71;
t14 = t160 * t264 + t161 * t263;
t1 = (t164 * t264 + t182) * t161 + (t160 * t230 + (t149 - t263) * t164 + t252 * t223 - t238) * t160;
t2 = [(t10 * t38 + t11 * t37) * t271 + (t23 * t58 + t24 * t57) * t272 + (t41 * t62 + t42 * t61) * t273 + (t83 * t97 + t84 * t96) * t274 + t175; m(6) * (t10 * t40 + t11 * t39 + t12 * t38 + t13 * t37) + m(5) * (t23 * t60 + t24 * t59 + t25 * t58 + t26 * t57) + m(4) * (t41 * t66 + t42 * t65 + t55 * t62 + t56 * t61) + m(3) * (t117 * t84 + t118 * t83 + t94 * t97 + t95 * t96) + t175; (t12 * t40 + t13 * t39) * t271 + (t25 * t60 + t26 * t59) * t272 + (t55 * t66 + t56 * t65) * t273 + (t117 * t95 + t118 * t94) * t274 + t175; m(6) * (t160 * t217 + t161 * t218) + m(5) * (t160 * t213 + t161 * t214) + m(4) * ((t164 * t61 - t41) * t161 + (t164 * t62 + t42) * t160); m(6) * (t160 * t215 + t161 * t216) + m(5) * (t160 * t211 + t161 * t212) + m(4) * ((t164 * t65 - t55) * t161 + (t56 + t295) * t160); 0; m(6) * (-t10 * t63 - t11 * t64 + t27 * t37 + t28 * t38) + (t160 * t214 - t161 * t213) * t266 + (-t160 * t58 - t161 * t57) * t267 + t174; m(6) * (-t12 * t63 - t13 * t64 + t27 * t39 + t28 * t40) + (-t160 * t60 - t161 * t59) * t267 + (t160 * t212 - t161 * t211) * t266 + t174; m(6) * (t160 * t210 + t161 * t209); ((t160 * t80 + t161 * t82) * (((-t82 + t147) * t164 + t225) * t160 + (-t161 * t190 + t164 * t80 + t240) * t161) + t235 * t116 * t105) * t272 + t160 * ((t160 * t45 + (t21 + t280) * t164) * t160 + (t22 * t164 + (t233 * t73 + t234 * t77) * t161 + (-t46 + (-qJD(4) * t74 + t164 * t77) * t159 + (-qJD(4) * t78 - t164 * t73) * t158) * t160) * t161) - t161 * ((t161 * t46 + (t18 + t278) * t164) * t161 + (t17 * t164 + (-t233 * t74 - t234 * t78) * t160 + (-t45 + (qJD(4) * t73 + t164 * t78) * t159 + (qJD(4) * t77 - t164 * t74) * t158) * t161) * t160) + (t1 * t14 - t27 * t64 - t28 * t63) * t271 + t160 * ((t160 * t47 + (t19 - t279) * t164) * t160 + (t20 * t164 + (-t233 * t67 + t234 * t75) * t161 + (-t48 + (qJD(4) * t68 + t164 * t75) * t159 + (-qJD(4) * t76 + t164 * t67) * t158) * t160) * t161) - t161 * ((t161 * t48 + (t16 - t277) * t164) * t161 + (t15 * t164 + (t233 * t68 - t234 * t76) * t160 + (-t47 + (-qJD(4) * t67 + t164 * t76) * t159 + (qJD(4) * t75 + t164 * t68) * t158) * t161) * t160) + ((-t15 - t17) * t161 + (t16 + t18) * t160) * t244 + ((-t19 - t21) * t161 + (t20 + t22) * t160) * t243; m(6) * ((t160 * t38 + t161 * t37) * t233 + (-t160 * t218 + t161 * t217) * t158); m(6) * ((t160 * t40 + t161 * t39) * t233 + (-t160 * t216 + t161 * t215) * t158); 0; m(6) * ((-t1 + (-t160 * t63 - t161 * t64) * qJD(4)) * t159 + (qJD(4) * t14 - t160 * t209 + t210 * t161) * t158); (-0.1e1 + t235) * t158 * t233 * t271;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
