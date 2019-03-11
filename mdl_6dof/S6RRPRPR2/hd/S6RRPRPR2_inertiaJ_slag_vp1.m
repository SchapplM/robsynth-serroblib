% Calculate joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:05
% EndTime: 2019-03-09 10:12:17
% DurationCPUTime: 5.05s
% Computational Cost: add. (6411->361), mult. (5827->529), div. (0->0), fcn. (5916->10), ass. (0->177)
t298 = Icges(5,4) + Icges(6,6);
t297 = Icges(5,1) + Icges(6,2);
t296 = -Icges(5,2) - Icges(6,3);
t175 = qJ(2) + pkin(10);
t166 = qJ(4) + t175;
t162 = cos(t166);
t295 = t298 * t162;
t161 = sin(t166);
t294 = t298 * t161;
t293 = Icges(6,4) - Icges(5,5);
t292 = Icges(6,5) - Icges(5,6);
t291 = t296 * t161 + t295;
t290 = t297 * t162 - t294;
t181 = sin(qJ(1));
t184 = cos(qJ(1));
t289 = t181 * t291 + t184 * t292;
t288 = -t181 * t292 + t184 * t291;
t287 = t290 * t181 + t184 * t293;
t286 = -t181 * t293 + t290 * t184;
t285 = Icges(6,1) + Icges(5,3);
t284 = t292 * t161 - t162 * t293;
t283 = t296 * t162 - t294;
t282 = t297 * t161 + t295;
t281 = t285 * t181 + t284 * t184;
t280 = -t284 * t181 + t285 * t184;
t279 = Icges(3,3) + Icges(4,3);
t164 = sin(t175);
t165 = cos(t175);
t180 = sin(qJ(2));
t183 = cos(qJ(2));
t278 = Icges(3,5) * t183 + Icges(4,5) * t165 - Icges(3,6) * t180 - Icges(4,6) * t164;
t277 = t289 * t161 - t287 * t162;
t276 = t288 * t161 - t286 * t162;
t176 = t181 ^ 2;
t275 = t181 * pkin(7);
t242 = t162 * t184;
t182 = cos(qJ(6));
t238 = t184 * t182;
t179 = sin(qJ(6));
t241 = t181 * t179;
t121 = t161 * t238 - t241;
t239 = t184 * t179;
t240 = t181 * t182;
t122 = t161 * t239 + t240;
t64 = t122 * rSges(7,1) + t121 * rSges(7,2) + rSges(7,3) * t242;
t274 = t181 * pkin(5) + pkin(9) * t242 + t64;
t273 = -t161 * t293 - t292 * t162;
t272 = t161 * t283 + t282 * t162;
t220 = rSges(4,1) * t165 - rSges(4,2) * t164;
t81 = t161 * rSges(7,3) + (-rSges(7,1) * t179 - rSges(7,2) * t182) * t162;
t271 = -pkin(9) * t161 - t81;
t270 = -t278 * t181 + t279 * t184;
t269 = t279 * t181 + t278 * t184;
t177 = t184 ^ 2;
t123 = t161 * t240 + t239;
t124 = t161 * t241 - t238;
t243 = t162 * t181;
t58 = Icges(7,5) * t122 + Icges(7,6) * t121 + Icges(7,3) * t242;
t60 = Icges(7,4) * t122 + Icges(7,2) * t121 + Icges(7,6) * t242;
t62 = Icges(7,1) * t122 + Icges(7,4) * t121 + Icges(7,5) * t242;
t19 = t123 * t60 + t124 * t62 + t58 * t243;
t59 = Icges(7,5) * t124 + Icges(7,6) * t123 + Icges(7,3) * t243;
t61 = Icges(7,4) * t124 + Icges(7,2) * t123 + Icges(7,6) * t243;
t63 = Icges(7,1) * t124 + Icges(7,4) * t123 + Icges(7,5) * t243;
t20 = t123 * t61 + t124 * t63 + t59 * t243;
t9 = t19 * t181 - t20 * t184;
t268 = -t9 + t280 * t177 + (t276 * t181 + (-t277 + t281) * t184) * t181;
t267 = m(6) / 0.2e1;
t266 = m(7) / 0.2e1;
t265 = t181 / 0.2e1;
t264 = -t184 / 0.2e1;
t263 = pkin(2) * t180;
t178 = -qJ(3) - pkin(7);
t163 = t183 * pkin(2) + pkin(1);
t244 = t161 * t184;
t193 = rSges(5,1) * t242 - rSges(5,2) * t244 + t181 * rSges(5,3);
t219 = rSges(5,1) * t162 - rSges(5,2) * t161;
t53 = t181 * (-t184 * rSges(5,3) + t219 * t181) + t184 * t193;
t157 = t184 * t163;
t173 = t184 * pkin(7);
t261 = t181 * (t173 + (-pkin(1) + t163) * t181) + t184 * (-t184 * pkin(1) + t157 - t275);
t260 = rSges(3,1) * t183;
t258 = rSges(3,2) * t180;
t256 = t184 * rSges(3,3);
t25 = t161 * t58 + (-t179 * t62 - t182 * t60) * t162;
t255 = t25 * t181;
t26 = t161 * t59 + (-t179 * t63 - t182 * t61) * t162;
t254 = t26 * t184;
t253 = Icges(3,4) * t180;
t252 = Icges(3,4) * t183;
t251 = Icges(4,4) * t164;
t250 = Icges(4,4) * t165;
t245 = qJ(5) * t161;
t235 = pkin(4) * t242 + qJ(5) * t244;
t237 = t176 * (pkin(4) * t162 + t245) + t184 * t235;
t134 = t161 * pkin(4) - t162 * qJ(5);
t135 = -t161 * rSges(6,2) - t162 * rSges(6,3);
t236 = -t134 - t135;
t234 = t181 * rSges(3,3) + t184 * t260;
t231 = t176 + t177;
t17 = t121 * t60 + t122 * t62 + t58 * t242;
t18 = t121 * t61 + t122 * t63 + t59 * t242;
t8 = t17 * t181 - t18 * t184;
t230 = (t8 + t281 * t176 + ((-t276 + t280) * t181 + t277 * t184) * t184) * t181;
t229 = t267 + t266;
t228 = -t164 * rSges(4,1) - t165 * rSges(4,2) - t263;
t142 = pkin(3) * t165 + t163;
t137 = t184 * t142;
t227 = t184 * (t137 - t157) + t261 + (t142 - t163) * t176;
t174 = -pkin(8) + t178;
t226 = -t181 * t174 + t137;
t192 = t181 * rSges(6,1) - rSges(6,2) * t242 + rSges(6,3) * t244;
t36 = t181 * (-t184 * rSges(6,1) + (-rSges(6,2) * t162 + rSges(6,3) * t161) * t181) + t184 * t192 + t237;
t78 = Icges(7,3) * t161 + (-Icges(7,5) * t179 - Icges(7,6) * t182) * t162;
t79 = Icges(7,6) * t161 + (-Icges(7,4) * t179 - Icges(7,2) * t182) * t162;
t80 = Icges(7,5) * t161 + (-Icges(7,1) * t179 - Icges(7,4) * t182) * t162;
t29 = t121 * t79 + t122 * t80 + t78 * t242;
t3 = t29 * t161 + (t17 * t184 + t18 * t181) * t162;
t30 = t123 * t79 + t124 * t80 + t78 * t243;
t4 = t30 * t161 + (t181 * t20 + t184 * t19) * t162;
t225 = t4 * t264 + t3 * t265 + t161 * (-t254 + t255) / 0.2e1 + t9 * t243 / 0.2e1 + t8 * t242 / 0.2e1;
t224 = -t134 + t271;
t222 = -pkin(3) * t164 - t263;
t221 = -t258 + t260;
t218 = -t124 * rSges(7,1) - t123 * rSges(7,2);
t213 = -t179 * t80 - t182 * t79;
t212 = Icges(3,1) * t183 - t253;
t211 = Icges(4,1) * t165 - t251;
t209 = -Icges(3,2) * t180 + t252;
t208 = -Icges(4,2) * t164 + t250;
t194 = t181 * rSges(4,3) + t184 * t220;
t191 = -t134 + t222;
t136 = t161 * rSges(5,1) + t162 * rSges(5,2);
t190 = -t136 + t222;
t189 = t226 + t235;
t65 = rSges(7,3) * t243 - t218;
t21 = t237 + t274 * t184 + (-t184 * pkin(5) + pkin(9) * t243 + t65) * t181;
t188 = -t135 + t191;
t187 = t268 * t184 + t230;
t186 = t191 + t271;
t185 = t255 / 0.2e1 - t254 / 0.2e1 + (t286 * t161 + t288 * t162 + t273 * t181 + t272 * t184 + t29) * t265 + (t287 * t161 + t289 * t162 + t272 * t181 - t273 * t184 + t30) * t264;
t155 = t184 * rSges(2,1) - t181 * rSges(2,2);
t154 = -t181 * rSges(2,1) - t184 * rSges(2,2);
t153 = t180 * rSges(3,1) + t183 * rSges(3,2);
t106 = t228 * t184;
t105 = t228 * t181;
t89 = t275 + (pkin(1) - t258) * t184 + t234;
t88 = t256 + t173 + (-pkin(1) - t221) * t181;
t77 = t190 * t184;
t76 = t190 * t181;
t75 = t161 * t78;
t74 = -t181 * t178 + t157 + t194;
t73 = (rSges(4,3) - t178) * t184 + (-t163 - t220) * t181;
t72 = t236 * t184;
t71 = t236 * t181;
t68 = t184 * (-t184 * t258 + t234) + (t221 * t181 - t256) * t181;
t67 = t193 + t226;
t66 = (rSges(5,3) - t174) * t184 + (-t142 - t219) * t181;
t55 = t188 * t184;
t54 = t188 * t181;
t52 = t224 * t184;
t51 = t224 * t181;
t50 = t189 + t192;
t49 = (rSges(6,1) - t174) * t184 + (-t142 + (rSges(6,2) - pkin(4)) * t162 + (-rSges(6,3) - qJ(5)) * t161) * t181;
t40 = t186 * t184;
t39 = t186 * t181;
t38 = t161 * t64 - t81 * t242;
t37 = -t161 * t65 + t81 * t243;
t35 = t184 * t194 + (-t184 * rSges(4,3) + t220 * t181) * t181 + t261;
t34 = t189 + t274;
t33 = (pkin(5) - t174) * t184 + (-t245 - t142 + (-rSges(7,3) - pkin(4) - pkin(9)) * t162) * t181 + t218;
t32 = (t213 * t162 + t75) * t161;
t31 = (-t181 * t64 + t184 * t65) * t162;
t22 = t227 + t53;
t12 = t36 + t227;
t11 = t21 + t227;
t1 = [t165 * (Icges(4,2) * t165 + t251) + t164 * (Icges(4,1) * t164 + t250) + t183 * (Icges(3,2) * t183 + t253) + t180 * (Icges(3,1) * t180 + t252) + Icges(2,3) + t75 + t282 * t161 + (t213 - t283) * t162 + m(7) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2) + m(5) * (t66 ^ 2 + t67 ^ 2) + m(4) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t88 ^ 2 + t89 ^ 2) + m(2) * (t154 ^ 2 + t155 ^ 2); m(3) * (-t181 * t89 - t184 * t88) * t153 + t185 + m(4) * (t105 * t74 + t106 * t73) + m(5) * (t66 * t77 + t67 * t76) + m(6) * (t49 * t55 + t50 * t54) + m(7) * (t33 * t40 + t34 * t39) + (t165 * (Icges(4,6) * t181 + t208 * t184) + t164 * (Icges(4,5) * t181 + t211 * t184) + t183 * (Icges(3,6) * t181 + t209 * t184) + t180 * (Icges(3,5) * t181 + t212 * t184)) * t265 + (t165 * (-Icges(4,6) * t184 + t208 * t181) + t164 * (-Icges(4,5) * t184 + t211 * t181) + t183 * (-Icges(3,6) * t184 + t209 * t181) + t180 * (-Icges(3,5) * t184 + t212 * t181)) * t264 + (Icges(3,5) * t180 + Icges(4,5) * t164 + Icges(3,6) * t183 + Icges(4,6) * t165) * (t176 / 0.2e1 + t177 / 0.2e1); m(7) * (t11 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t12 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t22 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t105 ^ 2 + t106 ^ 2 + t35 ^ 2) + m(3) * (t231 * t153 ^ 2 + t68 ^ 2) + t230 + t269 * t181 * t176 + (t270 * t177 + (t270 * t181 + t269 * t184) * t181 + t268) * t184; m(7) * (t181 * t33 - t184 * t34) + m(6) * (t181 * t49 - t184 * t50) + m(5) * (t181 * t66 - t184 * t67) + m(4) * (t181 * t73 - t184 * t74); m(7) * (t181 * t40 - t184 * t39) + m(6) * (t181 * t55 - t184 * t54) + m(5) * (t181 * t77 - t184 * t76) + m(4) * (-t184 * t105 + t181 * t106); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t229) * t231; t185 + m(7) * (t33 * t52 + t34 * t51) + m(6) * (t49 * t72 + t50 * t71) + m(5) * (-t181 * t67 - t184 * t66) * t136; m(7) * (t11 * t21 + t39 * t51 + t40 * t52) + m(6) * (t12 * t36 + t54 * t71 + t55 * t72) + m(5) * (t53 * t22 + (-t181 * t76 - t184 * t77) * t136) + t187; m(6) * (t72 * t181 - t71 * t184) + m(7) * (t52 * t181 - t51 * t184); m(7) * (t21 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t36 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t231 * t136 ^ 2 + t53 ^ 2) + t187; 0.2e1 * ((t181 * t34 + t184 * t33) * t266 + (t181 * t50 + t184 * t49) * t267) * t161; m(7) * (-t162 * t11 + (t181 * t39 + t184 * t40) * t161) + m(6) * (-t162 * t12 + (t181 * t54 + t184 * t55) * t161); 0; m(7) * (-t162 * t21 + (t181 * t51 + t184 * t52) * t161) + m(6) * (-t162 * t36 + (t181 * t71 + t184 * t72) * t161); 0.2e1 * t229 * (t161 ^ 2 * t231 + t162 ^ 2); m(7) * (t33 * t37 + t34 * t38) + t32 + ((t29 / 0.2e1 + t25 / 0.2e1) * t184 + (t26 / 0.2e1 + t30 / 0.2e1) * t181) * t162; m(7) * (t11 * t31 + t37 * t40 + t38 * t39) + t225; m(7) * (t37 * t181 - t38 * t184); m(7) * (t21 * t31 + t37 * t52 + t38 * t51) + t225; m(7) * (-t31 * t162 + (t181 * t38 + t184 * t37) * t161); t161 * t32 + m(7) * (t31 ^ 2 + t37 ^ 2 + t38 ^ 2) + (t184 * t3 + t181 * t4 + t161 * (t181 * t26 + t184 * t25)) * t162;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
