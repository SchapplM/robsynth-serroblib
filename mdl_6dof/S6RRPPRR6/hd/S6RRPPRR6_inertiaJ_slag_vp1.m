% Calculate joint inertia matrix for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:07
% EndTime: 2019-03-09 09:13:14
% DurationCPUTime: 3.76s
% Computational Cost: add. (7087->400), mult. (10221->595), div. (0->0), fcn. (11897->10), ass. (0->201)
t285 = Icges(3,1) + Icges(4,1);
t283 = Icges(4,4) + Icges(3,5);
t170 = sin(qJ(2));
t284 = (Icges(3,4) - Icges(4,5)) * t170;
t282 = Icges(4,2) + Icges(3,3);
t173 = cos(qJ(2));
t281 = t283 * t173 + (-Icges(3,6) + Icges(4,6)) * t170;
t280 = t285 * t173 - t284;
t171 = sin(qJ(1));
t279 = -t171 / 0.2e1;
t271 = t171 / 0.2e1;
t174 = cos(qJ(1));
t278 = -t174 / 0.2e1;
t277 = t174 / 0.2e1;
t222 = pkin(10) + qJ(5);
t157 = sin(t222);
t215 = cos(t222);
t124 = t170 * t157 + t173 * t215;
t100 = t124 * t171;
t236 = Icges(6,3) * t174;
t237 = Icges(6,5) * t174;
t242 = Icges(6,1) * t100;
t203 = t170 * t215;
t125 = -t173 * t157 + t203;
t99 = t125 * t171;
t249 = Icges(6,6) * t99;
t254 = Icges(6,4) * t99;
t267 = t99 ^ 2;
t229 = t173 * t174;
t101 = t157 * t229 - t174 * t203;
t102 = t124 * t174;
t169 = sin(qJ(6));
t172 = cos(qJ(6));
t85 = -t102 * t169 - t171 * t172;
t86 = t102 * t172 - t169 * t171;
t39 = Icges(7,5) * t86 + Icges(7,6) * t85 + Icges(7,3) * t101;
t40 = Icges(7,4) * t86 + Icges(7,2) * t85 + Icges(7,6) * t101;
t41 = Icges(7,1) * t86 + Icges(7,4) * t85 + Icges(7,5) * t101;
t83 = -t100 * t169 + t172 * t174;
t84 = t100 * t172 + t169 * t174;
t11 = -t39 * t99 + t40 * t83 + t41 * t84;
t248 = Icges(7,6) * t99;
t250 = Icges(7,2) * t83;
t253 = Icges(7,4) * t84;
t178 = (-0.2e1 * t248 + t250 + 0.2e1 * t253) * t83;
t251 = Icges(7,5) * t99;
t255 = Icges(7,1) * t84;
t3 = t11 * t171 - (Icges(7,3) * t267 + (-0.2e1 * t251 + t255) * t84 + t178) * t174;
t57 = Icges(6,5) * t102 - Icges(6,6) * t101 - Icges(6,3) * t171;
t58 = Icges(6,4) * t102 - Icges(6,2) * t101 - Icges(6,6) * t171;
t59 = Icges(6,1) * t102 - Icges(6,4) * t101 - Icges(6,5) * t171;
t263 = -(t100 * t59 + t174 * t57 + t99 * t58) * t171 + (t267 * Icges(6,2) + (t236 + 0.2e1 * t249) * t174 + (0.2e1 * t237 + t242 + 0.2e1 * t254) * t100) * t174 - t3;
t179 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t174;
t180 = t237 + t242 + t254;
t12 = t101 * t39 + t40 * t85 + t41 * t86;
t247 = Icges(7,3) * t99;
t252 = Icges(7,5) * t84;
t181 = Icges(7,6) * t83 - t247 + t252;
t182 = -t248 + t250 + t253;
t183 = Icges(7,4) * t83 - t251 + t255;
t175 = t101 * t181 + t85 * t182 + t86 * t183;
t4 = t12 * t171 - t175 * t174;
t262 = (-t101 * t58 + t102 * t59 - t171 * t57) * t171 - (t102 * t180 - t101 * t179 - t171 * (Icges(6,5) * t100 + t236 + t249)) * t174 + t4;
t276 = t101 / 0.2e1;
t275 = t124 / 0.2e1;
t167 = cos(pkin(10));
t166 = sin(pkin(10));
t233 = t166 * t170;
t192 = t167 * t173 + t233;
t274 = t192 / 0.2e1;
t232 = t166 * t173;
t130 = t167 * t170 - t232;
t273 = -t130 / 0.2e1;
t270 = -t171 * t281 + t174 * t282;
t269 = t171 * t282 + t174 * t281;
t43 = t86 * rSges(7,1) + t85 * rSges(7,2) + t101 * rSges(7,3);
t258 = t102 * pkin(5) + pkin(9) * t101 + t43;
t261 = pkin(5) * t100;
t207 = -rSges(7,1) * t84 - rSges(7,2) * t83;
t42 = -rSges(7,3) * t99 - t207;
t15 = -(-pkin(9) * t99 + t261 + t42) * t171 - t258 * t174;
t268 = -t262 * t171 - t263 * t174;
t164 = t171 ^ 2;
t165 = t174 ^ 2;
t266 = m(4) / 0.2e1;
t265 = m(5) / 0.2e1;
t264 = m(7) / 0.2e1;
t154 = pkin(4) * t167 + pkin(3);
t260 = -pkin(3) + t154;
t51 = Icges(7,3) * t124 + (Icges(7,5) * t172 - Icges(7,6) * t169) * t125;
t53 = Icges(7,5) * t124 + (Icges(7,1) * t172 - Icges(7,4) * t169) * t125;
t257 = t125 * t172 * t53 + t124 * t51;
t54 = rSges(7,3) * t124 + (rSges(7,1) * t172 - rSges(7,2) * t169) * t125;
t256 = pkin(5) * t125 + pkin(9) * t124 + t54;
t52 = Icges(7,6) * t124 + (Icges(7,4) * t172 - Icges(7,2) * t169) * t125;
t246 = t169 * t52;
t245 = t174 * rSges(4,2);
t244 = t174 * rSges(3,3);
t243 = -rSges(5,3) - qJ(4);
t240 = Icges(3,4) * t173;
t238 = Icges(4,5) * t173;
t235 = qJ(3) * t170;
t234 = qJ(4) * t174;
t231 = t170 * t174;
t120 = t130 * t174;
t121 = t192 * t174;
t228 = t121 * rSges(5,1) + t120 * rSges(5,2);
t225 = pkin(2) * t229 + qJ(3) * t231;
t227 = t164 * (pkin(2) * t173 + t235) + t174 * t225;
t142 = pkin(2) * t170 - qJ(3) * t173;
t226 = -rSges(4,1) * t170 + rSges(4,3) * t173 - t142;
t224 = t174 * pkin(1) + t171 * pkin(7);
t223 = t164 + t165;
t221 = pkin(4) * t233;
t13 = t124 * t181 + (-t169 * t182 + t172 * t183) * t125;
t16 = -t51 * t99 + t52 * t83 + t53 * t84;
t220 = t13 / 0.2e1 + t16 / 0.2e1;
t14 = t124 * t39 + (-t169 * t40 + t172 * t41) * t125;
t17 = t101 * t51 + t52 * t85 + t53 * t86;
t219 = t17 / 0.2e1 + t14 / 0.2e1;
t168 = -pkin(8) - qJ(4);
t218 = t154 * t229 + t171 * t168 + t174 * t221;
t217 = rSges(4,1) * t229 + t171 * rSges(4,2) + rSges(4,3) * t231;
t216 = -pkin(3) * t170 - t142;
t152 = pkin(3) * t229;
t214 = -qJ(4) * t171 + t152;
t213 = t265 + m(6) / 0.2e1 + t264;
t212 = t171 * (t171 * t173 * pkin(3) + t234) + t174 * t214 + t227;
t211 = t224 + t225;
t61 = t102 * rSges(6,1) - t101 * rSges(6,2) - rSges(6,3) * t171;
t210 = -rSges(5,1) * t130 + rSges(5,2) * t192 + t216;
t209 = pkin(4) * t232 - t260 * t170 + t216;
t208 = Icges(5,5) * t273 + Icges(5,6) * t274 + t283 * t170 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(3,6) / 0.2e1) * t173;
t206 = rSges(3,1) * t173 - rSges(3,2) * t170;
t118 = t130 * t171;
t119 = t192 * t171;
t205 = t119 * rSges(5,1) + t118 * rSges(5,2);
t79 = rSges(6,1) * t125 - rSges(6,2) * t124;
t191 = t209 - t79;
t44 = t191 * t171;
t45 = t191 * t174;
t204 = t171 * t44 + t174 * t45;
t60 = t100 * rSges(6,1) + t99 * rSges(6,2) + rSges(6,3) * t174;
t34 = -t171 * t60 - t174 * t61;
t200 = -Icges(3,2) * t170 + t240;
t197 = Icges(4,3) * t170 + t238;
t190 = rSges(3,1) * t229 - rSges(3,2) * t231 + t171 * rSges(3,3);
t76 = Icges(6,5) * t125 - Icges(6,6) * t124;
t77 = Icges(6,4) * t125 - Icges(6,2) * t124;
t78 = Icges(6,1) * t125 - Icges(6,4) * t124;
t189 = -t124 * t179 / 0.2e1 + t125 * t180 / 0.2e1 + t100 * t78 / 0.2e1 + t76 * t277 + t99 * t77 / 0.2e1 + t220;
t188 = t77 * t276 - t102 * t78 / 0.2e1 + t76 * t271 + t58 * t275 - t125 * t59 / 0.2e1 - t219;
t156 = t174 * t168;
t187 = t171 * (-t234 - t156 + (t260 * t173 + t221) * t171) + t174 * (-t214 + t218) + t212;
t186 = t209 - t256;
t161 = t174 * pkin(7);
t177 = t156 + t161 + (-pkin(1) + (-pkin(2) - t154) * t173 + (-pkin(4) * t166 - qJ(3)) * t170) * t171;
t37 = t177 - t60;
t184 = t211 + t218;
t38 = t184 + t61;
t185 = m(6) * (t171 * t38 + t174 * t37);
t1 = t11 * t101 + t16 * t124 - (Icges(7,1) * t84 ^ 2 - (-t247 + 0.2e1 * t252) * t99 + t178) * t99;
t2 = t12 * t101 + t17 * t124 - t175 * t99;
t176 = t4 * t276 + (-t13 * t174 + t14 * t171) * t275 + t2 * t271 + t1 * t278 - t99 * t3 / 0.2e1;
t146 = rSges(2,1) * t174 - t171 * rSges(2,2);
t145 = -t171 * rSges(2,1) - rSges(2,2) * t174;
t144 = rSges(3,1) * t170 + rSges(3,2) * t173;
t94 = t226 * t174;
t93 = t226 * t171;
t92 = t190 + t224;
t91 = t244 + t161 + (-pkin(1) - t206) * t171;
t89 = Icges(5,1) * t130 - Icges(5,4) * t192;
t88 = Icges(5,4) * t130 - Icges(5,2) * t192;
t74 = t211 + t217;
t73 = t245 + t161 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t173 + (-rSges(4,3) - qJ(3)) * t170) * t171;
t70 = t174 * t190 + (t206 * t171 - t244) * t171;
t69 = Icges(5,1) * t121 + Icges(5,4) * t120 - Icges(5,5) * t171;
t68 = Icges(5,1) * t119 + Icges(5,4) * t118 + Icges(5,5) * t174;
t67 = Icges(5,4) * t121 + Icges(5,2) * t120 - Icges(5,6) * t171;
t66 = Icges(5,4) * t119 + Icges(5,2) * t118 + Icges(5,6) * t174;
t65 = Icges(5,5) * t121 + Icges(5,6) * t120 - Icges(5,3) * t171;
t64 = Icges(5,5) * t119 + Icges(5,6) * t118 + Icges(5,3) * t174;
t56 = t210 * t174;
t55 = t210 * t171;
t50 = t243 * t171 + t152 + t211 + t228;
t49 = t161 + t243 * t174 + (-t235 - pkin(1) + (-pkin(2) - pkin(3)) * t173) * t171 - t205;
t47 = t174 * t217 + (-t245 + (rSges(4,1) * t173 + rSges(4,3) * t170) * t171) * t171 + t227;
t36 = t256 * t174;
t35 = t256 * t171;
t31 = t186 * t174;
t30 = t186 * t171;
t29 = t171 * t205 + t174 * t228 + t212;
t26 = t184 + t258;
t25 = -t261 - (-rSges(7,3) - pkin(9)) * t99 + t177 + t207;
t24 = -t101 * t54 + t124 * t43;
t23 = -t124 * t42 - t54 * t99;
t20 = t101 * t42 + t43 * t99;
t19 = (-t125 * t246 + t257) * t124;
t18 = t187 - t34;
t6 = -t15 + t187;
t5 = [-t124 * t77 - t192 * t88 + t130 * t89 + Icges(2,3) + (t78 - t246) * t125 + m(7) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t49 ^ 2 + t50 ^ 2) + m(4) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t91 ^ 2 + t92 ^ 2) + m(2) * (t145 ^ 2 + t146 ^ 2) + t257 + ((Icges(3,2) + Icges(4,3)) * t173 + t284) * t173 + (t285 * t170 - t238 + t240) * t170; (-t118 * t88 / 0.2e1 - t119 * t89 / 0.2e1 + t66 * t274 + t68 * t273 + t208 * t174 + (Icges(3,6) * t277 + Icges(4,6) * t278 + t197 * t271 + t200 * t279) * t173 + (t277 * t283 + t280 * t279) * t170 - t189) * t174 + (t120 * t88 / 0.2e1 + t121 * t89 / 0.2e1 - t192 * t67 / 0.2e1 + t130 * t69 / 0.2e1 + t208 * t171 + (Icges(3,6) * t271 + Icges(4,6) * t279 + t197 * t278 + t200 * t277) * t173 + (t271 * t283 + t280 * t277) * t170 - t188) * t171 + m(7) * (t25 * t31 + t26 * t30) + m(6) * (t37 * t45 + t38 * t44) + m(5) * (t49 * t56 + t50 * t55) + m(4) * (t73 * t94 + t74 * t93) + m(3) * (-t171 * t92 - t174 * t91) * t144; m(7) * (t30 ^ 2 + t31 ^ 2 + t6 ^ 2) + m(6) * (t18 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t29 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(4) * (t47 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t223 * t144 ^ 2 + t70 ^ 2) + ((t118 * t66 + t119 * t68 + t174 * t64) * t174 + t270 * t165 + t263) * t174 + ((t120 * t67 + t121 * t69 - t171 * t65) * t171 + t269 * t164 + (-t118 * t67 - t119 * t69 - t120 * t66 - t121 * t68 + (t64 + t270) * t171 + (-t65 + t269) * t174) * t174 + t262) * t171; 0.2e1 * ((t171 * t26 + t174 * t25) * t264 + t185 / 0.2e1 + (t171 * t50 + t174 * t49) * t265 + (t171 * t74 + t174 * t73) * t266) * t170; m(7) * (-t173 * t6 + (t171 * t30 + t174 * t31) * t170) + m(6) * (t204 * t170 - t173 * t18) + m(5) * (-t173 * t29 + (t171 * t55 + t174 * t56) * t170) + m(4) * (-t173 * t47 + (t171 * t93 + t174 * t94) * t170); 0.2e1 * (t266 + t213) * (t223 * t170 ^ 2 + t173 ^ 2); m(7) * (-t171 * t25 + t174 * t26) + m(6) * (-t171 * t37 + t174 * t38) + m(5) * (-t171 * t49 + t174 * t50); m(7) * (-t171 * t31 + t174 * t30) + m(6) * (-t171 * t45 + t174 * t44) + m(5) * (-t171 * t56 + t174 * t55); 0; 0.2e1 * t213 * t223; t189 * t174 + t188 * t171 + m(7) * (t25 * t36 + t26 * t35) + t79 * t185; m(7) * (t15 * t6 + t30 * t35 + t31 * t36) + m(6) * (t34 * t18 + t204 * t79) + t268; m(6) * (t223 * t79 * t170 - t173 * t34) + m(7) * (-t15 * t173 + (t171 * t35 + t174 * t36) * t170); m(7) * (-t36 * t171 + t174 * t35); m(6) * (t223 * t79 ^ 2 + t34 ^ 2) + m(7) * (t15 ^ 2 + t35 ^ 2 + t36 ^ 2) - t268; m(7) * (t23 * t25 + t24 * t26) + t19 - t220 * t99 + t219 * t101; m(7) * (t20 * t6 + t23 * t31 + t24 * t30) + t176; m(7) * (-t20 * t173 + (t171 * t24 + t174 * t23) * t170); m(7) * (-t23 * t171 + t174 * t24); m(7) * (t15 * t20 + t23 * t36 + t24 * t35) - t176; t101 * t2 - t99 * t1 + t124 * (t14 * t101 - t13 * t99 + t19) + m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
