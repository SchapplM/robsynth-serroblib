% Calculate joint inertia matrix for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:32
% EndTime: 2019-03-09 10:03:40
% DurationCPUTime: 3.80s
% Computational Cost: add. (3457->451), mult. (8382->634), div. (0->0), fcn. (9027->6), ass. (0->204)
t279 = rSges(7,1) + pkin(5);
t285 = Icges(4,1) + Icges(3,3);
t188 = sin(qJ(2));
t191 = cos(qJ(2));
t284 = (-Icges(4,4) + Icges(3,5)) * t191 + (Icges(4,5) - Icges(3,6)) * t188;
t190 = cos(qJ(4));
t192 = cos(qJ(1));
t239 = t190 * t192;
t187 = sin(qJ(4));
t189 = sin(qJ(1));
t246 = t187 * t189;
t147 = -t188 * t239 + t246;
t243 = t189 * t190;
t244 = t188 * t192;
t148 = t187 * t244 + t243;
t238 = t191 * t192;
t254 = -rSges(7,3) - qJ(6);
t253 = t147 * rSges(7,2) + t279 * t148 + t238 * t254;
t283 = -t189 / 0.2e1;
t282 = t189 / 0.2e1;
t281 = -t192 / 0.2e1;
t280 = t192 / 0.2e1;
t121 = -Icges(7,5) * t188 + (-Icges(7,1) * t187 + Icges(7,4) * t190) * t191;
t122 = Icges(6,4) * t188 + (-Icges(6,1) * t187 + Icges(6,5) * t190) * t191;
t278 = -t122 - t121;
t112 = Icges(6,6) * t188 + (-Icges(6,5) * t187 + Icges(6,3) * t190) * t191;
t113 = Icges(5,3) * t188 + (-Icges(5,5) * t187 - Icges(5,6) * t190) * t191;
t116 = -Icges(7,6) * t188 + (-Icges(7,4) * t187 + Icges(7,2) * t190) * t191;
t117 = Icges(6,2) * t188 + (-Icges(6,4) * t187 + Icges(6,6) * t190) * t191;
t240 = t190 * t191;
t277 = (t112 + t116) * t240 + (t113 + t117) * t188;
t264 = m(7) / 0.2e1;
t265 = m(6) / 0.2e1;
t227 = t265 + t264;
t276 = 0.2e1 * t227;
t275 = t188 / 0.2e1;
t274 = -t284 * t189 + t285 * t192;
t273 = t285 * t189 + t284 * t192;
t67 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t238;
t73 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t238;
t79 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t238;
t19 = t147 * t73 + t148 * t79 - t238 * t67;
t149 = t187 * t192 + t188 * t243;
t151 = t188 * t246 - t239;
t242 = t189 * t191;
t68 = Icges(7,5) * t151 - Icges(7,6) * t149 - Icges(7,3) * t242;
t74 = Icges(7,4) * t151 - Icges(7,2) * t149 - Icges(7,6) * t242;
t80 = Icges(7,1) * t151 - Icges(7,4) * t149 - Icges(7,5) * t242;
t20 = t147 * t74 + t148 * t80 - t238 * t68;
t69 = Icges(6,5) * t148 + Icges(6,6) * t238 + Icges(6,3) * t147;
t75 = Icges(6,4) * t148 + Icges(6,2) * t238 + Icges(6,6) * t147;
t81 = Icges(6,1) * t148 + Icges(6,4) * t238 + Icges(6,5) * t147;
t21 = t147 * t69 + t148 * t81 + t238 * t75;
t70 = Icges(6,5) * t151 + Icges(6,6) * t242 - Icges(6,3) * t149;
t76 = Icges(6,4) * t151 + Icges(6,2) * t242 - Icges(6,6) * t149;
t82 = Icges(6,1) * t151 + Icges(6,4) * t242 - Icges(6,5) * t149;
t22 = t147 * t70 + t148 * t82 + t238 * t76;
t71 = Icges(5,5) * t148 - Icges(5,6) * t147 + Icges(5,3) * t238;
t77 = Icges(5,4) * t148 - Icges(5,2) * t147 + Icges(5,6) * t238;
t83 = Icges(5,1) * t148 - Icges(5,4) * t147 + Icges(5,5) * t238;
t23 = -t147 * t77 + t148 * t83 + t238 * t71;
t72 = Icges(5,5) * t151 + Icges(5,6) * t149 + Icges(5,3) * t242;
t78 = Icges(5,4) * t151 + Icges(5,2) * t149 + Icges(5,6) * t242;
t84 = Icges(5,1) * t151 + Icges(5,4) * t149 + Icges(5,5) * t242;
t24 = -t147 * t78 + t148 * t84 + t238 * t72;
t111 = -Icges(7,3) * t188 + (-Icges(7,5) * t187 + Icges(7,6) * t190) * t191;
t44 = -t111 * t238 + t147 * t116 + t148 * t121;
t45 = t147 * t112 + t117 * t238 + t148 * t122;
t118 = Icges(5,6) * t188 + (-Icges(5,4) * t187 - Icges(5,2) * t190) * t191;
t123 = Icges(5,5) * t188 + (-Icges(5,1) * t187 - Icges(5,4) * t190) * t191;
t46 = t113 * t238 - t147 * t118 + t148 * t123;
t271 = ((t19 + t21 + t23) * t192 + (t20 + t22 + t24) * t189) * t191 + (t44 + t45 + t46) * t188;
t25 = -t149 * t73 + t151 * t79 - t242 * t67;
t26 = -t149 * t74 + t151 * t80 - t242 * t68;
t27 = -t149 * t69 + t151 * t81 + t242 * t75;
t28 = -t149 * t70 + t151 * t82 + t242 * t76;
t29 = t149 * t77 + t151 * t83 + t242 * t71;
t30 = t149 * t78 + t151 * t84 + t242 * t72;
t47 = -t111 * t242 - t116 * t149 + t121 * t151;
t48 = -t112 * t149 + t117 * t242 + t122 * t151;
t49 = t113 * t242 + t118 * t149 + t123 * t151;
t270 = ((t25 + t27 + t29) * t192 + (t26 + t28 + t30) * t189) * t191 + (t47 + t48 + t49) * t188;
t32 = -t188 * t67 + (-t187 * t79 + t190 * t73) * t191;
t34 = t188 * t75 + (-t187 * t81 + t190 * t69) * t191;
t36 = t188 * t71 + (-t187 * t83 - t190 * t77) * t191;
t269 = t32 + t34 + t36;
t33 = -t188 * t68 + (-t187 * t80 + t190 * t74) * t191;
t35 = t188 * t76 + (-t187 * t82 + t190 * t70) * t191;
t37 = t188 * t72 + (-t187 * t84 - t190 * t78) * t191;
t268 = t33 + t35 + t37;
t184 = t189 ^ 2;
t186 = t192 ^ 2;
t267 = m(4) / 0.2e1;
t266 = m(5) / 0.2e1;
t263 = -pkin(2) - pkin(8);
t259 = m(7) * t191;
t258 = t149 * rSges(7,2);
t257 = t149 * rSges(6,3);
t256 = t192 * rSges(4,1);
t255 = t192 * rSges(3,3);
t252 = t279 * t151 + t254 * t242 - t258;
t251 = Icges(3,4) * t188;
t250 = Icges(3,4) * t191;
t249 = Icges(4,6) * t188;
t248 = Icges(4,6) * t191;
t247 = qJ(3) * t188;
t245 = t187 * t191;
t241 = t190 * t118;
t236 = (-rSges(7,1) * t187 + rSges(7,2) * t190) * t191 - pkin(5) * t245 + t254 * t188;
t232 = pkin(2) * t238 + qJ(3) * t244;
t235 = t184 * (pkin(2) * t191 + t247) + t192 * t232;
t163 = pkin(2) * t188 - qJ(3) * t191;
t233 = rSges(4,2) * t188 + rSges(4,3) * t191 - t163;
t231 = t192 * pkin(1) + t189 * pkin(7);
t230 = t189 * pkin(3) + pkin(8) * t238;
t180 = t192 * pkin(7);
t181 = t192 * pkin(3);
t229 = t180 + t181;
t228 = t184 + t186;
t226 = ((-t123 * t187 - t241) * t191 - t188 * t111 + t278 * t245 + t277) * t188;
t137 = t149 * qJ(5);
t225 = t137 + t229;
t86 = t148 * rSges(6,1) + rSges(6,2) * t238 + t147 * rSges(6,3);
t87 = t148 * rSges(5,1) - t147 * rSges(5,2) + rSges(5,3) * t238;
t224 = -Icges(4,4) * t188 / 0.2e1 + Icges(3,5) * t275 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t191;
t223 = -pkin(1) - t247;
t222 = -pkin(8) * t188 - t163;
t98 = t148 * pkin(4) + t147 * qJ(5);
t221 = t189 * (pkin(8) * t242 - t181) + t192 * t230 + t235;
t220 = t231 + t232;
t134 = t188 * rSges(5,3) + (-rSges(5,1) * t187 - rSges(5,2) * t190) * t191;
t219 = -t134 + t222;
t154 = (-pkin(4) * t187 + qJ(5) * t190) * t191;
t218 = -t154 + t222;
t217 = rSges(3,1) * t191 - rSges(3,2) * t188;
t216 = -t151 * rSges(5,1) - t149 * rSges(5,2);
t110 = t154 * t242;
t99 = t151 * pkin(4) - t137;
t17 = t110 + t236 * t242 + (-t99 - t252) * t188;
t95 = t188 * t98;
t18 = t95 + t253 * t188 + (-t154 - t236) * t238;
t215 = t17 * t192 + t18 * t189;
t38 = t258 + (-pkin(4) - t279) * t151 + ((-t254 + t263) * t191 + t223) * t189 + t225;
t197 = t220 + t230;
t193 = t197 + t98;
t39 = t193 + t253;
t214 = t189 * t39 + t192 * t38;
t196 = t218 - t236;
t55 = t196 * t189;
t56 = t196 * t192;
t213 = t189 * t55 + t192 * t56;
t212 = Icges(3,1) * t191 - t251;
t211 = -Icges(3,2) * t188 + t250;
t208 = -Icges(4,2) * t191 + t249;
t207 = Icges(4,3) * t188 - t248;
t202 = t147 * t192 - t149 * t189;
t133 = t188 * rSges(6,2) + (-rSges(6,1) * t187 + rSges(6,3) * t190) * t191;
t201 = -t133 + t218;
t200 = rSges(3,1) * t238 - rSges(3,2) * t244 + t189 * rSges(3,3);
t199 = t189 * rSges(4,1) - rSges(4,2) * t238 + rSges(4,3) * t244;
t198 = t189 * t99 + t192 * t98 + t221;
t195 = t34 / 0.2e1 + t32 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t36 / 0.2e1;
t194 = t37 / 0.2e1 + t35 / 0.2e1 + t33 / 0.2e1 + t49 / 0.2e1 + t48 / 0.2e1 + t47 / 0.2e1;
t185 = t191 ^ 2;
t183 = t188 ^ 2;
t167 = rSges(2,1) * t192 - t189 * rSges(2,2);
t166 = -t189 * rSges(2,1) - rSges(2,2) * t192;
t165 = rSges(3,1) * t188 + rSges(3,2) * t191;
t103 = t233 * t192;
t102 = t233 * t189;
t101 = t200 + t231;
t100 = t255 + t180 + (-pkin(1) - t217) * t189;
t93 = t99 * t238;
t92 = t199 + t220;
t91 = t256 + t180 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t191 + (-rSges(4,3) - qJ(3)) * t188) * t189;
t90 = rSges(5,3) * t242 - t216;
t89 = t151 * rSges(6,1) + rSges(6,2) * t242 - t257;
t66 = t219 * t192;
t65 = t219 * t189;
t64 = t192 * t200 + (t189 * t217 - t255) * t189;
t62 = t201 * t192;
t61 = t201 * t189;
t60 = -t134 * t238 + t188 * t87;
t59 = t134 * t242 - t188 * t90;
t58 = t197 + t87;
t57 = ((-rSges(5,3) + t263) * t191 + t223) * t189 + t216 + t229;
t51 = t192 * t199 + (-t256 + (-rSges(4,2) * t191 + rSges(4,3) * t188) * t189) * t189 + t235;
t50 = (-t189 * t87 + t192 * t90) * t191;
t43 = t193 + t86;
t42 = t257 + (-rSges(6,1) - pkin(4)) * t151 + ((-rSges(6,2) + t263) * t191 + t223) * t189 + t225;
t41 = t188 * t86 + t95 + (-t133 - t154) * t238;
t40 = t133 * t242 + t110 + (-t89 - t99) * t188;
t31 = t189 * t90 + t192 * t87 + t221;
t16 = t93 + (t192 * t89 + (-t86 - t98) * t189) * t191;
t15 = t189 * t89 + t192 * t86 + t198;
t14 = t93 + (t252 * t192 + (-t98 - t253) * t189) * t191;
t13 = t189 * t252 + t192 * t253 + t198;
t12 = t29 * t189 - t192 * t30;
t11 = t27 * t189 - t192 * t28;
t10 = t25 * t189 - t192 * t26;
t9 = t23 * t189 - t192 * t24;
t8 = t21 * t189 - t192 * t22;
t7 = t19 * t189 - t192 * t20;
t1 = [Icges(2,3) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(7) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t91 ^ 2 + t92 ^ 2) + m(3) * (t100 ^ 2 + t101 ^ 2) + m(2) * (t166 ^ 2 + t167 ^ 2) + (-t241 + t249 + t251 + (-t123 + t278) * t187 + (Icges(4,3) + Icges(3,2)) * t191) * t191 + (-t111 + t248 + t250 + (Icges(3,1) + Icges(4,2)) * t188) * t188 + t277; (t224 * t192 + (Icges(4,5) * t281 + Icges(3,6) * t280 + t207 * t282 + t211 * t283) * t191 + (Icges(4,4) * t281 + Icges(3,5) * t280 + t208 * t282 + t212 * t283) * t188 - t194) * t192 + (t224 * t189 + (Icges(4,5) * t283 + Icges(3,6) * t282 + t207 * t281 + t211 * t280) * t191 + (Icges(4,4) * t283 + Icges(3,5) * t282 + t208 * t281 + t212 * t280) * t188 + t195) * t189 + m(4) * (t102 * t92 + t103 * t91) + m(5) * (t57 * t66 + t58 * t65) + m(6) * (t42 * t62 + t43 * t61) + m(7) * (t38 * t56 + t39 * t55) + m(3) * (-t100 * t192 - t101 * t189) * t165; m(7) * (t13 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t15 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t31 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2 + t51 ^ 2) + m(3) * (t165 ^ 2 * t228 + t64 ^ 2) + (t274 * t186 - t10 - t11 - t12) * t192 + (t7 + t8 + t9 + t273 * t184 + (t274 * t189 + t273 * t192) * t192) * t189; 0.2e1 * ((t189 * t43 + t192 * t42) * t265 + t214 * t264 + (t189 * t58 + t192 * t57) * t266 + (t189 * t92 + t192 * t91) * t267) * t188; m(7) * (-t191 * t13 + t188 * t213) + m(6) * (-t191 * t15 + (t189 * t61 + t192 * t62) * t188) + m(5) * (-t191 * t31 + (t189 * t65 + t192 * t66) * t188) + m(4) * (-t191 * t51 + (t102 * t189 + t103 * t192) * t188); 0.2e1 * (t267 + t266 + t227) * (t183 * t228 + t185); m(6) * (t40 * t42 + t41 * t43) + m(7) * (t17 * t38 + t18 * t39) + m(5) * (t57 * t59 + t58 * t60) + (t189 * t194 + t192 * t195) * t191 + t226; m(7) * (t13 * t14 + t17 * t56 + t18 * t55) + m(6) * (t15 * t16 + t40 * t62 + t41 * t61) + m(5) * (t31 * t50 + t59 * t66 + t60 * t65) + ((t9 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t192 + (t12 / 0.2e1 + t11 / 0.2e1 + t10 / 0.2e1) * t189) * t191 + (t269 * t189 - t268 * t192) * t275 + t271 * t282 + t270 * t281; m(5) * (-t50 * t191 + (t189 * t60 + t192 * t59) * t188) + m(6) * (-t16 * t191 + (t189 * t41 + t192 * t40) * t188) + m(7) * (-t14 * t191 + t188 * t215); t226 * t188 + m(7) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t16 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t50 ^ 2 + t59 ^ 2 + t60 ^ 2) + (t271 * t192 + t270 * t189 + (t268 * t189 + t269 * t192) * t188) * t191; m(6) * (t147 * t42 - t149 * t43) + m(7) * (t147 * t38 - t149 * t39); m(7) * (t13 * t240 + t147 * t56 - t149 * t55) + m(6) * (t147 * t62 - t149 * t61 + t15 * t240); (-t185 * t190 + t188 * t202) * t276; m(7) * (t14 * t240 + t147 * t17 - t149 * t18) + m(6) * (t147 * t40 - t149 * t41 + t16 * t240); (t185 * t190 ^ 2 + t147 ^ 2 + t149 ^ 2) * t276; -t214 * t259; m(7) * (-t188 * t13 - t191 * t213); (0.1e1 - t228) * t188 * t259; m(7) * (-t188 * t14 - t191 * t215); (-t188 * t190 - t202) * t259; m(7) * (t185 * t228 + t183);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
