% Calculate joint inertia matrix for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:50
% EndTime: 2019-03-09 09:49:56
% DurationCPUTime: 3.64s
% Computational Cost: add. (6012->442), mult. (7934->634), div. (0->0), fcn. (8560->8), ass. (0->201)
t184 = qJ(2) + pkin(9);
t180 = cos(t184);
t190 = sin(qJ(1));
t191 = cos(qJ(4));
t236 = t190 * t191;
t188 = sin(qJ(4));
t193 = cos(qJ(1));
t238 = t188 * t193;
t153 = t180 * t238 - t236;
t235 = t191 * t193;
t237 = t190 * t188;
t154 = t180 * t235 + t237;
t179 = sin(t184);
t240 = t179 * t193;
t279 = rSges(7,3) + qJ(6);
t280 = rSges(7,1) + pkin(5);
t249 = t153 * rSges(7,2) + t280 * t154 - t240 * t279;
t282 = Icges(3,3) + Icges(4,3);
t189 = sin(qJ(2));
t192 = cos(qJ(2));
t281 = Icges(3,5) * t192 + Icges(4,5) * t180 - Icges(3,6) * t189 - Icges(4,6) * t179;
t115 = -Icges(5,3) * t180 + (Icges(5,5) * t191 - Icges(5,6) * t188) * t179;
t117 = -Icges(6,2) * t180 + (Icges(6,4) * t191 + Icges(6,6) * t188) * t179;
t278 = t115 + t117;
t226 = m(6) / 0.2e1 + m(7) / 0.2e1;
t277 = 0.2e1 * t226;
t276 = t179 / 0.2e1;
t275 = t180 / 0.2e1;
t274 = t189 / 0.2e1;
t273 = t192 / 0.2e1;
t272 = -t281 * t190 + t282 * t193;
t271 = t190 * t282 + t281 * t193;
t151 = t180 * t237 + t235;
t152 = t180 * t236 - t238;
t242 = t179 * t190;
t65 = Icges(7,5) * t152 + Icges(7,6) * t151 - Icges(7,3) * t242;
t71 = Icges(7,4) * t152 + Icges(7,2) * t151 - Icges(7,6) * t242;
t77 = Icges(7,1) * t152 + Icges(7,4) * t151 - Icges(7,5) * t242;
t20 = t151 * t71 + t152 * t77 - t242 * t65;
t66 = Icges(7,5) * t154 + Icges(7,6) * t153 - Icges(7,3) * t240;
t72 = Icges(7,4) * t154 + Icges(7,2) * t153 - Icges(7,6) * t240;
t78 = Icges(7,1) * t154 + Icges(7,4) * t153 - Icges(7,5) * t240;
t21 = t151 * t72 + t152 * t78 - t242 * t66;
t67 = Icges(6,5) * t152 + Icges(6,6) * t242 + Icges(6,3) * t151;
t73 = Icges(6,4) * t152 + Icges(6,2) * t242 + Icges(6,6) * t151;
t79 = Icges(6,1) * t152 + Icges(6,4) * t242 + Icges(6,5) * t151;
t22 = t151 * t67 + t152 * t79 + t242 * t73;
t68 = Icges(6,5) * t154 + Icges(6,6) * t240 + Icges(6,3) * t153;
t74 = Icges(6,4) * t154 + Icges(6,2) * t240 + Icges(6,6) * t153;
t80 = Icges(6,1) * t154 + Icges(6,4) * t240 + Icges(6,5) * t153;
t23 = t151 * t68 + t152 * t80 + t242 * t74;
t69 = Icges(5,5) * t152 - Icges(5,6) * t151 + Icges(5,3) * t242;
t75 = Icges(5,4) * t152 - Icges(5,2) * t151 + Icges(5,6) * t242;
t81 = Icges(5,1) * t152 - Icges(5,4) * t151 + Icges(5,5) * t242;
t24 = -t151 * t75 + t152 * t81 + t242 * t69;
t70 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t240;
t76 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t240;
t82 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t240;
t25 = -t151 * t76 + t152 * t82 + t242 * t70;
t113 = Icges(7,3) * t180 + (Icges(7,5) * t191 + Icges(7,6) * t188) * t179;
t116 = Icges(7,6) * t180 + (Icges(7,4) * t191 + Icges(7,2) * t188) * t179;
t119 = Icges(7,5) * t180 + (Icges(7,1) * t191 + Icges(7,4) * t188) * t179;
t42 = -t113 * t242 + t116 * t151 + t119 * t152;
t114 = -Icges(6,6) * t180 + (Icges(6,5) * t191 + Icges(6,3) * t188) * t179;
t120 = -Icges(6,4) * t180 + (Icges(6,1) * t191 + Icges(6,5) * t188) * t179;
t43 = t114 * t151 + t117 * t242 + t120 * t152;
t118 = -Icges(5,6) * t180 + (Icges(5,4) * t191 - Icges(5,2) * t188) * t179;
t121 = -Icges(5,5) * t180 + (Icges(5,1) * t191 - Icges(5,4) * t188) * t179;
t44 = t115 * t242 - t118 * t151 + t121 * t152;
t269 = (-t42 - t43 - t44) * t180 + ((t21 + t23 + t25) * t193 + (t20 + t22 + t24) * t190) * t179;
t26 = t153 * t71 + t154 * t77 - t240 * t65;
t27 = t153 * t72 + t154 * t78 - t240 * t66;
t28 = t153 * t67 + t154 * t79 + t240 * t73;
t29 = t153 * t68 + t154 * t80 + t240 * t74;
t30 = -t153 * t75 + t154 * t81 + t240 * t69;
t31 = -t153 * t76 + t154 * t82 + t240 * t70;
t45 = -t113 * t240 + t153 * t116 + t154 * t119;
t46 = t153 * t114 + t117 * t240 + t154 * t120;
t47 = t115 * t240 - t153 * t118 + t154 * t121;
t268 = (-t45 - t46 - t47) * t180 + ((t27 + t29 + t31) * t193 + (t26 + t28 + t30) * t190) * t179;
t32 = t180 * t65 + (t188 * t71 + t191 * t77) * t179;
t34 = -t180 * t73 + (t188 * t67 + t191 * t79) * t179;
t36 = -t180 * t69 + (-t188 * t75 + t191 * t81) * t179;
t267 = -t32 - t34 - t36;
t33 = t180 * t66 + (t188 * t72 + t191 * t78) * t179;
t35 = -t180 * t74 + (t188 * t68 + t191 * t80) * t179;
t37 = -t180 * t70 + (-t188 * t76 + t191 * t82) * t179;
t266 = t33 + t35 + t37;
t241 = t179 * t191;
t243 = t179 * t188;
t265 = t180 * t113 + (t114 + t116) * t243 + (t119 + t120 + t121) * t241;
t264 = t180 ^ 2;
t185 = t190 ^ 2;
t186 = t193 ^ 2;
t263 = -t180 / 0.2e1;
t260 = m(7) * t179;
t259 = pkin(2) * t189;
t258 = pkin(3) * t180;
t87 = t154 * rSges(6,1) + rSges(6,2) * t240 + t153 * rSges(6,3);
t97 = t154 * pkin(4) + t153 * qJ(5);
t257 = -t87 - t97;
t256 = rSges(3,1) * t192;
t255 = rSges(3,2) * t189;
t254 = t151 * rSges(7,2);
t253 = t151 * rSges(6,3);
t252 = t193 * rSges(3,3);
t250 = t152 * t280 - t279 * t242 + t254;
t142 = (pkin(4) * t191 + qJ(5) * t188) * t179;
t141 = t151 * qJ(5);
t96 = t152 * pkin(4) + t141;
t248 = t142 * t242 + t180 * t96;
t247 = Icges(3,4) * t189;
t246 = Icges(3,4) * t192;
t245 = Icges(4,4) * t179;
t244 = Icges(4,4) * t180;
t239 = t180 * t193;
t187 = -qJ(3) - pkin(7);
t234 = t193 * t187;
t177 = pkin(2) * t192 + pkin(1);
t173 = t193 * t177;
t183 = t193 * pkin(7);
t232 = t190 * (t234 + t183 + (-pkin(1) + t177) * t190) + t193 * (-t193 * pkin(1) + t173 + (-pkin(7) - t187) * t190);
t231 = (rSges(7,1) * t191 + rSges(7,2) * t188) * t179 + pkin(5) * t241 + t279 * t180;
t229 = pkin(3) * t239 + pkin(8) * t240;
t228 = t190 * rSges(3,3) + t193 * t256;
t227 = t185 + t186;
t225 = t118 * t243 + t180 * t278 - t265;
t224 = -t97 - t249;
t88 = t154 * rSges(5,1) - t153 * rSges(5,2) + rSges(5,3) * t240;
t223 = Icges(3,5) * t274 + Icges(4,5) * t276 + Icges(3,6) * t273 + Icges(4,6) * t275;
t222 = -rSges(4,1) * t179 - rSges(4,2) * t180 - t259;
t221 = -pkin(3) * t179 + pkin(8) * t180 - t259;
t220 = -t177 - t258;
t219 = -t141 - t234;
t218 = -t190 * t187 + t173;
t217 = t185 * (pkin(8) * t179 + t258) + t193 * t229 + t232;
t124 = -t180 * rSges(5,3) + (rSges(5,1) * t191 - rSges(5,2) * t188) * t179;
t216 = -t124 + t221;
t215 = -t142 + t221;
t214 = -t255 + t256;
t213 = rSges(4,1) * t180 - rSges(4,2) * t179;
t212 = -t152 * rSges(5,1) + t151 * rSges(5,2);
t211 = Icges(3,1) * t192 - t247;
t210 = Icges(4,1) * t180 - t245;
t209 = -Icges(3,2) * t189 + t246;
t208 = -Icges(4,2) * t179 + t244;
t123 = -t180 * rSges(6,2) + (rSges(6,1) * t191 + rSges(6,3) * t188) * t179;
t201 = -t123 + t215;
t200 = rSges(4,1) * t239 - rSges(4,2) * t240 + t190 * rSges(4,3);
t199 = t218 + t229;
t198 = t190 * t96 + t193 * t97 + t217;
t197 = t215 - t231;
t196 = t34 / 0.2e1 + t44 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1 + t32 / 0.2e1 + t36 / 0.2e1;
t195 = t35 / 0.2e1 + t47 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1 + t33 / 0.2e1 + t37 / 0.2e1;
t194 = t199 + t97;
t178 = t179 ^ 2;
t167 = rSges(2,1) * t193 - t190 * rSges(2,2);
t166 = -t190 * rSges(2,1) - rSges(2,2) * t193;
t165 = rSges(3,1) * t189 + rSges(3,2) * t192;
t126 = t222 * t193;
t125 = t222 * t190;
t111 = t190 * pkin(7) + (pkin(1) - t255) * t193 + t228;
t110 = t252 + t183 + (-pkin(1) - t214) * t190;
t99 = t200 + t218;
t98 = (rSges(4,3) - t187) * t193 + (-t177 - t213) * t190;
t90 = t193 * (-t193 * t255 + t228) + (t214 * t190 - t252) * t190;
t89 = t96 * t240;
t85 = rSges(5,3) * t242 - t212;
t84 = t152 * rSges(6,1) + rSges(6,2) * t242 + t253;
t64 = t216 * t193;
t63 = t216 * t190;
t62 = t201 * t193;
t61 = t201 * t190;
t60 = t199 + t88;
t59 = -t234 + ((-rSges(5,3) - pkin(8)) * t179 + t220) * t190 + t212;
t58 = -t124 * t240 - t180 * t88;
t57 = t124 * t242 + t180 * t85;
t56 = t197 * t193;
t55 = t197 * t190;
t51 = t193 * t200 + (-t193 * rSges(4,3) + t190 * t213) * t190 + t232;
t50 = (-t190 * t88 + t193 * t85) * t179;
t49 = t194 + t87;
t48 = -t253 + (-rSges(6,1) - pkin(4)) * t152 + ((-rSges(6,2) - pkin(8)) * t179 + t220) * t190 + t219;
t41 = t194 + t249;
t40 = -t254 + (-pkin(4) - t280) * t152 + ((-pkin(8) + t279) * t179 + t220) * t190 + t219;
t39 = t257 * t180 + (-t123 - t142) * t240;
t38 = t123 * t242 + t180 * t84 + t248;
t19 = t190 * t85 + t193 * t88 + t217;
t18 = t89 + (t190 * t257 + t193 * t84) * t179;
t17 = t224 * t180 + (-t142 - t231) * t240;
t16 = t180 * t250 + t231 * t242 + t248;
t15 = t89 + (t190 * t224 + t193 * t250) * t179;
t14 = t190 * t84 + t193 * t87 + t198;
t13 = t190 * t250 + t193 * t249 + t198;
t12 = t31 * t190 - t193 * t30;
t11 = t29 * t190 - t193 * t28;
t10 = t27 * t190 - t193 * t26;
t9 = t25 * t190 - t193 * t24;
t8 = t23 * t190 - t193 * t22;
t7 = t21 * t190 - t193 * t20;
t1 = [t192 * (Icges(3,2) * t192 + t247) + t189 * (Icges(3,1) * t189 + t246) + Icges(2,3) + (Icges(4,1) * t179 - t118 * t188 + t244) * t179 + (Icges(4,2) * t180 + t245 - t278) * t180 + m(6) * (t48 ^ 2 + t49 ^ 2) + m(7) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t110 ^ 2 + t111 ^ 2) + m(2) * (t166 ^ 2 + t167 ^ 2) + t265; ((-Icges(4,6) * t193 + t208 * t190) * t263 - t179 * (-Icges(4,5) * t193 + t210 * t190) / 0.2e1 - t192 * (-Icges(3,6) * t193 + t209 * t190) / 0.2e1 - t189 * (-Icges(3,5) * t193 + t211 * t190) / 0.2e1 + t223 * t193 - t196) * t193 + ((Icges(4,6) * t190 + t208 * t193) * t275 + (Icges(4,5) * t190 + t210 * t193) * t276 + (Icges(3,6) * t190 + t209 * t193) * t273 + (Icges(3,5) * t190 + t211 * t193) * t274 + t223 * t190 + t195) * t190 + m(5) * (t59 * t64 + t60 * t63) + m(6) * (t48 * t62 + t49 * t61) + m(7) * (t40 * t56 + t41 * t55) + m(4) * (t125 * t99 + t126 * t98) + m(3) * (-t110 * t193 - t111 * t190) * t165; m(6) * (t14 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(7) * (t13 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t19 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t125 ^ 2 + t126 ^ 2 + t51 ^ 2) + m(3) * (t165 ^ 2 * t227 + t90 ^ 2) + (t272 * t186 - t7 - t8 - t9) * t193 + (t10 + t11 + t12 + t271 * t185 + (t272 * t190 + t271 * t193) * t193) * t190; m(6) * (t190 * t48 - t193 * t49) + m(7) * (t190 * t40 - t193 * t41) + m(5) * (t190 * t59 - t193 * t60) + m(4) * (t190 * t98 - t193 * t99); m(6) * (t190 * t62 - t193 * t61) + m(7) * (t190 * t56 - t193 * t55) + m(5) * (t190 * t64 - t193 * t63) + m(4) * (-t125 * t193 + t190 * t126); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t226) * t227; t225 * t180 + m(6) * (t38 * t48 + t39 * t49) + m(7) * (t16 * t40 + t17 * t41) + m(5) * (t57 * t59 + t58 * t60) + (t190 * t196 + t193 * t195) * t179; m(6) * (t14 * t18 + t38 * t62 + t39 * t61) + m(7) * (t13 * t15 + t16 * t56 + t17 * t55) + m(5) * (t19 * t50 + t57 * t64 + t58 * t63) + ((t12 / 0.2e1 + t11 / 0.2e1 + t10 / 0.2e1) * t193 + (t9 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t190) * t179 + (t266 * t190 + t267 * t193) * t263 + t268 * t190 / 0.2e1 - t269 * t193 / 0.2e1; m(5) * (t57 * t190 - t193 * t58) + m(6) * (t38 * t190 - t193 * t39) + m(7) * (t16 * t190 - t17 * t193); m(7) * (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(6) * (t18 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t50 ^ 2 + t57 ^ 2 + t58 ^ 2) - t225 * t264 + (t268 * t193 + t269 * t190 + (t267 * t190 - t266 * t193) * t180) * t179; m(6) * (t151 * t49 + t153 * t48) + m(7) * (t151 * t41 + t153 * t40); m(6) * (t14 * t243 + t151 * t61 + t153 * t62) + m(7) * (t13 * t243 + t151 * t55 + t153 * t56); (-t151 * t193 + t153 * t190) * t277; m(7) * (t15 * t243 + t151 * t17 + t153 * t16) + m(6) * (t151 * t39 + t153 * t38 + t18 * t243); (t178 * t188 ^ 2 + t151 ^ 2 + t153 ^ 2) * t277; (-t190 * t41 - t193 * t40) * t260; m(7) * (t180 * t13 + (-t190 * t55 - t193 * t56) * t179); 0; m(7) * (t180 * t15 + (-t16 * t193 - t17 * t190) * t179); (-t151 * t190 - t153 * t193 + t180 * t188) * t260; m(7) * (t178 * t227 + t264);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
