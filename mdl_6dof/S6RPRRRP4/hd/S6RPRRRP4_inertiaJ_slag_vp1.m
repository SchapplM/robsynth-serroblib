% Calculate joint inertia matrix for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:59
% EndTime: 2019-03-09 06:07:07
% DurationCPUTime: 3.82s
% Computational Cost: add. (9463->382), mult. (8242->558), div. (0->0), fcn. (8746->10), ass. (0->186)
t183 = pkin(10) + qJ(3);
t178 = qJ(4) + t183;
t173 = cos(t178);
t192 = cos(qJ(5));
t193 = cos(qJ(1));
t232 = t192 * t193;
t190 = sin(qJ(5));
t191 = sin(qJ(1));
t234 = t191 * t190;
t143 = -t173 * t234 - t232;
t233 = t191 * t192;
t235 = t190 * t193;
t144 = t173 * t233 - t235;
t172 = sin(t178);
t239 = t172 * t191;
t86 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t239;
t88 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t239;
t285 = t86 + t88;
t145 = -t173 * t235 + t233;
t146 = t173 * t232 + t234;
t237 = t172 * t193;
t87 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t237;
t89 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t237;
t284 = t87 + t89;
t90 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t239;
t92 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t239;
t283 = t90 + t92;
t91 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t237;
t93 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t237;
t282 = t91 + t93;
t94 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t239;
t96 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t239;
t281 = t94 + t96;
t95 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t237;
t97 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t237;
t280 = t95 + t97;
t188 = -qJ(6) - pkin(9);
t279 = rSges(7,3) - t188;
t278 = t283 * t143 + t281 * t144 + t285 * t239;
t277 = t282 * t143 + t280 * t144 + t284 * t239;
t276 = t283 * t145 + t281 * t146 + t285 * t237;
t275 = t282 * t145 + t280 * t146 + t284 * t237;
t111 = -Icges(7,3) * t173 + (Icges(7,5) * t192 - Icges(7,6) * t190) * t172;
t113 = -Icges(7,6) * t173 + (Icges(7,4) * t192 - Icges(7,2) * t190) * t172;
t115 = -Icges(7,5) * t173 + (Icges(7,1) * t192 - Icges(7,4) * t190) * t172;
t51 = t111 * t239 + t113 * t143 + t115 * t144;
t112 = -Icges(6,3) * t173 + (Icges(6,5) * t192 - Icges(6,6) * t190) * t172;
t114 = -Icges(6,6) * t173 + (Icges(6,4) * t192 - Icges(6,2) * t190) * t172;
t116 = -Icges(6,5) * t173 + (Icges(6,1) * t192 - Icges(6,4) * t190) * t172;
t52 = t112 * t239 + t114 * t143 + t116 * t144;
t274 = -t51 - t52;
t53 = t111 * t237 + t145 * t113 + t146 * t115;
t54 = t112 * t237 + t145 * t114 + t146 * t116;
t273 = -t53 - t54;
t175 = pkin(5) * t192 + pkin(4);
t236 = t173 * t193;
t272 = t146 * rSges(7,1) + t145 * rSges(7,2) + pkin(5) * t234 + t175 * t236 + t237 * t279;
t184 = t191 ^ 2;
t271 = t274 * t173 + (t191 * t278 + t277 * t193) * t172;
t270 = t273 * t173 + (t191 * t276 + t193 * t275) * t172;
t269 = t277 * t191 - t193 * t278;
t268 = t191 * t275 - t193 * t276;
t212 = -t144 * rSges(7,1) - t143 * rSges(7,2);
t255 = pkin(9) + t188;
t256 = -pkin(4) + t175;
t253 = -pkin(5) * t235 + (-t172 * t255 + t173 * t256) * t191 + rSges(7,3) * t239 - t212;
t227 = pkin(4) * t236 + pkin(9) * t237;
t267 = -t227 + t272;
t266 = -t111 - t112;
t265 = (-t113 - t114) * t190;
t264 = (t115 + t116) * t172 * t192;
t206 = Icges(5,5) * t173 - Icges(5,6) * t172;
t125 = -Icges(5,3) * t193 + t191 * t206;
t126 = Icges(5,3) * t191 + t193 * t206;
t185 = t193 ^ 2;
t241 = Icges(5,4) * t173;
t208 = -Icges(5,2) * t172 + t241;
t128 = Icges(5,6) * t191 + t193 * t208;
t242 = Icges(5,4) * t172;
t210 = Icges(5,1) * t173 - t242;
t130 = Icges(5,5) * t191 + t193 * t210;
t204 = -t128 * t172 + t130 * t173;
t127 = -Icges(5,6) * t193 + t191 * t208;
t129 = -Icges(5,5) * t193 + t191 * t210;
t205 = t127 * t172 - t129 * t173;
t263 = -t185 * t125 - (t204 * t191 + (-t126 + t205) * t193) * t191 - t269;
t262 = t173 ^ 2;
t260 = t191 / 0.2e1;
t259 = -t193 / 0.2e1;
t176 = sin(t183);
t258 = pkin(3) * t176;
t257 = pkin(4) * t173;
t189 = -pkin(7) - qJ(2);
t187 = cos(pkin(10));
t174 = t187 * pkin(2) + pkin(1);
t254 = t172 * t265 + t266 * t173 + t264;
t177 = cos(t183);
t252 = rSges(4,1) * t177;
t251 = rSges(4,2) * t176;
t41 = -t173 * t86 + (-t190 * t90 + t192 * t94) * t172;
t250 = t41 * t193;
t42 = -t173 * t87 + (-t190 * t91 + t192 * t95) * t172;
t249 = t42 * t191;
t43 = -t173 * t88 + (-t190 * t92 + t192 * t96) * t172;
t248 = t43 * t193;
t44 = -t173 * t89 + (-t190 * t93 + t192 * t97) * t172;
t247 = t44 * t191;
t246 = rSges(3,3) + qJ(2);
t244 = Icges(4,4) * t176;
t243 = Icges(4,4) * t177;
t160 = pkin(3) * t177 + t174;
t153 = t193 * t160;
t231 = t193 * (-t193 * t174 + t153) + (t160 - t174) * t184;
t230 = (t255 - rSges(7,3)) * t173 + (rSges(7,1) * t192 - rSges(7,2) * t190 + t256) * t172;
t118 = -t173 * rSges(6,3) + (rSges(6,1) * t192 - rSges(6,2) * t190) * t172;
t152 = t172 * pkin(4) - t173 * pkin(9);
t229 = -t118 - t152;
t199 = rSges(5,1) * t236 - rSges(5,2) * t237 + t191 * rSges(5,3);
t214 = rSges(5,1) * t173 - rSges(5,2) * t172;
t74 = t191 * (-t193 * rSges(5,3) + t191 * t214) + t193 * t199;
t228 = t184 * (pkin(9) * t172 + t257) + t193 * t227;
t226 = t191 * rSges(4,3) + t193 * t252;
t224 = t184 + t185;
t223 = (t184 * t126 + (t205 * t193 + (-t125 + t204) * t191) * t193 + t268) * t191;
t222 = -t152 - t230;
t101 = t146 * rSges(6,1) + t145 * rSges(6,2) + rSges(6,3) * t237;
t151 = rSges(5,1) * t172 + rSges(5,2) * t173;
t219 = -t151 - t258;
t218 = -t152 - t258;
t182 = -pkin(8) + t189;
t217 = -t191 * t182 + t153;
t213 = -t144 * rSges(6,1) - t143 * rSges(6,2);
t99 = rSges(6,3) * t239 - t213;
t45 = t193 * t101 + t191 * t99 + t228;
t216 = -t118 + t218;
t215 = -t251 + t252;
t211 = Icges(4,1) * t177 - t244;
t209 = -Icges(4,2) * t176 + t243;
t207 = Icges(4,5) * t177 - Icges(4,6) * t176;
t149 = Icges(5,2) * t173 + t242;
t150 = Icges(5,1) * t172 + t241;
t201 = -t149 * t172 + t150 * t173;
t200 = t218 - t230;
t22 = t253 * t191 + t267 * t193 + t228;
t186 = sin(pkin(10));
t198 = rSges(3,1) * t187 - rSges(3,2) * t186 + pkin(1);
t196 = t263 * t193 + t223;
t195 = -(t249 - t250 + t247 - t248) * t173 / 0.2e1 + t270 * t260 + t271 * t259 + t269 * t239 / 0.2e1 + t268 * t237 / 0.2e1;
t148 = Icges(5,5) * t172 + Icges(5,6) * t173;
t194 = -t250 / 0.2e1 + t249 / 0.2e1 - t248 / 0.2e1 + t247 / 0.2e1 + (t128 * t173 + t130 * t172 + t191 * t148 + t193 * t201 - t273) * t260 + (t127 * t173 + t129 * t172 - t193 * t148 + t191 * t201 - t274) * t259;
t167 = rSges(2,1) * t193 - t191 * rSges(2,2);
t166 = -t191 * rSges(2,1) - rSges(2,2) * t193;
t159 = rSges(4,1) * t176 + rSges(4,2) * t177;
t134 = Icges(4,3) * t191 + t193 * t207;
t133 = -Icges(4,3) * t193 + t191 * t207;
t124 = t191 * t246 + t193 * t198;
t123 = -t191 * t198 + t193 * t246;
t120 = t219 * t193;
t119 = t219 * t191;
t109 = -t191 * t189 + (t174 - t251) * t193 + t226;
t108 = (rSges(4,3) - t189) * t193 + (-t174 - t215) * t191;
t103 = t199 + t217;
t102 = (rSges(5,3) - t182) * t193 + (-t160 - t214) * t191;
t85 = t229 * t193;
t84 = t229 * t191;
t83 = t193 * (-t193 * t251 + t226) + (-t193 * rSges(4,3) + t191 * t215) * t191;
t73 = t216 * t193;
t72 = t216 * t191;
t67 = t222 * t193;
t66 = t222 * t191;
t65 = t217 + t101 + t227;
t64 = -t193 * t182 + (-t257 - t160 + (-rSges(6,3) - pkin(9)) * t172) * t191 + t213;
t63 = -t173 * t101 - t118 * t237;
t62 = t118 * t239 + t173 * t99;
t61 = t200 * t193;
t60 = t200 * t191;
t59 = t217 + t272;
t58 = (pkin(5) * t190 - t182) * t193 + (-t172 * t279 - t173 * t175 - t160) * t191 + t212;
t55 = (-t101 * t191 + t193 * t99) * t172;
t46 = t74 + t231;
t28 = -t173 * t267 - t230 * t237;
t27 = t173 * t253 + t230 * t239;
t24 = t45 + t231;
t23 = (-t191 * t267 + t193 * t253) * t172;
t21 = t22 + t231;
t1 = [Icges(3,2) * t187 ^ 2 + t177 * (Icges(4,2) * t177 + t244) + t176 * (Icges(4,1) * t176 + t243) + Icges(2,3) + (Icges(3,1) * t186 + 0.2e1 * Icges(3,4) * t187) * t186 + (t149 + t266) * t173 + (t150 + t265) * t172 + m(7) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2) + m(3) * (t123 ^ 2 + t124 ^ 2) + m(2) * (t166 ^ 2 + t167 ^ 2) + t264; m(7) * (t191 * t58 - t193 * t59) + m(6) * (t191 * t64 - t193 * t65) + m(5) * (t191 * t102 - t103 * t193) + m(4) * (t191 * t108 - t109 * t193) + m(3) * (t191 * t123 - t124 * t193); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t224; m(4) * (-t108 * t193 - t109 * t191) * t159 + t194 + m(7) * (t58 * t61 + t59 * t60) + m(6) * (t64 * t73 + t65 * t72) + m(5) * (t102 * t120 + t103 * t119) + (t185 / 0.2e1 + t184 / 0.2e1) * (Icges(4,5) * t176 + Icges(4,6) * t177) + (t177 * (Icges(4,6) * t191 + t193 * t209) + t176 * (Icges(4,5) * t191 + t193 * t211)) * t260 + (t177 * (-Icges(4,6) * t193 + t191 * t209) + t176 * (-Icges(4,5) * t193 + t191 * t211)) * t259; m(5) * (-t119 * t193 + t120 * t191) + m(6) * (t73 * t191 - t193 * t72) + m(7) * (t61 * t191 - t193 * t60); m(7) * (t21 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t24 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t46 ^ 2) + t191 * t184 * t134 + m(4) * (t159 ^ 2 * t224 + t83 ^ 2) + t223 + (-t185 * t133 + (-t191 * t133 + t193 * t134) * t191 + t263) * t193; m(7) * (t58 * t67 + t59 * t66) + m(6) * (t64 * t85 + t65 * t84) + m(5) * (-t102 * t193 - t103 * t191) * t151 + t194; m(6) * (t85 * t191 - t193 * t84) + m(7) * (t67 * t191 - t193 * t66); m(7) * (t21 * t22 + t60 * t66 + t61 * t67) + m(6) * (t24 * t45 + t72 * t84 + t73 * t85) + m(5) * (t74 * t46 + (-t119 * t191 - t120 * t193) * t151) + t196; m(7) * (t22 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t45 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t151 ^ 2 * t224 + t74 ^ 2) + t196; -t254 * t173 + m(7) * (t27 * t58 + t28 * t59) + m(6) * (t62 * t64 + t63 * t65) + ((t42 / 0.2e1 + t54 / 0.2e1 + t53 / 0.2e1 + t44 / 0.2e1) * t193 + (t52 / 0.2e1 + t51 / 0.2e1 + t43 / 0.2e1 + t41 / 0.2e1) * t191) * t172; m(6) * (t62 * t191 - t193 * t63) + m(7) * (t27 * t191 - t193 * t28); m(7) * (t21 * t23 + t27 * t61 + t28 * t60) + m(6) * (t24 * t55 + t62 * t73 + t63 * t72) + t195; m(7) * (t22 * t23 + t27 * t67 + t28 * t66) + m(6) * (t45 * t55 + t62 * t85 + t63 * t84) + t195; m(7) * (t23 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t55 ^ 2 + t62 ^ 2 + t63 ^ 2) + t254 * t262 + (t270 * t193 + t271 * t191 + ((-t42 - t44) * t193 + (-t41 - t43) * t191) * t173) * t172; m(7) * (t191 * t59 + t193 * t58) * t172; 0; m(7) * (-t173 * t21 + (t191 * t60 + t193 * t61) * t172); m(7) * (-t173 * t22 + (t191 * t66 + t193 * t67) * t172); m(7) * (-t173 * t23 + (t191 * t28 + t193 * t27) * t172); m(7) * (t172 ^ 2 * t224 + t262);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
