% Calculate joint inertia matrix for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:54
% EndTime: 2019-03-09 05:59:03
% DurationCPUTime: 4.32s
% Computational Cost: add. (12296->420), mult. (11210->581), div. (0->0), fcn. (12071->10), ass. (0->219)
t216 = qJ(1) + pkin(10);
t209 = sin(t216);
t210 = cos(t216);
t218 = qJ(4) + qJ(5);
t212 = cos(t218);
t211 = sin(t218);
t223 = cos(qJ(3));
t272 = t211 * t223;
t156 = -t209 * t272 - t210 * t212;
t270 = t212 * t223;
t157 = t209 * t270 - t210 * t211;
t220 = sin(qJ(3));
t278 = t209 * t220;
t101 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t278;
t103 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t278;
t324 = t101 + t103;
t158 = t209 * t212 - t210 * t272;
t159 = t209 * t211 + t210 * t270;
t275 = t210 * t220;
t102 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t275;
t104 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t275;
t323 = t102 + t104;
t105 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t278;
t107 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t278;
t322 = t105 + t107;
t106 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t275;
t108 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t275;
t321 = t106 + t108;
t109 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t278;
t111 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t278;
t320 = t109 + t111;
t110 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t275;
t112 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t275;
t319 = t110 + t112;
t318 = t322 * t156 + t320 * t157 + t324 * t278;
t317 = t321 * t156 + t319 * t157 + t323 * t278;
t316 = t322 * t158 + t320 * t159 + t324 * t275;
t315 = t321 * t158 + t319 * t159 + t323 * t275;
t160 = -Icges(7,3) * t223 + (Icges(7,5) * t212 - Icges(7,6) * t211) * t220;
t162 = -Icges(7,6) * t223 + (Icges(7,4) * t212 - Icges(7,2) * t211) * t220;
t164 = -Icges(7,5) * t223 + (Icges(7,1) * t212 - Icges(7,4) * t211) * t220;
t68 = t156 * t162 + t157 * t164 + t160 * t278;
t161 = -Icges(6,3) * t223 + (Icges(6,5) * t212 - Icges(6,6) * t211) * t220;
t163 = -Icges(6,6) * t223 + (Icges(6,4) * t212 - Icges(6,2) * t211) * t220;
t165 = -Icges(6,5) * t223 + (Icges(6,1) * t212 - Icges(6,4) * t211) * t220;
t69 = t156 * t163 + t157 * t165 + t161 * t278;
t314 = -t69 - t68;
t70 = t158 * t162 + t159 * t164 + t160 * t275;
t71 = t158 * t163 + t159 * t165 + t161 * t275;
t313 = -t70 - t71;
t308 = (t164 + t165) * t212 * t220;
t310 = (-t162 - t163) * t211;
t311 = -t160 - t161;
t299 = t220 * t310 + t223 * t311 + t308;
t312 = t299 * t223;
t309 = Icges(4,5) * t220;
t219 = sin(qJ(4));
t184 = t219 * pkin(4) + pkin(5) * t211;
t307 = -t157 * rSges(7,1) - t156 * rSges(7,2) + t210 * t184;
t306 = t309 / 0.2e1;
t305 = t314 * t223 + (t209 * t318 + t317 * t210) * t220;
t304 = t313 * t223 + (t209 * t316 + t210 * t315) * t220;
t303 = t317 * t209 - t210 * t318;
t302 = t209 * t315 - t210 * t316;
t54 = -t223 * t101 + (-t105 * t211 + t109 * t212) * t220;
t56 = -t223 * t103 + (-t107 * t211 + t111 * t212) * t220;
t301 = -t54 - t56;
t55 = -t223 * t102 + (-t106 * t211 + t110 * t212) * t220;
t57 = -t223 * t104 + (-t108 * t211 + t112 * t212) * t220;
t300 = t55 + t57;
t225 = -pkin(9) - pkin(8);
t266 = t220 * t225;
t276 = t210 * t219;
t255 = pkin(4) * t276 + t209 * t266;
t222 = cos(qJ(4));
t206 = t222 * pkin(4) + pkin(3);
t183 = pkin(5) * t212 + t206;
t256 = t183 - t206;
t215 = -qJ(6) + t225;
t269 = t215 * t220;
t298 = rSges(7,3) * t278 + (t256 * t223 - t269) * t209 + t255 - t307;
t252 = t215 - t225;
t297 = (t252 - rSges(7,3)) * t223 + (rSges(7,1) * t212 - rSges(7,2) * t211 + t256) * t220;
t274 = t210 * t223;
t296 = t159 * rSges(7,1) + t158 * rSges(7,2) + rSges(7,3) * t275 + t183 * t274 + t209 * t184;
t207 = t209 ^ 2;
t208 = t210 ^ 2;
t295 = t209 / 0.2e1;
t294 = -t210 / 0.2e1;
t293 = t210 / 0.2e1;
t292 = -t223 / 0.2e1;
t190 = t220 * rSges(4,1) + t223 * rSges(4,2);
t291 = m(4) * t190;
t290 = pkin(3) * t223;
t289 = pkin(8) * t220;
t221 = sin(qJ(1));
t288 = t221 * pkin(1);
t287 = -pkin(3) + t206;
t286 = t312 + (t301 * t209 - t300 * t210) * t220;
t285 = t298 * t275;
t283 = t210 * rSges(4,3);
t279 = t209 * t219;
t257 = -pkin(4) * t279 - t206 * t274;
t282 = -t252 * t275 + t257 + t296;
t237 = -t157 * rSges(6,1) - t156 * rSges(6,2);
t115 = rSges(6,3) * t278 - t237;
t167 = -t223 * rSges(6,3) + (rSges(6,1) * t212 - rSges(6,2) * t211) * t220;
t81 = t223 * t115 + t167 * t278;
t280 = Icges(4,4) * t223;
t176 = -Icges(5,6) * t223 + (Icges(5,4) * t222 - Icges(5,2) * t219) * t220;
t268 = t219 * t176;
t267 = t219 * t223;
t265 = t222 * t223;
t117 = t159 * rSges(6,1) + t158 * rSges(6,2) + rSges(6,3) * t275;
t229 = -t210 * t266 - t257;
t254 = pkin(3) * t274 + pkin(8) * t275;
t130 = t229 - t254;
t264 = -t117 - t130;
t129 = (t287 * t223 - t289) * t209 - t255;
t155 = (pkin(8) + t225) * t223 + t287 * t220;
t263 = t223 * t129 + t155 * t278;
t261 = t207 * (t289 + t290) + t210 * t254;
t260 = -t155 - t167;
t178 = -t223 * rSges(5,3) + (rSges(5,1) * t222 - rSges(5,2) * t219) * t220;
t197 = t220 * pkin(3) - t223 * pkin(8);
t258 = -t178 - t197;
t253 = t207 + t208;
t171 = -t209 * t267 - t210 * t222;
t172 = t209 * t265 - t276;
t121 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t278;
t123 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t278;
t125 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t278;
t60 = -t223 * t121 + (-t123 * t219 + t125 * t222) * t220;
t175 = -Icges(5,3) * t223 + (Icges(5,5) * t222 - Icges(5,6) * t219) * t220;
t177 = -Icges(5,5) * t223 + (Icges(5,1) * t222 - Icges(5,4) * t219) * t220;
t79 = t171 * t176 + t172 * t177 + t175 * t278;
t251 = t79 / 0.2e1 + t60 / 0.2e1;
t173 = t209 * t222 - t210 * t267;
t174 = t210 * t265 + t279;
t122 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t275;
t124 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t275;
t126 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t275;
t61 = -t223 * t122 + (-t124 * t219 + t126 * t222) * t220;
t80 = t173 * t176 + t174 * t177 + t175 * t275;
t250 = t80 / 0.2e1 + t61 / 0.2e1;
t249 = -t130 - t282;
t248 = t304 * t275 + t305 * t278;
t247 = -t155 - t297;
t246 = -t197 + t260;
t128 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t275;
t224 = cos(qJ(1));
t214 = t224 * pkin(1);
t245 = t210 * pkin(2) + t209 * pkin(7) + t214;
t244 = t278 / 0.2e1;
t243 = t275 / 0.2e1;
t242 = t210 * pkin(7) - t288;
t44 = t298 * t223 + t297 * t278;
t241 = t209 * t129 + t210 * t130 + t261;
t240 = -t197 + t247;
t239 = rSges(4,1) * t223 - rSges(4,2) * t220;
t238 = -t172 * rSges(5,1) - t171 * rSges(5,2);
t234 = -Icges(4,2) * t220 + t280;
t233 = Icges(4,5) * t223 - Icges(4,6) * t220;
t230 = rSges(4,1) * t274 - rSges(4,2) * t275 + t209 * rSges(4,3);
t228 = t286 * t223 + t248;
t227 = (-t301 - t314) * t244 + (t300 - t313) * t243;
t226 = t304 * t295 + t305 * t294 + (t300 * t209 + t301 * t210) * t292 + t303 * t244 + t302 * t243;
t192 = t224 * rSges(2,1) - t221 * rSges(2,2);
t191 = -t221 * rSges(2,1) - t224 * rSges(2,2);
t187 = Icges(4,6) * t223 + t309;
t180 = t210 * rSges(3,1) - t209 * rSges(3,2) + t214;
t179 = -t209 * rSges(3,1) - t210 * rSges(3,2) - t288;
t148 = t220 * t222 * t177;
t143 = Icges(4,3) * t209 + t233 * t210;
t142 = -Icges(4,3) * t210 + t233 * t209;
t135 = t258 * t210;
t134 = t258 * t209;
t133 = t230 + t245;
t132 = t283 + (-pkin(2) - t239) * t209 + t242;
t127 = rSges(5,3) * t278 - t238;
t113 = t129 * t275;
t100 = t210 * t230 + (t239 * t209 - t283) * t209;
t97 = t115 * t275;
t95 = t246 * t210;
t94 = t246 * t209;
t93 = -t223 * t175 - t220 * t268 + t148;
t90 = -t223 * t128 - t178 * t275;
t89 = t223 * t127 + t178 * t278;
t87 = t245 + t128 + t254;
t86 = (-t290 - pkin(2) + (-rSges(5,3) - pkin(8)) * t220) * t209 + t238 + t242;
t82 = -t223 * t117 - t167 * t275;
t78 = t229 + t245 + t117;
t77 = (-rSges(6,3) * t220 - t206 * t223 - pkin(2)) * t209 + t237 + t242 + t255;
t76 = t240 * t210;
t75 = t240 * t209;
t74 = (t127 * t210 - t128 * t209) * t220;
t73 = -t210 * t269 + t245 + t296;
t72 = (-t183 * t223 - pkin(2) + (-rSges(7,3) + t215) * t220) * t209 + t242 + t307;
t63 = -t117 * t278 + t97;
t62 = t209 * t127 + t210 * t128 + t261;
t59 = t264 * t223 + t260 * t275;
t58 = t263 + t81;
t53 = t122 * t275 + t173 * t124 + t174 * t126;
t52 = t121 * t275 + t173 * t123 + t174 * t125;
t51 = t122 * t278 + t171 * t124 + t172 * t126;
t50 = t121 * t278 + t171 * t123 + t172 * t125;
t45 = -t282 * t223 - t275 * t297;
t35 = t264 * t278 + t113 + t97;
t34 = t209 * t115 + t210 * t117 + t241;
t33 = t249 * t223 + t247 * t275;
t32 = t44 + t263;
t31 = -t282 * t278 + t285;
t28 = t53 * t209 - t52 * t210;
t27 = t51 * t209 - t50 * t210;
t26 = t249 * t278 + t113 + t285;
t25 = t298 * t209 + t282 * t210 + t241;
t14 = -t80 * t223 + (t209 * t52 + t210 * t53) * t220;
t13 = -t79 * t223 + (t209 * t50 + t210 * t51) * t220;
t1 = [Icges(2,3) + Icges(3,3) + t148 + (Icges(4,4) * t220 + Icges(4,2) * t223 - t175 + t311) * t223 + (Icges(4,1) * t220 - t268 + t280 + t310) * t220 + m(7) * (t72 ^ 2 + t73 ^ 2) + m(6) * (t77 ^ 2 + t78 ^ 2) + m(5) * (t86 ^ 2 + t87 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2) + m(3) * (t179 ^ 2 + t180 ^ 2) + m(2) * (t191 ^ 2 + t192 ^ 2) + t308; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t72 * t76 + t73 * t75) + m(6) * (t77 * t95 + t78 * t94) + m(5) * (t134 * t87 + t135 * t86) + (-t56 / 0.2e1 - t54 / 0.2e1 - t69 / 0.2e1 - t68 / 0.2e1 - t132 * t291 + (-Icges(4,6) * t210 + t234 * t209) * t292 + t210 * t306 + t187 * t293 - t251) * t210 + (t71 / 0.2e1 + t70 / 0.2e1 + t57 / 0.2e1 + t55 / 0.2e1 - t133 * t291 + t223 * (Icges(4,6) * t209 + t234 * t210) / 0.2e1 + t209 * t306 + t187 * t295 + t250) * t209; m(4) * t100 + m(5) * t62 + m(6) * t34 + m(7) * t25; m(7) * (t25 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t34 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(5) * (t134 ^ 2 + t135 ^ 2 + t62 ^ 2) + m(4) * (t253 * t190 ^ 2 + t100 ^ 2) + (-t208 * t142 - t27 - t303) * t210 + (t207 * t143 + t28 + (-t209 * t142 + t210 * t143) * t210 + t302) * t209; (-t93 - t299) * t223 + m(7) * (t32 * t72 + t33 * t73) + m(6) * (t58 * t77 + t59 * t78) + m(5) * (t86 * t89 + t87 * t90) + (t251 * t209 + t250 * t210) * t220 + t227; m(5) * t74 + m(6) * t35 + m(7) * t26; (t61 * t209 - t60 * t210) * t292 + t13 * t294 + t14 * t295 + t226 + m(7) * (t26 * t25 + t32 * t76 + t33 * t75) + m(6) * (t35 * t34 + t58 * t95 + t59 * t94) + m(5) * (t134 * t90 + t135 * t89 + t62 * t74) + (t27 * t295 + t28 * t293) * t220; (t93 * t223 + t286) * t223 + (t210 * t14 + t209 * t13 - t223 * (t209 * t60 + t210 * t61)) * t220 + m(7) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t74 ^ 2 + t89 ^ 2 + t90 ^ 2) + t248; -t312 + m(7) * (t44 * t72 + t45 * t73) + m(6) * (t77 * t81 + t78 * t82) + t227; m(6) * t63 + m(7) * t31; m(7) * (t31 * t25 + t44 * t76 + t45 * t75) + m(6) * (t63 * t34 + t81 * t95 + t82 * t94) + t226; m(7) * (t26 * t31 + t32 * t44 + t33 * t45) + m(6) * (t63 * t35 + t58 * t81 + t59 * t82) + t228; m(7) * (t31 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t63 ^ 2 + t81 ^ 2 + t82 ^ 2) + t228; m(7) * (t209 * t73 + t210 * t72) * t220; -m(7) * t223; m(7) * (-t223 * t25 + (t209 * t75 + t210 * t76) * t220); m(7) * (-t223 * t26 + (t209 * t33 + t210 * t32) * t220); m(7) * (-t223 * t31 + (t209 * t45 + t210 * t44) * t220); m(7) * (t253 * t220 ^ 2 + t223 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
