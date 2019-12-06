% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:49
% EndTime: 2019-12-05 15:05:03
% DurationCPUTime: 7.34s
% Computational Cost: add. (17919->451), mult. (25625->738), div. (0->0), fcn. (30894->10), ass. (0->234)
t207 = qJ(3) + pkin(9);
t204 = sin(t207);
t205 = cos(t207);
t208 = sin(pkin(8));
t176 = (-Icges(5,5) * t204 - Icges(5,6) * t205) * t208;
t214 = sin(qJ(3));
t216 = cos(qJ(3));
t194 = (-Icges(4,5) * t214 - Icges(4,6) * t216) * t208;
t319 = -t176 - t194;
t316 = -2 * Icges(4,4);
t314 = -2 * Icges(5,4);
t313 = 2 * Icges(5,4);
t318 = Icges(4,1) - Icges(4,2);
t317 = -Icges(5,1) + Icges(5,2);
t315 = 2 * Icges(4,4);
t211 = cos(pkin(7));
t209 = sin(pkin(7));
t210 = cos(pkin(8));
t254 = t209 * t210;
t170 = t204 * t254 + t205 * t211;
t249 = t211 * t204;
t171 = t205 * t254 - t249;
t247 = t211 * t216;
t252 = t209 * t214;
t190 = -t210 * t252 - t247;
t248 = t211 * t214;
t251 = t209 * t216;
t191 = t210 * t251 - t248;
t311 = Icges(4,5) * t190 - Icges(5,5) * t170 - Icges(4,6) * t191 - Icges(5,6) * t171;
t172 = -t209 * t205 + t210 * t249;
t250 = t210 * t211;
t173 = t204 * t209 + t205 * t250;
t192 = -t210 * t248 + t251;
t193 = t210 * t247 + t252;
t310 = Icges(4,5) * t192 - Icges(5,5) * t172 - Icges(4,6) * t193 - Icges(5,6) * t173;
t257 = t208 * t211;
t258 = t208 * t209;
t299 = t209 * (Icges(4,5) * t258 + t190 * t315 + t318 * t191) + t211 * (Icges(4,5) * t257 + t192 * t315 + t318 * t193);
t283 = pkin(3) * t216;
t221 = qJ(4) * t208 + t210 * t283;
t136 = pkin(3) * t252 + t211 * t221;
t312 = -rSges(5,1) * t173 + rSges(5,2) * t172 - rSges(5,3) * t257 - t136;
t307 = Icges(5,6) * t210 + (t317 * t204 + t205 * t314) * t208;
t306 = Icges(5,5) * t210 + (t204 * t313 + t317 * t205) * t208;
t305 = Icges(4,6) * t210 + (-t318 * t214 + t216 * t316) * t208;
t304 = -Icges(4,5) * t210 + (t214 * t316 + t318 * t216) * t208;
t149 = rSges(4,1) * t190 - rSges(4,2) * t191;
t197 = (-rSges(4,1) * t214 - rSges(4,2) * t216) * t208;
t105 = t149 * t210 + t197 * t258;
t150 = rSges(4,1) * t192 - rSges(4,2) * t193;
t106 = -t150 * t210 - t197 * t257;
t303 = t105 * t209 - t106 * t211;
t302 = (-Icges(5,6) * t257 + t317 * t172 + t173 * t314) * t211 + (-Icges(5,6) * t258 + t317 * t170 + t171 * t314) * t209;
t301 = (-Icges(5,5) * t257 + t172 * t313 + t317 * t173) * t211 + (-Icges(5,5) * t258 + t170 * t313 + t317 * t171) * t209;
t300 = (-Icges(4,6) * t257 + t318 * t192 + t193 * t316) * t211 + (-Icges(4,6) * t258 + t318 * t190 + t191 * t316) * t209;
t296 = 2 * qJD(3);
t295 = m(4) / 0.2e1;
t294 = m(5) / 0.2e1;
t293 = m(6) / 0.2e1;
t237 = pkin(3) * t248;
t135 = t209 * t221 - t237;
t113 = t135 * t257;
t213 = sin(qJ(5));
t215 = cos(qJ(5));
t255 = t208 * t215;
t159 = -t173 * t213 + t211 * t255;
t256 = t208 * t213;
t160 = t173 * t215 + t211 * t256;
t87 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t172;
t234 = -pkin(4) * t173 - pkin(6) * t172 - t136 - t87;
t157 = -t171 * t213 + t209 * t255;
t158 = t171 * t215 + t209 * t256;
t86 = rSges(6,1) * t158 + rSges(6,2) * t157 + rSges(6,3) * t170;
t273 = pkin(4) * t171 + pkin(6) * t170 + t86;
t45 = t113 + (t209 * t234 + t211 * t273) * t208;
t162 = -qJ(4) * t210 + t208 * t283;
t243 = t210 * t135 + t162 * t258;
t184 = -t205 * t256 - t210 * t215;
t185 = t205 * t255 - t210 * t213;
t259 = t204 * t208;
t117 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t259;
t244 = t117 + (pkin(4) * t205 + pkin(6) * t204) * t208;
t58 = t210 * t273 + t244 * t258 + t243;
t59 = t234 * t210 + (-t162 - t244) * t257;
t103 = rSges(6,1) * t157 - rSges(6,2) * t158;
t104 = rSges(6,1) * t159 - rSges(6,2) * t160;
t72 = (t103 * t211 - t104 * t209) * t208;
t142 = rSges(6,1) * t184 - rSges(6,2) * t185;
t76 = t103 * t210 + t142 * t258;
t77 = -t104 * t210 - t142 * t257;
t292 = m(6) * (t45 * t72 + t58 * t76 + t59 * t77);
t291 = t170 / 0.2e1;
t290 = t171 / 0.2e1;
t289 = t172 / 0.2e1;
t288 = t173 / 0.2e1;
t287 = t204 / 0.2e1;
t286 = t209 / 0.2e1;
t285 = -t210 / 0.2e1;
t284 = t211 / 0.2e1;
t151 = Icges(6,4) * t157;
t84 = Icges(6,1) * t158 + Icges(6,5) * t170 + t151;
t281 = -Icges(6,2) * t158 + t151 + t84;
t280 = m(6) * qJD(5);
t80 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t170;
t264 = Icges(6,4) * t158;
t82 = Icges(6,2) * t157 + Icges(6,6) * t170 + t264;
t53 = t184 * t82 + t185 * t84 + t259 * t80;
t81 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t172;
t263 = Icges(6,4) * t160;
t83 = Icges(6,2) * t159 + Icges(6,6) * t172 + t263;
t152 = Icges(6,4) * t159;
t85 = Icges(6,1) * t160 + Icges(6,5) * t172 + t152;
t54 = t184 * t83 + t185 * t85 + t259 * t81;
t114 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t259;
t262 = Icges(6,4) * t185;
t115 = Icges(6,2) * t184 + Icges(6,6) * t259 + t262;
t180 = Icges(6,4) * t184;
t116 = Icges(6,1) * t185 + Icges(6,5) * t259 + t180;
t67 = t114 * t259 + t115 * t184 + t116 * t185;
t279 = t205 * (t170 * t53 + t172 * t54 + t259 * t67);
t125 = -rSges(5,1) * t170 - rSges(5,2) * t171;
t188 = t190 * pkin(3);
t169 = t210 * t188;
t181 = (-rSges(5,1) * t204 - rSges(5,2) * t205) * t208;
t206 = t208 ^ 2;
t238 = pkin(3) * t206 * t214;
t78 = t125 * t210 + t169 + (t181 * t208 - t238) * t209;
t278 = t209 * t78;
t277 = -Icges(6,2) * t160 + t152 + t85;
t276 = Icges(6,1) * t157 - t264 - t82;
t275 = Icges(6,1) * t159 - t263 - t83;
t231 = -rSges(6,1) * t215 + rSges(6,2) * t213;
t95 = rSges(6,3) * t171 + t170 * t231;
t274 = -pkin(4) * t170 + pkin(6) * t171 + t95;
t253 = t209 * t211;
t246 = Icges(6,1) * t184 - t115 - t262;
t245 = -Icges(6,2) * t185 + t116 + t180;
t189 = t192 * pkin(3);
t242 = rSges(5,1) * t172 + rSges(5,2) * t173 - t189;
t156 = (rSges(6,3) * t205 + t204 * t231) * t208;
t241 = t156 + (-pkin(4) * t204 + pkin(6) * t205) * t208;
t240 = qJD(3) * t208;
t239 = qJD(3) * t210;
t97 = Icges(6,5) * t157 - Icges(6,6) * t158;
t38 = t184 * t281 + t185 * t276 + t259 * t97;
t98 = Icges(6,5) * t159 - Icges(6,6) * t160;
t39 = t184 * t277 + t185 * t275 + t259 * t98;
t139 = Icges(6,5) * t184 - Icges(6,6) * t185;
t60 = t139 * t259 + t184 * t245 + t185 * t246;
t18 = -t210 * t60 + (t209 * t38 + t211 * t39) * t208;
t228 = -Icges(6,5) * t215 + Icges(6,6) * t213;
t227 = Icges(6,3) * t171 + t170 * t228 + t213 * t82 - t215 * t84;
t229 = -Icges(6,4) * t215 + Icges(6,2) * t213;
t91 = Icges(6,6) * t171 + t170 * t229;
t230 = -Icges(6,1) * t215 + Icges(6,4) * t213;
t93 = Icges(6,5) * t171 + t170 * t230;
t32 = t184 * t91 + t185 * t93 + (t204 * t227 + t205 * t80) * t208;
t226 = Icges(6,3) * t173 + t172 * t228 + t213 * t83 - t215 * t85;
t92 = Icges(6,6) * t173 + t172 * t229;
t94 = Icges(6,5) * t173 + t172 * t230;
t33 = t184 * t92 + t185 * t94 + (t204 * t226 + t205 * t81) * t208;
t154 = (Icges(6,6) * t205 + t204 * t229) * t208;
t155 = (Icges(6,5) * t205 + t204 * t230) * t208;
t223 = t115 * t213 - t116 * t215 + (Icges(6,3) * t205 + t204 * t228) * t208;
t52 = t154 * t184 + t155 * t185 + (t114 * t205 + t204 * t223) * t208;
t5 = t170 * t32 + t171 * t53 + t172 * t33 + t173 * t54 + (t204 * t52 + t205 * t67) * t208;
t236 = t18 / 0.2e1 - t5 / 0.2e1;
t96 = rSges(6,3) * t173 + t172 * t231;
t235 = pkin(4) * t172 - pkin(6) * t173 - t189 - t96;
t34 = t157 * t281 + t158 * t276 + t170 * t97;
t35 = t157 * t277 + t158 * t275 + t170 * t98;
t46 = t139 * t170 + t157 * t245 + t158 * t246;
t11 = -t210 * t46 + (t209 * t34 + t211 * t35) * t208;
t36 = t159 * t281 + t160 * t276 + t172 * t97;
t37 = t159 * t277 + t160 * t275 + t172 * t98;
t47 = t139 * t172 + t159 * t245 + t160 * t246;
t12 = -t210 * t47 + (t209 * t36 + t211 * t37) * t208;
t225 = t11 * t286 + t12 * t284;
t224 = m(6) * (t209 * t76 - t211 * t77);
t56 = t117 * t171 + t156 * t170 + (-t204 * t95 - t205 * t86) * t208;
t57 = -t117 * t173 - t156 * t172 + (t204 * t96 + t205 * t87) * t208;
t222 = (t209 * t56 - t211 * t57) * t293;
t218 = (-t210 * t72 + (t209 * t77 + t211 * t76) * t208) * t293;
t42 = -t170 * t96 - t171 * t87 + t172 * t95 + t173 * t86;
t219 = m(6) * (-t210 * t42 + (t209 * t57 + t211 * t56) * t208);
t15 = t218 - t219 / 0.2e1;
t23 = t222 - t224 / 0.2e1;
t30 = 0.2e1 * (t42 / 0.4e1 - t72 / 0.4e1) * m(6);
t220 = -t30 * qJD(1) - t23 * qJD(2) + t15 * qJD(4);
t25 = t157 * t91 + t158 * t93 + t170 * t227 + t171 * t80;
t26 = t157 * t92 + t158 * t94 + t170 * t226 + t171 * t81;
t43 = t114 * t171 + t154 * t157 + t155 * t158 + t170 * t223;
t48 = t157 * t82 + t158 * t84 + t170 * t80;
t49 = t157 * t83 + t158 * t85 + t170 * t81;
t62 = t114 * t170 + t115 * t157 + t116 * t158;
t3 = t170 * t25 + t171 * t48 + t172 * t26 + t173 * t49 + (t204 * t43 + t205 * t62) * t208;
t27 = t159 * t91 + t160 * t93 + t172 * t227 + t173 * t80;
t28 = t159 * t92 + t160 * t94 + t172 * t226 + t173 * t81;
t44 = t114 * t173 + t154 * t159 + t155 * t160 + t172 * t223;
t50 = t159 * t82 + t160 * t84 + t172 * t80;
t51 = t159 * t83 + t160 * t85 + t172 * t81;
t63 = t114 * t172 + t115 * t159 + t116 * t160;
t4 = t170 * t27 + t171 * t50 + t172 * t28 + t173 * t51 + (t204 * t44 + t205 * t63) * t208;
t64 = -t170 * t87 + t172 * t86;
t69 = t117 * t170 - t259 * t86;
t70 = -t117 * t172 + t259 * t87;
t217 = m(6) * (t42 * t64 + t56 * t69 + t57 * t70) + (t170 * t48 + t172 * t49 + t259 * t62) * t290 + (t170 * t50 + t172 * t51 + t259 * t63) * t288 + t289 * t4 + t291 * t3;
t198 = t206 * t237;
t166 = t188 * t257;
t165 = -rSges(5,3) * t210 + (rSges(5,1) * t205 - rSges(5,2) * t204) * t208;
t138 = rSges(4,1) * t193 + rSges(4,2) * t192 + rSges(4,3) * t257;
t137 = rSges(4,1) * t191 + rSges(4,2) * t190 + rSges(4,3) * t258;
t111 = rSges(5,1) * t171 - rSges(5,2) * t170 + rSges(5,3) * t258;
t88 = (t149 * t211 - t150 * t209) * t208;
t79 = -t181 * t257 + t210 * t242 + t198;
t75 = t104 * t259 - t142 * t172;
t74 = -t103 * t259 + t142 * t170;
t73 = t166 + (t125 * t211 + t209 * t242) * t208;
t68 = t103 * t172 - t104 * t170;
t66 = t210 * t235 - t241 * t257 + t198;
t65 = t169 + t274 * t210 + (t208 * t241 - t238) * t209;
t55 = t166 + (t209 * t235 + t211 * t274) * t208;
t31 = (t42 + t72) * t293;
t24 = t224 / 0.2e1 + t222;
t17 = t170 * t38 + t172 * t39 + t259 * t60;
t14 = t218 + t219 / 0.2e1;
t13 = -t210 * t52 + (t209 * t32 + t211 * t33) * t208;
t10 = t170 * t36 + t172 * t37 + t259 * t47;
t9 = t170 * t34 + t172 * t35 + t259 * t46;
t7 = -t210 * t44 + (t209 * t27 + t211 * t28) * t208;
t6 = -t210 * t43 + (t209 * t25 + t211 * t26) * t208;
t2 = t18 * t285 + t208 * t225 + t292;
t1 = (t279 / 0.2e1 + t5 * t287) * t208 + t217;
t8 = [0, 0, t31 * qJD(5) + (t293 * t55 + t294 * t73 + t295 * t88) * t296, 0, t31 * qJD(3) + t280 * t68; 0, 0, t24 * qJD(5) + (t303 * t295 + (-t211 * t79 + t278) * t294 + (t209 * t65 - t211 * t66) * t293) * t296, 0, t24 * qJD(3) + (t209 * t74 - t211 * t75) * t280; -t30 * qJD(5), -t23 * qJD(5), t2 * qJD(5) + 0.4e1 * (m(6) * (t45 * t55 + t58 * t65 + t59 * t66) / 0.4e1 + m(5) * (t113 * t73 + t243 * t78) / 0.4e1) * qJD(3) + (-t13 / 0.2e1 + (-t176 / 0.2e1 - t194 / 0.2e1) * t210 ^ 2 + m(4) * (t137 * t105 - t138 * t106) + m(5) * (t111 * t78 + t312 * t79)) * t239 + (m(4) * ((t137 * t211 - t138 * t209) * t88 + t303 * (-t210 * rSges(4,3) + (rSges(4,1) * t216 - rSges(4,2) * t214) * t208)) + m(5) * (t165 * t278 + t312 * t73 * t209 + ((-t162 - t165) * t79 + t111 * t73) * t211) + ((t301 * t204 + t302 * t205 - t214 * t299 + t300 * t216) * t208 - t311 * t254 - t310 * t250 + (-t204 * t306 - t205 * t307 + t214 * t304 - t216 * t305) * t210) * t285 + (t6 + (t301 * t170 + t302 * t171 + t299 * t190 + t300 * t191 + (t209 ^ 2 * t311 + t253 * t310) * t208 + t319 * t254) * t208) * t286 + (t7 + (t301 * t172 + t302 * t173 + t299 * t192 + t300 * t193 + (t211 ^ 2 * t310 + t253 * t311) * t208 + t319 * t250) * t208) * t284 - (t170 * t306 + t171 * t307 + t190 * t304 + t305 * t191) * t254 / 0.2e1 - (t172 * t306 + t173 * t307 + t192 * t304 + t193 * t305) * t250 / 0.2e1) * t240, t15 * qJD(5), t2 * qJD(3) + t220 + (t17 * t285 + t11 * t291 + t12 * t289 + m(6) * (t45 * t68 + t58 * t74 + t59 * t75 + t64 * t72 + t69 * t76 + t70 * t77) - t217 + (t10 * t284 + t9 * t286 - t279 / 0.2e1 + t236 * t204) * t208) * qJD(5); 0, 0, t14 * qJD(5) + ((-t210 * t55 + (t209 * t66 + t211 * t65) * t208) * t293 + (-t210 * t73 + (t209 * t79 + t211 * t78) * t208) * t294) * t296, 0, t14 * qJD(3) + (-t210 * t68 + (t209 * t75 + t211 * t74) * t208) * t280; t30 * qJD(3), t23 * qJD(3), (t6 * t291 + t7 * t289 + m(6) * (t42 * t45 + t55 * t64 + t56 * t58 + t57 * t59 + t65 * t69 + t66 * t70) - t292) * qJD(3) + t1 * qJD(5) + (-t171 * t62 / 0.2e1 - t173 * t63 / 0.2e1 + t236) * t239 + (t4 * t284 + t3 * t286 + (t209 * t48 + t211 * t49) * t290 + (t209 * t50 + t211 * t51) * t288 + t13 * t287 + (t67 * t285 + (t209 * t53 + t211 * t54) * t208 / 0.2e1) * t205 - t225) * t240 - t220, -t15 * qJD(3), t1 * qJD(3) + (m(6) * (t64 * t68 + t69 * t74 + t70 * t75) + t10 * t289 + t9 * t291 + t17 * t259 / 0.2e1) * qJD(5);];
Cq = t8;
