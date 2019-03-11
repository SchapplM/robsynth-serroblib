% Calculate kinetic energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:19
% EndTime: 2019-03-09 15:28:22
% DurationCPUTime: 3.32s
% Computational Cost: add. (1485->287), mult. (1745->408), div. (0->0), fcn. (1527->8), ass. (0->153)
t322 = Icges(4,4) + Icges(6,4) - Icges(5,5);
t321 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t320 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t227 = qJ(2) + qJ(3);
t224 = cos(t227);
t319 = t322 * t224;
t223 = sin(t227);
t318 = t322 * t223;
t317 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t316 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t315 = t320 * t223 - t319;
t314 = t321 * t224 - t318;
t313 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t230 = sin(qJ(1));
t233 = cos(qJ(1));
t312 = t315 * t230 + t316 * t233;
t311 = -t316 * t230 + t315 * t233;
t310 = t314 * t230 - t317 * t233;
t309 = t317 * t230 + t314 * t233;
t308 = -t320 * t224 - t318;
t307 = t321 * t223 + t319;
t306 = -t316 * t223 + t317 * t224;
t192 = V_base(5) + (-qJD(2) - qJD(3)) * t233;
t217 = qJD(2) * t230 + V_base(4);
t193 = qJD(3) * t230 + t217;
t219 = V_base(6) + qJD(1);
t305 = (t308 * t223 + t307 * t224) * t219 + (t311 * t223 + t309 * t224) * t193 + (t312 * t223 + t310 * t224) * t192;
t304 = (t317 * t223 + t316 * t224) * t219 + (t313 * t230 + t306 * t233) * t193 + (t306 * t230 - t313 * t233) * t192;
t229 = sin(qJ(2));
t300 = pkin(2) * t229;
t299 = pkin(4) * t223;
t232 = cos(qJ(2));
t298 = pkin(2) * t232;
t296 = Icges(2,4) * t230;
t295 = Icges(3,4) * t229;
t294 = Icges(3,4) * t232;
t287 = t224 * t230;
t286 = t224 * t233;
t228 = sin(qJ(6));
t285 = t228 * t230;
t284 = t228 * t233;
t231 = cos(qJ(6));
t283 = t230 * t231;
t282 = t231 * t233;
t128 = -pkin(8) * t233 + t230 * t298;
t214 = t230 * pkin(1) - t233 * pkin(7);
t281 = -t128 - t214;
t265 = pkin(3) * t224 + qJ(4) * t223;
t168 = t265 * t233;
t177 = pkin(4) * t286 - qJ(5) * t230;
t280 = -t168 - t177;
t279 = qJD(4) * t223;
t278 = qJD(6) * t224;
t277 = V_base(5) * pkin(6) + V_base(1);
t167 = t265 * t230;
t274 = -t167 + t281;
t187 = pkin(3) * t223 - qJ(4) * t224;
t273 = -t187 - t299;
t216 = -qJD(2) * t233 + V_base(5);
t272 = t216 * t300 + t277;
t176 = pkin(4) * t287 + qJ(5) * t233;
t271 = -t176 + t274;
t270 = pkin(5) * t223 + pkin(9) * t224;
t269 = rSges(3,1) * t232 - rSges(3,2) * t229;
t268 = rSges(4,1) * t224 - rSges(4,2) * t223;
t267 = rSges(5,1) * t224 + rSges(5,3) * t223;
t266 = rSges(6,1) * t223 - rSges(6,2) * t224;
t264 = Icges(3,1) * t232 - t295;
t260 = -Icges(3,2) * t229 + t294;
t256 = Icges(3,5) * t232 - Icges(3,6) * t229;
t252 = t192 * t187 + t233 * t279 + t272;
t215 = t233 * pkin(1) + t230 * pkin(7);
t251 = -V_base(4) * pkin(6) + t219 * t215 + V_base(2);
t250 = V_base(4) * t214 - t215 * V_base(5) + V_base(3);
t246 = (-Icges(3,3) * t233 + t230 * t256) * t216 + (Icges(3,3) * t230 + t233 * t256) * t217 + (Icges(3,5) * t229 + Icges(3,6) * t232) * t219;
t245 = -qJD(5) * t230 + t192 * t299 + t252;
t129 = pkin(8) * t230 + t233 * t298;
t244 = t217 * t128 - t129 * t216 + t250;
t243 = t219 * t129 - t217 * t300 + t251;
t242 = t219 * t168 + t230 * t279 + t243;
t241 = -qJD(4) * t224 + t193 * t167 + t244;
t240 = qJD(5) * t233 + t219 * t177 + t242;
t239 = t193 * t176 + t241;
t159 = -Icges(3,6) * t233 + t230 * t260;
t160 = Icges(3,6) * t230 + t233 * t260;
t161 = -Icges(3,5) * t233 + t230 * t264;
t162 = Icges(3,5) * t230 + t233 * t264;
t203 = Icges(3,2) * t232 + t295;
t206 = Icges(3,1) * t229 + t294;
t235 = (-t160 * t229 + t162 * t232) * t217 + (-t159 * t229 + t161 * t232) * t216 + (-t203 * t229 + t206 * t232) * t219;
t225 = Icges(2,4) * t233;
t211 = rSges(2,1) * t233 - rSges(2,2) * t230;
t210 = rSges(2,1) * t230 + rSges(2,2) * t233;
t209 = rSges(3,1) * t229 + rSges(3,2) * t232;
t208 = Icges(2,1) * t233 - t296;
t207 = Icges(2,1) * t230 + t225;
t205 = -Icges(2,2) * t230 + t225;
t204 = Icges(2,2) * t233 + t296;
t199 = qJD(6) * t223 + t219;
t198 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t197 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t196 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t191 = -pkin(5) * t224 + pkin(9) * t223;
t190 = -rSges(6,1) * t224 - rSges(6,2) * t223;
t189 = rSges(4,1) * t223 + rSges(4,2) * t224;
t188 = rSges(5,1) * t223 - rSges(5,3) * t224;
t174 = t223 * t282 - t285;
t173 = -t223 * t284 - t283;
t172 = t223 * t283 + t284;
t171 = -t223 * t285 + t282;
t170 = t270 * t233;
t169 = t270 * t230;
t165 = rSges(3,3) * t230 + t233 * t269;
t164 = -rSges(3,3) * t233 + t230 * t269;
t156 = t233 * t278 + t193;
t155 = t230 * t278 + t192;
t154 = rSges(4,3) * t230 + t233 * t268;
t153 = rSges(5,2) * t230 + t233 * t267;
t152 = -rSges(6,3) * t230 + t233 * t266;
t151 = -rSges(4,3) * t233 + t230 * t268;
t150 = -rSges(5,2) * t233 + t230 * t267;
t149 = rSges(6,3) * t233 + t230 * t266;
t127 = rSges(7,3) * t223 + (-rSges(7,1) * t231 + rSges(7,2) * t228) * t224;
t126 = Icges(7,5) * t223 + (-Icges(7,1) * t231 + Icges(7,4) * t228) * t224;
t125 = Icges(7,6) * t223 + (-Icges(7,4) * t231 + Icges(7,2) * t228) * t224;
t124 = Icges(7,3) * t223 + (-Icges(7,5) * t231 + Icges(7,6) * t228) * t224;
t123 = V_base(5) * rSges(2,3) - t210 * t219 + t277;
t122 = t211 * t219 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t119 = t210 * V_base(4) - t211 * V_base(5) + V_base(3);
t115 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t286;
t114 = rSges(7,1) * t172 + rSges(7,2) * t171 + rSges(7,3) * t287;
t113 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t286;
t112 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t287;
t111 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t286;
t110 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t287;
t109 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t286;
t108 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t287;
t107 = t209 * t216 + (-t164 - t214) * t219 + t277;
t106 = t165 * t219 - t209 * t217 + t251;
t105 = t164 * t217 - t165 * t216 + t250;
t104 = t189 * t192 + (-t151 + t281) * t219 + t272;
t103 = t154 * t219 - t189 * t193 + t243;
t102 = t151 * t193 - t154 * t192 + t244;
t101 = t188 * t192 + (-t150 + t274) * t219 + t252;
t100 = t153 * t219 + (-t187 - t188) * t193 + t242;
t99 = t190 * t192 + (-t149 + t271) * t219 + t245;
t98 = t152 * t219 + (-t190 + t273) * t193 + t240;
t97 = t150 * t193 + (-t153 - t168) * t192 + t241;
t96 = t149 * t193 + (-t152 + t280) * t192 + t239;
t95 = -t114 * t199 + t127 * t155 + t191 * t192 + (-t169 + t271) * t219 + t245;
t94 = t115 * t199 - t127 * t156 + t170 * t219 + (-t191 + t273) * t193 + t240;
t93 = t114 * t156 - t115 * t155 + t169 * t193 + (-t170 + t280) * t192 + t239;
t1 = t199 * ((t108 * t155 + t109 * t156 + t124 * t199) * t223 + ((t111 * t228 - t113 * t231) * t156 + (t110 * t228 - t112 * t231) * t155 + (t125 * t228 - t126 * t231) * t199) * t224) / 0.2e1 + m(1) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(2) * (t119 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t96 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(7) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + t217 * (t246 * t230 + t235 * t233) / 0.2e1 + t216 * (t230 * t235 - t246 * t233) / 0.2e1 + t156 * ((t109 * t286 + t111 * t173 + t113 * t174) * t156 + (t108 * t286 + t110 * t173 + t112 * t174) * t155 + (t124 * t286 + t125 * t173 + t126 * t174) * t199) / 0.2e1 + t155 * ((t109 * t287 + t111 * t171 + t113 * t172) * t156 + (t108 * t287 + t110 * t171 + t112 * t172) * t155 + (t124 * t287 + t125 * t171 + t126 * t172) * t199) / 0.2e1 + ((-t204 * t230 + t207 * t233 + Icges(1,4)) * V_base(5) + (-t205 * t230 + t208 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t204 * t233 + t207 * t230 + Icges(1,2)) * V_base(5) + (t205 * t233 + t208 * t230 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t305 * t230 - t304 * t233) * t192 / 0.2e1 + (t304 * t230 + t305 * t233) * t193 / 0.2e1 + ((t160 * t232 + t162 * t229) * t217 + (t159 * t232 + t161 * t229) * t216 + (t309 * t223 - t311 * t224) * t193 + (t310 * t223 - t312 * t224) * t192 + (t203 * t232 + t206 * t229 + t307 * t223 - t308 * t224 + Icges(2,3)) * t219) * t219 / 0.2e1 + V_base(4) * t219 * (Icges(2,5) * t233 - Icges(2,6) * t230) + V_base(5) * t219 * (Icges(2,5) * t230 + Icges(2,6) * t233) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
