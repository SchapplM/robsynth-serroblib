% Calculate kinetic energy for
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:08
% EndTime: 2019-03-09 02:56:11
% DurationCPUTime: 3.10s
% Computational Cost: add. (1178->286), mult. (1464->380), div. (0->0), fcn. (1246->8), ass. (0->149)
t323 = -Icges(5,4) - Icges(6,6);
t322 = Icges(5,1) + Icges(6,2);
t321 = Icges(5,2) + Icges(6,3);
t211 = qJ(3) + pkin(9);
t201 = cos(t211);
t320 = t323 * t201;
t200 = sin(t211);
t319 = t323 * t200;
t318 = -Icges(6,4) + Icges(5,5);
t317 = Icges(6,5) - Icges(5,6);
t316 = t321 * t201 - t319;
t315 = t322 * t200 - t320;
t314 = Icges(2,4) + Icges(3,6);
t313 = Icges(2,1) + Icges(3,2);
t312 = -Icges(3,4) + Icges(2,5);
t311 = Icges(3,5) - Icges(2,6);
t310 = Icges(2,2) + Icges(3,3);
t215 = sin(qJ(1));
t218 = cos(qJ(1));
t309 = t215 * t316 - t218 * t317;
t308 = -t215 * t317 - t218 * t316;
t307 = t315 * t215 + t218 * t318;
t306 = t215 * t318 - t315 * t218;
t305 = t321 * t200 + t320;
t304 = t322 * t201 + t319;
t303 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t214 = sin(qJ(3));
t217 = cos(qJ(3));
t302 = Icges(4,5) * t214 + Icges(4,6) * t217 + t200 * t318 - t317 * t201;
t301 = t314 * t218;
t300 = t314 * t215;
t279 = Icges(4,4) * t214;
t247 = Icges(4,2) * t217 + t279;
t138 = Icges(4,6) * t218 + t215 * t247;
t139 = Icges(4,6) * t215 - t218 * t247;
t278 = Icges(4,4) * t217;
t249 = Icges(4,1) * t214 + t278;
t140 = Icges(4,5) * t218 + t215 * t249;
t141 = Icges(4,5) * t215 - t218 * t249;
t179 = -Icges(4,2) * t214 + t278;
t184 = Icges(4,1) * t217 - t279;
t196 = qJD(3) * t215 + V_base(5);
t197 = qJD(3) * t218 + V_base(4);
t202 = V_base(6) + qJD(1);
t299 = t196 * (t139 * t217 + t141 * t214 + t306 * t200 + t308 * t201) + t197 * (t138 * t217 + t140 * t214 + t307 * t200 + t309 * t201) + (t179 * t217 + t184 * t214 + t200 * t304 - t201 * t305) * t202;
t298 = -t310 * t218 - t300;
t297 = t310 * t215 - t301;
t296 = t313 * t215 + t301;
t295 = t313 * t218 - t300;
t292 = (Icges(4,5) * t217 - Icges(4,6) * t214 + t317 * t200 + t201 * t318) * t202 + (t215 * t302 + t218 * t303) * t197 + (t215 * t303 - t218 * t302) * t196;
t285 = pkin(3) * t214;
t284 = pkin(3) * t217;
t283 = pkin(7) * t218;
t282 = t215 * pkin(7);
t271 = t200 * t215;
t270 = t200 * t218;
t213 = sin(qJ(6));
t269 = t213 * t218;
t268 = t215 * t213;
t216 = cos(qJ(6));
t267 = t215 * t216;
t266 = t216 * t218;
t250 = pkin(4) * t200 - qJ(5) * t201;
t142 = t250 * t215;
t147 = qJ(4) * t218 + t215 * t285;
t265 = -t142 - t147;
t264 = qJD(5) * t215;
t263 = qJD(6) * t200;
t188 = t215 * pkin(1) - qJ(2) * t218;
t262 = V_base(4) * t188 + V_base(3);
t261 = V_base(5) * pkin(6) + V_base(1);
t162 = pkin(4) * t201 + qJ(5) * t200;
t258 = -t162 - t284;
t257 = -t188 - t282;
t256 = qJD(2) * t215 + t261;
t146 = qJ(4) * t215 - t218 * t285;
t255 = -t146 + t257;
t254 = V_base(5) * pkin(2) + t256;
t253 = rSges(4,1) * t214 + rSges(4,2) * t217;
t252 = rSges(5,1) * t200 + rSges(5,2) * t201;
t251 = rSges(6,2) * t200 + rSges(6,3) * t201;
t192 = pkin(1) * t218 + t215 * qJ(2);
t231 = -qJD(2) * t218 + t202 * t192 + V_base(2);
t143 = t250 * t218;
t230 = t143 + t255;
t229 = qJD(4) * t218 + t196 * t284 + t254;
t228 = t196 * t162 + t229;
t224 = V_base(4) * t282 + (-t192 - t283) * V_base(5) + t262;
t223 = t202 * t283 + (-pkin(2) - pkin(6)) * V_base(4) + t231;
t222 = t197 * t146 + t224;
t221 = qJD(4) * t215 + t202 * t147 + t223;
t220 = qJD(5) * t200 - t197 * t143 + t222;
t219 = qJD(5) * t218 * t201 + t202 * t142 + t221;
t194 = rSges(2,1) * t218 - t215 * rSges(2,2);
t193 = -rSges(3,2) * t218 + t215 * rSges(3,3);
t191 = rSges(4,1) * t217 - rSges(4,2) * t214;
t190 = t215 * rSges(2,1) + rSges(2,2) * t218;
t189 = -t215 * rSges(3,2) - rSges(3,3) * t218;
t171 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t170 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t169 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t168 = qJD(6) * t201 + t202;
t165 = pkin(5) * t218 + pkin(8) * t271;
t164 = rSges(5,1) * t201 - rSges(5,2) * t200;
t163 = -rSges(6,2) * t201 + rSges(6,3) * t200;
t161 = t215 * pkin(5) - pkin(8) * t270;
t153 = -t201 * t268 + t266;
t152 = -t201 * t267 - t269;
t151 = t201 * t269 + t267;
t150 = t201 * t266 - t268;
t149 = t215 * t263 + t197;
t148 = -t218 * t263 + t196;
t145 = t215 * rSges(4,3) - t218 * t253;
t144 = rSges(4,3) * t218 + t215 * t253;
t133 = rSges(6,1) * t218 - t215 * t251;
t132 = t215 * rSges(6,1) + t218 * t251;
t131 = t215 * rSges(5,3) - t218 * t252;
t130 = rSges(5,3) * t218 + t215 * t252;
t116 = rSges(7,3) * t201 + (rSges(7,1) * t213 + rSges(7,2) * t216) * t200;
t115 = V_base(5) * rSges(2,3) - t190 * t202 + t261;
t114 = t194 * t202 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t113 = Icges(7,5) * t201 + (Icges(7,1) * t213 + Icges(7,4) * t216) * t200;
t112 = Icges(7,6) * t201 + (Icges(7,4) * t213 + Icges(7,2) * t216) * t200;
t111 = Icges(7,3) * t201 + (Icges(7,5) * t213 + Icges(7,6) * t216) * t200;
t108 = t190 * V_base(4) - t194 * V_base(5) + V_base(3);
t107 = V_base(5) * rSges(3,1) + (-t188 - t189) * t202 + t256;
t106 = t202 * t193 + (-rSges(3,1) - pkin(6)) * V_base(4) + t231;
t105 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t271;
t104 = t151 * rSges(7,1) + t150 * rSges(7,2) - rSges(7,3) * t270;
t103 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t271;
t102 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t270;
t101 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t271;
t100 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t270;
t99 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t271;
t98 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t270;
t97 = t189 * V_base(4) + (-t192 - t193) * V_base(5) + t262;
t96 = t191 * t196 + (-t145 + t257) * t202 + t254;
t95 = t202 * t144 - t197 * t191 + t223;
t94 = -t196 * t144 + t197 * t145 + t224;
t93 = t164 * t196 + (-t131 + t255) * t202 + t229;
t92 = t202 * t130 + (-t164 - t284) * t197 + t221;
t91 = t197 * t131 + (-t130 - t147) * t196 + t222;
t90 = -t201 * t264 + t163 * t196 + (-t132 + t230) * t202 + t228;
t89 = t202 * t133 + (-t163 + t258) * t197 + t219;
t88 = t197 * t132 + (-t133 + t265) * t196 + t220;
t87 = -t104 * t168 + t116 * t148 + (pkin(8) * t196 - t264) * t201 + (-t161 + t230) * t202 + t228;
t86 = t168 * t105 - t149 * t116 + t202 * t165 + (-pkin(8) * t201 + t258) * t197 + t219;
t85 = t149 * t104 - t148 * t105 + t197 * t161 + (-t165 + t265) * t196 + t220;
t1 = t149 * ((t101 * t152 + t103 * t153 + t99 * t271) * t149 + (t100 * t152 + t102 * t153 + t271 * t98) * t148 + (t111 * t271 + t112 * t152 + t113 * t153) * t168) / 0.2e1 + t148 * ((t150 * t101 + t151 * t103 - t270 * t99) * t149 + (t150 * t100 + t151 * t102 - t98 * t270) * t148 + (-t111 * t270 + t150 * t112 + t151 * t113) * t168) / 0.2e1 + m(2) * (t108 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t168 * ((t111 * t168 + t148 * t98 + t149 * t99) * t201 + ((t101 * t216 + t103 * t213) * t149 + (t100 * t216 + t102 * t213) * t148 + (t112 * t216 + t113 * t213) * t168) * t200) / 0.2e1 + m(1) * (t169 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + (t292 * t215 - t299 * t218) * t196 / 0.2e1 + (t299 * t215 + t292 * t218) * t197 / 0.2e1 + ((t215 * t298 + t296 * t218 + Icges(1,4)) * V_base(5) + (t297 * t215 + t295 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t215 - t298 * t218 + Icges(1,2)) * V_base(5) + (t215 * t295 - t297 * t218 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t138 * t214 + t140 * t217 - t309 * t200 + t307 * t201) * t197 + (-t139 * t214 + t141 * t217 - t308 * t200 + t306 * t201) * t196 + (-t179 * t214 + t184 * t217 + t200 * t305 + t201 * t304 + Icges(3,1) + Icges(2,3)) * t202) * t202 / 0.2e1 + t202 * V_base(5) * (t312 * t215 - t311 * t218) + t202 * V_base(4) * (t311 * t215 + t312 * t218) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
