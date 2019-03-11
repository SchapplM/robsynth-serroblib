% Calculate kinetic energy for
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:11
% EndTime: 2019-03-09 02:15:14
% DurationCPUTime: 2.76s
% Computational Cost: add. (1339->285), mult. (1702->383), div. (0->0), fcn. (1550->8), ass. (0->151)
t325 = Icges(2,4) + Icges(3,6);
t324 = Icges(2,1) + Icges(3,2);
t323 = Icges(6,1) + Icges(7,1);
t322 = -Icges(3,4) + Icges(2,5);
t321 = Icges(6,4) - Icges(7,5);
t320 = Icges(7,4) + Icges(6,5);
t319 = Icges(3,5) - Icges(2,6);
t318 = Icges(2,2) + Icges(3,3);
t317 = Icges(6,2) + Icges(7,3);
t316 = Icges(7,6) - Icges(6,6);
t315 = Icges(6,3) + Icges(7,2);
t314 = rSges(7,1) + pkin(5);
t313 = rSges(7,3) + qJ(6);
t232 = cos(qJ(1));
t312 = t325 * t232;
t230 = sin(qJ(1));
t311 = t325 * t230;
t225 = pkin(9) + qJ(4);
t214 = sin(t225);
t231 = cos(qJ(5));
t273 = t232 * t231;
t229 = sin(qJ(5));
t275 = t230 * t229;
t171 = t214 * t275 - t273;
t274 = t230 * t231;
t277 = t229 * t232;
t172 = t214 * t274 + t277;
t215 = cos(t225);
t279 = t215 * t230;
t310 = t317 * t171 - t321 * t172 - t316 * t279;
t173 = t214 * t277 + t274;
t174 = -t214 * t273 + t275;
t278 = t215 * t232;
t309 = -t317 * t173 - t321 * t174 + t316 * t278;
t308 = t316 * t171 + t320 * t172 - t315 * t279;
t307 = -t316 * t173 + t320 * t174 + t315 * t278;
t306 = -t321 * t171 + t323 * t172 - t320 * t279;
t305 = t321 * t173 + t323 * t174 + t320 * t278;
t304 = (t317 * t229 - t321 * t231) * t215 + t316 * t214;
t303 = (t316 * t229 + t320 * t231) * t215 + t315 * t214;
t302 = (-t321 * t229 + t323 * t231) * t215 + t320 * t214;
t301 = t324 * t230 + t312;
t300 = t324 * t232 - t311;
t299 = t322 * t230 - t319 * t232;
t298 = t319 * t230 + t322 * t232;
t227 = cos(pkin(9));
t226 = sin(pkin(9));
t286 = Icges(4,4) * t226;
t253 = Icges(4,2) * t227 + t286;
t159 = Icges(4,6) * t230 - t232 * t253;
t285 = Icges(4,4) * t227;
t255 = Icges(4,1) * t226 + t285;
t161 = Icges(4,5) * t230 - t232 * t255;
t297 = t159 * t227 + t161 * t226 - t318 * t232 - t311;
t158 = Icges(4,6) * t232 + t230 * t253;
t160 = Icges(4,5) * t232 + t230 * t255;
t296 = t158 * t227 + t160 * t226 + t318 * t230 - t312;
t284 = Icges(5,4) * t214;
t252 = Icges(5,2) * t215 + t284;
t147 = Icges(5,6) * t232 + t230 * t252;
t148 = Icges(5,6) * t230 - t232 * t252;
t283 = Icges(5,4) * t215;
t254 = Icges(5,1) * t214 + t283;
t149 = Icges(5,5) * t232 + t230 * t254;
t150 = Icges(5,5) * t230 - t232 * t254;
t177 = -Icges(5,2) * t214 + t283;
t178 = Icges(5,1) * t215 - t284;
t209 = qJD(4) * t230 + V_base(5);
t210 = qJD(4) * t232 + V_base(4);
t216 = V_base(6) + qJD(1);
t295 = (t147 * t215 + t149 * t214) * t210 + (t148 * t215 + t150 * t214) * t209 + (t177 * t215 + t178 * t214) * t216;
t294 = -pkin(2) - pkin(6);
t289 = pkin(3) * t226;
t288 = pkin(3) * t227;
t280 = qJ(3) * t232;
t276 = t230 * qJ(3);
t271 = -rSges(7,2) * t279 + t313 * t171 + t172 * t314;
t270 = rSges(7,2) * t278 - t313 * t173 + t174 * t314;
t269 = rSges(7,2) * t214 + (t313 * t229 + t231 * t314) * t215;
t268 = qJD(5) * t215;
t202 = t230 * pkin(1) - qJ(2) * t232;
t267 = V_base(4) * t202 + V_base(3);
t266 = V_base(5) * pkin(6) + V_base(1);
t263 = V_base(4) * t276 + t267;
t262 = qJD(2) * t230 + t266;
t261 = -t202 - t276;
t205 = pkin(1) * t232 + t230 * qJ(2);
t260 = -t205 - t280;
t259 = pkin(4) * t214 - pkin(8) * t215;
t167 = pkin(7) * t230 - t232 * t289;
t258 = -t167 + t261;
t257 = rSges(4,1) * t226 + rSges(4,2) * t227;
t256 = rSges(5,1) * t214 + rSges(5,2) * t215;
t251 = Icges(4,5) * t226 + Icges(4,6) * t227;
t250 = Icges(5,5) * t214 + Icges(5,6) * t215;
t187 = -Icges(4,2) * t226 + t285;
t188 = Icges(4,1) * t227 - t286;
t244 = t187 * t227 + t188 * t226;
t243 = -qJD(2) * t232 + t216 * t205 + V_base(2);
t242 = V_base(5) * pkin(2) + qJD(3) * t232 + t262;
t241 = V_base(5) * t288 + t242;
t240 = qJD(3) * t230 + t216 * t280 + t243;
t239 = (Icges(5,3) * t232 + t230 * t250) * t210 + (Icges(5,3) * t230 - t232 * t250) * t209 + (Icges(5,5) * t215 - Icges(5,6) * t214) * t216;
t238 = (Icges(4,3) * t232 + t230 * t251) * V_base(4) + (Icges(4,3) * t230 - t232 * t251) * V_base(5) + (Icges(4,5) * t227 - Icges(4,6) * t226) * t216;
t168 = pkin(7) * t232 + t230 * t289;
t237 = V_base(4) * t167 + (-t168 + t260) * V_base(5) + t263;
t236 = t216 * t168 + (-t288 + t294) * V_base(4) + t240;
t166 = t259 * t232;
t180 = pkin(4) * t215 + pkin(8) * t214;
t235 = t209 * t180 + (t166 + t258) * t216 + t241;
t165 = t259 * t230;
t234 = -t209 * t165 - t210 * t166 + t237;
t233 = t216 * t165 - t210 * t180 + t236;
t207 = rSges(2,1) * t232 - t230 * rSges(2,2);
t206 = -rSges(3,2) * t232 + t230 * rSges(3,3);
t204 = t230 * rSges(2,1) + rSges(2,2) * t232;
t203 = -t230 * rSges(3,2) - rSges(3,3) * t232;
t189 = rSges(4,1) * t227 - rSges(4,2) * t226;
t185 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t184 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t183 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t182 = qJD(5) * t214 + t216;
t179 = rSges(5,1) * t215 - rSges(5,2) * t214;
t170 = -t230 * t268 + t210;
t169 = t232 * t268 + t209;
t163 = t230 * rSges(4,3) - t232 * t257;
t162 = rSges(4,3) * t232 + t230 * t257;
t153 = t230 * rSges(5,3) - t232 * t256;
t152 = rSges(5,3) * t232 + t230 * t256;
t143 = rSges(6,3) * t214 + (rSges(6,1) * t231 - rSges(6,2) * t229) * t215;
t141 = V_base(5) * rSges(2,3) - t204 * t216 + t266;
t140 = t207 * t216 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t132 = t204 * V_base(4) - t207 * V_base(5) + V_base(3);
t129 = V_base(5) * rSges(3,1) + (-t202 - t203) * t216 + t262;
t128 = t216 * t206 + (-rSges(3,1) - pkin(6)) * V_base(4) + t243;
t127 = t174 * rSges(6,1) + t173 * rSges(6,2) + rSges(6,3) * t278;
t125 = rSges(6,1) * t172 - rSges(6,2) * t171 - rSges(6,3) * t279;
t111 = t203 * V_base(4) + (-t205 - t206) * V_base(5) + t267;
t110 = t189 * V_base(5) + (-t163 + t261) * t216 + t242;
t109 = t216 * t162 + (-t189 + t294) * V_base(4) + t240;
t108 = V_base(4) * t163 + (-t162 + t260) * V_base(5) + t263;
t107 = t179 * t209 + (-t153 + t258) * t216 + t241;
t106 = t216 * t152 - t210 * t179 + t236;
t105 = -t209 * t152 + t210 * t153 + t237;
t104 = -t127 * t182 + t143 * t169 + t235;
t103 = t182 * t125 - t170 * t143 + t233;
t102 = -t169 * t125 + t170 * t127 + t234;
t101 = qJD(6) * t171 + t169 * t269 - t182 * t270 + t235;
t100 = -qJD(6) * t173 - t170 * t269 + t182 * t271 + t233;
t99 = qJD(6) * t215 * t229 - t169 * t271 + t170 * t270 + t234;
t1 = m(1) * (t183 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(3) * (t111 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t210 * (t295 * t230 + t239 * t232) / 0.2e1 + t209 * (t239 * t230 - t295 * t232) / 0.2e1 + ((-t304 * t173 + t302 * t174 + t303 * t278) * t182 + (-t310 * t173 + t306 * t174 + t308 * t278) * t170 + (-t309 * t173 + t305 * t174 + t307 * t278) * t169) * t169 / 0.2e1 + ((t304 * t171 + t302 * t172 - t303 * t279) * t182 + (t310 * t171 + t306 * t172 - t308 * t279) * t170 + (t309 * t171 + t305 * t172 - t307 * t279) * t169) * t170 / 0.2e1 + (((t304 * t229 + t302 * t231) * t182 + (t310 * t229 + t306 * t231) * t170 + (t309 * t229 + t305 * t231) * t169) * t215 + (t307 * t169 + t308 * t170 + t303 * t182) * t214) * t182 / 0.2e1 + ((-t147 * t214 + t149 * t215) * t210 + (-t148 * t214 + t150 * t215) * t209 + (-t159 * t226 + t161 * t227 + t299) * V_base(5) + (-t158 * t226 + t160 * t227 + t298) * V_base(4) + (-t177 * t214 + t178 * t215 - t187 * t226 + t188 * t227 + Icges(3,1) + Icges(2,3)) * t216) * t216 / 0.2e1 + (t238 * t232 + (t230 * t244 + t298) * t216 + (t297 * t230 + t301 * t232 + Icges(1,4)) * V_base(5) + (t296 * t230 + t300 * t232 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t238 * t230 + (-t232 * t244 + t299) * t216 + (t301 * t230 - t297 * t232 + Icges(1,2)) * V_base(5) + (t300 * t230 - t296 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
