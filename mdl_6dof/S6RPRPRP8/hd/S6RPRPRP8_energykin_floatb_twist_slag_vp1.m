% Calculate kinetic energy for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:15
% EndTime: 2019-03-09 03:24:18
% DurationCPUTime: 2.75s
% Computational Cost: add. (1361->279), mult. (1724->372), div. (0->0), fcn. (1572->8), ass. (0->145)
t331 = Icges(2,4) + Icges(3,6);
t330 = Icges(2,1) + Icges(3,2);
t329 = Icges(6,1) + Icges(7,1);
t328 = -Icges(3,4) + Icges(2,5);
t327 = Icges(6,4) - Icges(7,5);
t326 = Icges(7,4) + Icges(6,5);
t325 = Icges(3,5) - Icges(2,6);
t324 = Icges(2,2) + Icges(3,3);
t323 = Icges(6,2) + Icges(7,3);
t322 = Icges(7,6) - Icges(6,6);
t321 = Icges(4,3) + Icges(5,3);
t320 = Icges(6,3) + Icges(7,2);
t225 = qJ(3) + pkin(9);
t214 = sin(t225);
t215 = cos(t225);
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t319 = Icges(4,5) * t228 + Icges(5,5) * t214 + Icges(4,6) * t231 + Icges(5,6) * t215;
t318 = rSges(7,1) + pkin(5);
t317 = rSges(7,3) + qJ(6);
t232 = cos(qJ(1));
t316 = t331 * t232;
t229 = sin(qJ(1));
t315 = t331 * t229;
t280 = Icges(5,4) * t214;
t252 = Icges(5,2) * t215 + t280;
t148 = Icges(5,6) * t232 + t229 * t252;
t149 = Icges(5,6) * t229 - t232 * t252;
t279 = Icges(5,4) * t215;
t254 = Icges(5,1) * t214 + t279;
t150 = Icges(5,5) * t232 + t229 * t254;
t151 = Icges(5,5) * t229 - t232 * t254;
t282 = Icges(4,4) * t228;
t253 = Icges(4,2) * t231 + t282;
t158 = Icges(4,6) * t232 + t229 * t253;
t159 = Icges(4,6) * t229 - t232 * t253;
t281 = Icges(4,4) * t231;
t255 = Icges(4,1) * t228 + t281;
t160 = Icges(4,5) * t232 + t229 * t255;
t161 = Icges(4,5) * t229 - t232 * t255;
t177 = -Icges(5,2) * t214 + t279;
t178 = Icges(5,1) * t215 - t280;
t194 = -Icges(4,2) * t228 + t281;
t199 = Icges(4,1) * t231 - t282;
t210 = qJD(3) * t229 + V_base(5);
t211 = qJD(3) * t232 + V_base(4);
t216 = V_base(6) + qJD(1);
t314 = t210 * (t149 * t215 + t151 * t214 + t159 * t231 + t161 * t228) + t211 * (t148 * t215 + t150 * t214 + t158 * t231 + t160 * t228) + t216 * (t177 * t215 + t178 * t214 + t194 * t231 + t199 * t228);
t230 = cos(qJ(5));
t271 = t232 * t230;
t227 = sin(qJ(5));
t273 = t229 * t227;
t171 = t214 * t273 - t271;
t272 = t229 * t230;
t274 = t227 * t232;
t172 = t214 * t272 + t274;
t276 = t215 * t229;
t313 = t323 * t171 - t327 * t172 - t322 * t276;
t173 = t214 * t274 + t272;
t174 = -t214 * t271 + t273;
t275 = t215 * t232;
t312 = -t323 * t173 - t327 * t174 + t322 * t275;
t311 = t322 * t171 + t326 * t172 - t320 * t276;
t310 = -t322 * t173 + t326 * t174 + t320 * t275;
t309 = -t327 * t171 + t329 * t172 - t326 * t276;
t308 = t327 * t173 + t329 * t174 + t326 * t275;
t307 = (t323 * t227 - t327 * t230) * t215 + t322 * t214;
t306 = (t322 * t227 + t326 * t230) * t215 + t320 * t214;
t305 = (-t327 * t227 + t329 * t230) * t215 + t326 * t214;
t304 = -t324 * t232 - t315;
t303 = t324 * t229 - t316;
t302 = t330 * t229 + t316;
t301 = t330 * t232 - t315;
t298 = (Icges(4,5) * t231 + Icges(5,5) * t215 - Icges(4,6) * t228 - Icges(5,6) * t214) * t216 + (t319 * t229 + t321 * t232) * t211 + (t321 * t229 - t319 * t232) * t210;
t288 = pkin(3) * t228;
t287 = pkin(3) * t231;
t286 = pkin(7) * t232;
t285 = t229 * pkin(7);
t270 = -rSges(7,2) * t276 + t317 * t171 + t318 * t172;
t269 = rSges(7,2) * t275 - t317 * t173 + t318 * t174;
t268 = rSges(7,2) * t214 + (t317 * t227 + t318 * t230) * t215;
t267 = qJD(5) * t215;
t202 = t229 * pkin(1) - qJ(2) * t232;
t266 = V_base(4) * t202 + V_base(3);
t265 = V_base(5) * pkin(6) + V_base(1);
t262 = -t202 - t285;
t261 = qJD(2) * t229 + t265;
t167 = qJ(4) * t229 - t232 * t288;
t260 = -t167 + t262;
t259 = V_base(5) * pkin(2) + t261;
t258 = pkin(4) * t214 - pkin(8) * t215;
t257 = rSges(4,1) * t228 + rSges(4,2) * t231;
t256 = rSges(5,1) * t214 + rSges(5,2) * t215;
t206 = pkin(1) * t232 + t229 * qJ(2);
t243 = -qJD(2) * t232 + t216 * t206 + V_base(2);
t242 = qJD(4) * t232 + t210 * t287 + t259;
t239 = V_base(4) * t285 + (-t206 - t286) * V_base(5) + t266;
t238 = t216 * t286 + (-pkin(2) - pkin(6)) * V_base(4) + t243;
t237 = t211 * t167 + t239;
t168 = qJ(4) * t232 + t229 * t288;
t236 = qJD(4) * t229 + t216 * t168 + t238;
t166 = t258 * t232;
t180 = pkin(4) * t215 + pkin(8) * t214;
t235 = t210 * t180 + (t166 + t260) * t216 + t242;
t165 = t258 * t229;
t234 = -t211 * t166 + (-t165 - t168) * t210 + t237;
t233 = t216 * t165 + (-t180 - t287) * t211 + t236;
t208 = rSges(2,1) * t232 - t229 * rSges(2,2);
t207 = -rSges(3,2) * t232 + t229 * rSges(3,3);
t205 = rSges(4,1) * t231 - rSges(4,2) * t228;
t204 = t229 * rSges(2,1) + rSges(2,2) * t232;
t203 = -t229 * rSges(3,2) - rSges(3,3) * t232;
t186 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t185 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t184 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t183 = qJD(5) * t214 + t216;
t179 = rSges(5,1) * t215 - rSges(5,2) * t214;
t170 = -t229 * t267 + t211;
t169 = t232 * t267 + t210;
t164 = t229 * rSges(4,3) - t232 * t257;
t163 = rSges(4,3) * t232 + t229 * t257;
t153 = t229 * rSges(5,3) - t232 * t256;
t152 = rSges(5,3) * t232 + t229 * t256;
t144 = rSges(6,3) * t214 + (rSges(6,1) * t230 - rSges(6,2) * t227) * t215;
t142 = V_base(5) * rSges(2,3) - t204 * t216 + t265;
t141 = t208 * t216 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t132 = t204 * V_base(4) - t208 * V_base(5) + V_base(3);
t129 = V_base(5) * rSges(3,1) + (-t202 - t203) * t216 + t261;
t128 = t216 * t207 + (-rSges(3,1) - pkin(6)) * V_base(4) + t243;
t127 = t174 * rSges(6,1) + t173 * rSges(6,2) + rSges(6,3) * t275;
t125 = rSges(6,1) * t172 - rSges(6,2) * t171 - rSges(6,3) * t276;
t111 = t203 * V_base(4) + (-t206 - t207) * V_base(5) + t266;
t110 = t205 * t210 + (-t164 + t262) * t216 + t259;
t109 = t216 * t163 - t211 * t205 + t238;
t108 = -t210 * t163 + t211 * t164 + t239;
t107 = t179 * t210 + (-t153 + t260) * t216 + t242;
t106 = t216 * t152 + (-t179 - t287) * t211 + t236;
t105 = t211 * t153 + (-t152 - t168) * t210 + t237;
t104 = -t127 * t183 + t144 * t169 + t235;
t103 = t183 * t125 - t170 * t144 + t233;
t102 = -t169 * t125 + t170 * t127 + t234;
t101 = qJD(6) * t171 + t169 * t268 - t183 * t269 + t235;
t100 = -qJD(6) * t173 - t170 * t268 + t183 * t270 + t233;
t99 = qJD(6) * t215 * t227 - t169 * t270 + t170 * t269 + t234;
t1 = m(1) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(2) * (t132 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t111 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + ((-t173 * t307 + t174 * t305 + t275 * t306) * t183 + (-t173 * t313 + t309 * t174 + t311 * t275) * t170 + (-t312 * t173 + t308 * t174 + t310 * t275) * t169) * t169 / 0.2e1 + ((t171 * t307 + t172 * t305 - t276 * t306) * t183 + (t313 * t171 + t309 * t172 - t311 * t276) * t170 + (t171 * t312 + t172 * t308 - t310 * t276) * t169) * t170 / 0.2e1 + (((t227 * t307 + t230 * t305) * t183 + (t227 * t313 + t309 * t230) * t170 + (t227 * t312 + t230 * t308) * t169) * t215 + (t169 * t310 + t170 * t311 + t306 * t183) * t214) * t183 / 0.2e1 + (t298 * t229 - t314 * t232) * t210 / 0.2e1 + (t314 * t229 + t298 * t232) * t211 / 0.2e1 + ((t229 * t304 + t232 * t302 + Icges(1,4)) * V_base(5) + (t303 * t229 + t301 * t232 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t302 * t229 - t304 * t232 + Icges(1,2)) * V_base(5) + (t229 * t301 - t232 * t303 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t148 * t214 + t150 * t215 - t158 * t228 + t160 * t231) * t211 + (-t149 * t214 + t151 * t215 - t159 * t228 + t161 * t231) * t210 + (-t177 * t214 + t178 * t215 - t194 * t228 + t199 * t231 + Icges(3,1) + Icges(2,3)) * t216) * t216 / 0.2e1 + t216 * V_base(5) * (t328 * t229 - t325 * t232) + t216 * V_base(4) * (t325 * t229 + t328 * t232) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
