% Calculate kinetic energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:11:58
% EndTime: 2019-03-09 10:12:02
% DurationCPUTime: 3.63s
% Computational Cost: add. (1875->313), mult. (1756->438), div. (0->0), fcn. (1538->10), ass. (0->165)
t336 = Icges(5,4) + Icges(6,6);
t335 = Icges(5,1) + Icges(6,2);
t334 = -Icges(5,2) - Icges(6,3);
t235 = qJ(2) + pkin(10);
t226 = qJ(4) + t235;
t222 = cos(t226);
t333 = t336 * t222;
t221 = sin(t226);
t332 = t336 * t221;
t331 = Icges(6,4) - Icges(5,5);
t330 = Icges(6,5) - Icges(5,6);
t329 = t334 * t221 + t333;
t328 = t335 * t222 - t332;
t327 = Icges(6,1) + Icges(5,3);
t326 = Icges(3,3) + Icges(4,3);
t239 = sin(qJ(1));
t242 = cos(qJ(1));
t325 = t329 * t239 + t330 * t242;
t324 = -t330 * t239 + t329 * t242;
t323 = t328 * t239 + t331 * t242;
t322 = -t331 * t239 + t328 * t242;
t321 = t334 * t222 - t332;
t320 = t335 * t221 + t333;
t319 = t330 * t221 - t331 * t222;
t224 = sin(t235);
t225 = cos(t235);
t238 = sin(qJ(2));
t241 = cos(qJ(2));
t318 = Icges(3,5) * t241 + Icges(4,5) * t225 - Icges(3,6) * t238 - Icges(4,6) * t224;
t299 = Icges(4,4) * t225;
t268 = -Icges(4,2) * t224 + t299;
t155 = -Icges(4,6) * t242 + t239 * t268;
t156 = Icges(4,6) * t239 + t242 * t268;
t300 = Icges(4,4) * t224;
t271 = Icges(4,1) * t225 - t300;
t157 = -Icges(4,5) * t242 + t239 * t271;
t158 = Icges(4,5) * t239 + t242 * t271;
t301 = Icges(3,4) * t241;
t269 = -Icges(3,2) * t238 + t301;
t167 = -Icges(3,6) * t242 + t239 * t269;
t168 = Icges(3,6) * t239 + t242 * t269;
t302 = Icges(3,4) * t238;
t272 = Icges(3,1) * t241 - t302;
t169 = -Icges(3,5) * t242 + t239 * t272;
t170 = Icges(3,5) * t239 + t242 * t272;
t191 = Icges(4,2) * t225 + t300;
t192 = Icges(4,1) * t224 + t299;
t208 = Icges(3,2) * t241 + t302;
t211 = Icges(3,1) * t238 + t301;
t219 = -qJD(2) * t242 + V_base(5);
t220 = qJD(2) * t239 + V_base(4);
t227 = V_base(6) + qJD(1);
t317 = (-t191 * t224 + t192 * t225 - t208 * t238 + t211 * t241) * t227 + (-t156 * t224 + t158 * t225 - t168 * t238 + t170 * t241) * t220 + (-t155 * t224 + t157 * t225 - t167 * t238 + t169 * t241) * t219;
t195 = V_base(5) + (-qJD(2) - qJD(4)) * t242;
t196 = qJD(4) * t239 + t220;
t316 = (t221 * t321 + t222 * t320) * t227 + (-t221 * t324 + t222 * t322) * t196 + (-t221 * t325 + t222 * t323) * t195;
t315 = (Icges(3,5) * t238 + Icges(4,5) * t224 + Icges(3,6) * t241 + Icges(4,6) * t225) * t227 + (t239 * t326 + t242 * t318) * t220 + (t239 * t318 - t242 * t326) * t219;
t314 = (-t331 * t221 - t330 * t222) * t227 + (t239 * t327 + t319 * t242) * t196 + (t319 * t239 - t242 * t327) * t195;
t308 = pkin(2) * t238;
t307 = pkin(3) * t224;
t306 = pkin(9) * t221;
t305 = t241 * pkin(2);
t303 = Icges(2,4) * t239;
t294 = t222 * t239;
t293 = t222 * t242;
t237 = sin(qJ(6));
t292 = t237 * t242;
t291 = t239 * t237;
t240 = cos(qJ(6));
t290 = t239 * t240;
t289 = t240 * t242;
t151 = -qJ(3) * t242 + t239 * t305;
t217 = t239 * pkin(1) - pkin(7) * t242;
t288 = -t151 - t217;
t287 = pkin(3) * t225;
t285 = qJD(5) * t221;
t284 = qJD(6) * t222;
t283 = V_base(5) * pkin(6) + V_base(1);
t121 = -pkin(8) * t242 + t239 * t287;
t280 = -t121 + t288;
t273 = pkin(4) * t222 + qJ(5) * t221;
t163 = t273 * t239;
t279 = -t163 + t280;
t278 = qJD(3) * t239 + t219 * t308 + t283;
t277 = rSges(3,1) * t241 - rSges(3,2) * t238;
t276 = rSges(4,1) * t225 - rSges(4,2) * t224;
t275 = rSges(5,1) * t222 - rSges(5,2) * t221;
t274 = -rSges(6,2) * t222 + rSges(6,3) * t221;
t260 = t219 * t307 + t278;
t218 = pkin(1) * t242 + t239 * pkin(7);
t259 = -V_base(4) * pkin(6) + t227 * t218 + V_base(2);
t258 = V_base(4) * t217 - t218 * V_base(5) + V_base(3);
t257 = t220 * t151 + t258;
t183 = pkin(4) * t221 - qJ(5) * t222;
t256 = t195 * t183 + t242 * t285 + t260;
t152 = qJ(3) * t239 + t242 * t305;
t251 = -qJD(3) * t242 + t227 * t152 + t259;
t122 = pkin(8) * t239 + t242 * t287;
t250 = t220 * t121 + (-t122 - t152) * t219 + t257;
t249 = -qJD(5) * t222 + t196 * t163 + t250;
t248 = t227 * t122 + (-t307 - t308) * t220 + t251;
t164 = t273 * t242;
t247 = t227 * t164 + t239 * t285 + t248;
t231 = Icges(2,4) * t242;
t216 = rSges(2,1) * t242 - t239 * rSges(2,2);
t215 = t239 * rSges(2,1) + rSges(2,2) * t242;
t214 = rSges(3,1) * t238 + rSges(3,2) * t241;
t213 = Icges(2,1) * t242 - t303;
t212 = Icges(2,1) * t239 + t231;
t210 = -Icges(2,2) * t239 + t231;
t209 = Icges(2,2) * t242 + t303;
t202 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t201 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t200 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t194 = qJD(6) * t221 + t227;
t193 = rSges(4,1) * t224 + rSges(4,2) * t225;
t187 = -pkin(5) * t242 + pkin(9) * t294;
t186 = t239 * pkin(5) + pkin(9) * t293;
t185 = rSges(5,1) * t221 + rSges(5,2) * t222;
t184 = -rSges(6,2) * t221 - rSges(6,3) * t222;
t176 = t221 * t291 - t289;
t175 = t221 * t290 + t292;
t174 = t221 * t292 + t290;
t173 = t221 * t289 - t291;
t172 = t239 * rSges(3,3) + t242 * t277;
t171 = -rSges(3,3) * t242 + t239 * t277;
t162 = t239 * rSges(4,3) + t242 * t276;
t161 = -rSges(4,3) * t242 + t239 * t276;
t160 = t242 * t284 + t196;
t159 = t239 * t284 + t195;
t150 = -rSges(6,1) * t242 + t239 * t274;
t149 = t239 * rSges(6,1) + t242 * t274;
t148 = t239 * rSges(5,3) + t242 * t275;
t147 = -rSges(5,3) * t242 + t239 * t275;
t146 = V_base(5) * rSges(2,3) - t215 * t227 + t283;
t145 = t216 * t227 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t131 = t215 * V_base(4) - t216 * V_base(5) + V_base(3);
t130 = rSges(7,3) * t221 + (-rSges(7,1) * t237 - rSges(7,2) * t240) * t222;
t127 = Icges(7,5) * t221 + (-Icges(7,1) * t237 - Icges(7,4) * t240) * t222;
t126 = Icges(7,6) * t221 + (-Icges(7,4) * t237 - Icges(7,2) * t240) * t222;
t125 = Icges(7,3) * t221 + (-Icges(7,5) * t237 - Icges(7,6) * t240) * t222;
t118 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t294;
t117 = t174 * rSges(7,1) + t173 * rSges(7,2) + rSges(7,3) * t293;
t116 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t294;
t115 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t293;
t114 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t294;
t113 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t293;
t112 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t294;
t111 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t293;
t110 = t214 * t219 + (-t171 - t217) * t227 + t283;
t109 = t172 * t227 - t214 * t220 + t259;
t108 = t171 * t220 - t172 * t219 + t258;
t107 = t193 * t219 + (-t161 + t288) * t227 + t278;
t106 = t227 * t162 + (-t193 - t308) * t220 + t251;
t105 = t161 * t220 + (-t152 - t162) * t219 + t257;
t104 = t185 * t195 + (-t147 + t280) * t227 + t260;
t103 = t227 * t148 - t196 * t185 + t248;
t102 = t184 * t195 + (-t150 + t279) * t227 + t256;
t101 = t227 * t149 + (-t183 - t184) * t196 + t247;
t100 = t147 * t196 - t148 * t195 + t250;
t99 = t150 * t196 + (-t149 - t164) * t195 + t249;
t98 = t195 * t306 - t118 * t194 + t130 * t159 + (-t187 + t279) * t227 + t256;
t97 = t194 * t117 - t160 * t130 + t227 * t186 + (-t183 - t306) * t196 + t247;
t96 = -t117 * t159 + t118 * t160 + t187 * t196 + (-t164 - t186) * t195 + t249;
t1 = m(1) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(2) * (t131 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + t194 * ((t111 * t160 + t112 * t159 + t125 * t194) * t221 + ((-t113 * t240 - t115 * t237) * t160 + (-t114 * t240 - t116 * t237) * t159 + (-t126 * t240 - t127 * t237) * t194) * t222) / 0.2e1 + t160 * ((t111 * t293 + t173 * t113 + t174 * t115) * t160 + (t112 * t293 + t173 * t114 + t174 * t116) * t159 + (t125 * t293 + t173 * t126 + t174 * t127) * t194) / 0.2e1 + t159 * ((t111 * t294 + t113 * t175 + t115 * t176) * t160 + (t112 * t294 + t175 * t114 + t176 * t116) * t159 + (t125 * t294 + t126 * t175 + t127 * t176) * t194) / 0.2e1 + (t239 * t316 - t314 * t242) * t195 / 0.2e1 + (t239 * t314 + t316 * t242) * t196 / 0.2e1 + (t239 * t317 - t315 * t242) * t219 / 0.2e1 + (t315 * t239 + t317 * t242) * t220 / 0.2e1 + ((-t239 * t209 + t212 * t242 + Icges(1,4)) * V_base(5) + (-t239 * t210 + t213 * t242 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t209 * t242 + t239 * t212 + Icges(1,2)) * V_base(5) + (t210 * t242 + t239 * t213 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t156 * t225 + t158 * t224 + t168 * t241 + t170 * t238) * t220 + (t155 * t225 + t157 * t224 + t167 * t241 + t169 * t238) * t219 + (t221 * t322 + t222 * t324) * t196 + (t221 * t323 + t222 * t325) * t195 + (t191 * t225 + t192 * t224 + t208 * t241 + t211 * t238 + t221 * t320 - t321 * t222 + Icges(2,3)) * t227) * t227 / 0.2e1 + t227 * V_base(4) * (Icges(2,5) * t242 - Icges(2,6) * t239) + t227 * V_base(5) * (Icges(2,5) * t239 + Icges(2,6) * t242) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
