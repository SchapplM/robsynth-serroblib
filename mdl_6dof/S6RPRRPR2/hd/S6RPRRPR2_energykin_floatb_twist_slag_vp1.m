% Calculate kinetic energy for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:09
% EndTime: 2019-03-09 05:00:12
% DurationCPUTime: 3.58s
% Computational Cost: add. (2375->364), mult. (2092->519), div. (0->0), fcn. (1984->12), ass. (0->169)
t312 = -Icges(6,3) - Icges(5,3);
t246 = qJ(4) + pkin(11);
t235 = sin(t246);
t237 = cos(t246);
t247 = qJ(1) + pkin(10);
t238 = cos(t247);
t236 = sin(t247);
t253 = cos(qJ(3));
t290 = t236 * t253;
t168 = -t235 * t290 - t237 * t238;
t169 = -t235 * t238 + t237 * t290;
t252 = cos(qJ(4));
t249 = sin(qJ(4));
t286 = t249 * t253;
t180 = -t236 * t286 - t238 * t252;
t285 = t252 * t253;
t289 = t238 * t249;
t181 = t236 * t285 - t289;
t250 = sin(qJ(3));
t291 = t236 * t250;
t311 = Icges(5,5) * t181 + Icges(6,5) * t169 + Icges(5,6) * t180 + Icges(6,6) * t168 - t291 * t312;
t287 = t238 * t253;
t170 = -t235 * t287 + t236 * t237;
t171 = t235 * t236 + t237 * t287;
t182 = t236 * t252 - t238 * t286;
t292 = t236 * t249;
t183 = t238 * t285 + t292;
t288 = t238 * t250;
t310 = Icges(5,5) * t183 + Icges(6,5) * t171 + Icges(5,6) * t182 + Icges(6,6) * t170 - t288 * t312;
t309 = t312 * t253 + (Icges(5,5) * t252 + Icges(6,5) * t237 - Icges(5,6) * t249 - Icges(6,6) * t235) * t250;
t251 = sin(qJ(1));
t302 = pkin(1) * t251;
t254 = cos(qJ(1));
t301 = pkin(1) * t254;
t299 = t252 * pkin(4);
t298 = -pkin(6) - qJ(2);
t296 = Icges(2,4) * t251;
t295 = Icges(3,4) * t236;
t294 = Icges(4,4) * t250;
t293 = Icges(4,4) * t253;
t284 = pkin(5) * t237;
t282 = qJD(4) * t250;
t281 = qJD(5) * t250;
t280 = qJD(6) * t250;
t240 = V_base(6) + qJD(1);
t279 = t240 * t301 + V_base(2);
t278 = V_base(5) * pkin(6) + V_base(1);
t209 = qJD(3) * t236 + V_base(4);
t201 = pkin(2) * t236 - pkin(7) * t238;
t275 = -t201 - t302;
t274 = pkin(5) * t235;
t273 = V_base(5) * qJ(2) + t278;
t272 = t302 * V_base(4) + qJD(2) + V_base(3);
t179 = t238 * t282 + t209;
t271 = pkin(3) * t253 + pkin(8) * t250;
t208 = -qJD(3) * t238 + V_base(5);
t270 = rSges(4,1) * t253 - rSges(4,2) * t250;
t269 = Icges(4,1) * t253 - t294;
t268 = -Icges(4,2) * t250 + t293;
t267 = Icges(4,5) * t253 - Icges(4,6) * t250;
t178 = t236 * t282 + t208;
t266 = qJ(5) * t250 + t253 * t299;
t265 = (-Icges(4,3) * t238 + t236 * t267) * t208 + (Icges(4,3) * t236 + t238 * t267) * t209 + (Icges(4,5) * t250 + Icges(4,6) * t253) * t240;
t264 = pkin(9) * t250 + t253 * t284;
t202 = pkin(2) * t238 + pkin(7) * t236;
t263 = t202 * t240 + t298 * V_base(4) + t279;
t189 = t271 * t236;
t227 = pkin(3) * t250 - pkin(8) * t253;
t262 = t208 * t227 + (-t189 + t275) * t240 + t273;
t261 = V_base(4) * t201 + (-t202 - t301) * V_base(5) + t272;
t190 = t271 * t238;
t260 = t190 * t240 - t209 * t227 + t263;
t167 = -qJ(5) * t253 + t250 * t299;
t259 = t167 * t178 + t238 * t281 + t262;
t137 = pkin(4) * t292 + t238 * t266;
t223 = -qJD(4) * t253 + t240;
t258 = t137 * t223 + t236 * t281 + t260;
t257 = t189 * t209 - t190 * t208 + t261;
t136 = -pkin(4) * t289 + t236 * t266;
t256 = -qJD(5) * t253 + t136 * t179 + t257;
t158 = -Icges(4,6) * t238 + t236 * t268;
t159 = Icges(4,6) * t236 + t238 * t268;
t160 = -Icges(4,5) * t238 + t236 * t269;
t161 = Icges(4,5) * t236 + t238 * t269;
t213 = Icges(4,2) * t253 + t294;
t216 = Icges(4,1) * t250 + t293;
t255 = (-t159 * t250 + t161 * t253) * t209 + (-t158 * t250 + t160 * t253) * t208 + (-t213 * t250 + t216 * t253) * t240;
t242 = Icges(2,4) * t254;
t239 = qJ(6) + t246;
t233 = cos(t239);
t232 = sin(t239);
t231 = Icges(3,4) * t238;
t226 = rSges(2,1) * t254 - rSges(2,2) * t251;
t225 = rSges(2,1) * t251 + rSges(2,2) * t254;
t224 = rSges(4,1) * t250 + rSges(4,2) * t253;
t218 = Icges(2,1) * t254 - t296;
t217 = Icges(2,1) * t251 + t242;
t215 = -Icges(2,2) * t251 + t242;
t214 = Icges(2,2) * t254 + t296;
t206 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t205 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t204 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t200 = rSges(3,1) * t238 - rSges(3,2) * t236;
t199 = rSges(3,1) * t236 + rSges(3,2) * t238;
t198 = (-qJD(4) - qJD(6)) * t253 + t240;
t197 = Icges(3,1) * t238 - t295;
t196 = Icges(3,1) * t236 + t231;
t195 = -Icges(3,2) * t236 + t231;
t194 = Icges(3,2) * t238 + t295;
t187 = -rSges(5,3) * t253 + (rSges(5,1) * t252 - rSges(5,2) * t249) * t250;
t186 = -Icges(5,5) * t253 + (Icges(5,1) * t252 - Icges(5,4) * t249) * t250;
t185 = -Icges(5,6) * t253 + (Icges(5,4) * t252 - Icges(5,2) * t249) * t250;
t176 = -rSges(6,3) * t253 + (rSges(6,1) * t237 - rSges(6,2) * t235) * t250;
t175 = -Icges(6,5) * t253 + (Icges(6,1) * t237 - Icges(6,4) * t235) * t250;
t174 = -Icges(6,6) * t253 + (Icges(6,4) * t237 - Icges(6,2) * t235) * t250;
t166 = rSges(4,3) * t236 + t238 * t270;
t165 = -rSges(4,3) * t238 + t236 * t270;
t164 = -rSges(7,3) * t253 + (rSges(7,1) * t233 - rSges(7,2) * t232) * t250;
t163 = V_base(5) * rSges(2,3) - t225 * t240 + t278;
t162 = t226 * t240 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = -Icges(7,5) * t253 + (Icges(7,1) * t233 - Icges(7,4) * t232) * t250;
t154 = -Icges(7,6) * t253 + (Icges(7,4) * t233 - Icges(7,2) * t232) * t250;
t153 = -Icges(7,3) * t253 + (Icges(7,5) * t233 - Icges(7,6) * t232) * t250;
t152 = t232 * t236 + t233 * t287;
t151 = -t232 * t287 + t233 * t236;
t150 = -t232 * t238 + t233 * t290;
t149 = -t232 * t290 - t233 * t238;
t148 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t147 = t238 * t280 + t179;
t146 = t236 * t280 + t178;
t144 = -pkin(9) * t253 + t250 * t284;
t143 = V_base(5) * rSges(3,3) + (-t199 - t302) * t240 + t273;
t142 = t200 * t240 + (-rSges(3,3) + t298) * V_base(4) + t279;
t140 = V_base(4) * t199 + (-t200 - t301) * V_base(5) + t272;
t139 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t288;
t138 = rSges(5,1) * t181 + rSges(5,2) * t180 + rSges(5,3) * t291;
t135 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t288;
t134 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t291;
t133 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t288;
t132 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t291;
t129 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t288;
t128 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t291;
t127 = Icges(6,1) * t171 + Icges(6,4) * t170 + Icges(6,5) * t288;
t126 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t291;
t125 = Icges(6,4) * t171 + Icges(6,2) * t170 + Icges(6,6) * t288;
t124 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t291;
t120 = rSges(7,1) * t152 + rSges(7,2) * t151 + rSges(7,3) * t288;
t119 = rSges(7,1) * t150 + rSges(7,2) * t149 + rSges(7,3) * t291;
t118 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t288;
t117 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t291;
t116 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t288;
t115 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t291;
t114 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t288;
t113 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t291;
t111 = t236 * t274 + t238 * t264;
t110 = t236 * t264 - t238 * t274;
t109 = t208 * t224 + (-t165 + t275) * t240 + t273;
t108 = t166 * t240 - t209 * t224 + t263;
t107 = t165 * t209 - t166 * t208 + t261;
t106 = -t138 * t223 + t178 * t187 + t262;
t105 = t139 * t223 - t179 * t187 + t260;
t104 = t138 * t179 - t139 * t178 + t257;
t103 = t176 * t178 + (-t128 - t136) * t223 + t259;
t102 = t129 * t223 + (-t167 - t176) * t179 + t258;
t101 = t179 * t128 + (-t129 - t137) * t178 + t256;
t100 = -t119 * t198 + t144 * t178 + t146 * t164 + (-t110 - t136) * t223 + t259;
t99 = t111 * t223 + t120 * t198 - t147 * t164 + (-t144 - t167) * t179 + t258;
t98 = t179 * t110 + t147 * t119 - t146 * t120 + (-t111 - t137) * t178 + t256;
t1 = t147 * ((t114 * t288 + t151 * t116 + t152 * t118) * t147 + (t113 * t288 + t115 * t151 + t117 * t152) * t146 + (t151 * t154 + t152 * t155 + t153 * t288) * t198) / 0.2e1 + t146 * ((t114 * t291 + t116 * t149 + t118 * t150) * t147 + (t113 * t291 + t149 * t115 + t150 * t117) * t146 + (t149 * t154 + t150 * t155 + t153 * t291) * t198) / 0.2e1 + m(2) * (t148 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t140 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t209 * (t265 * t236 + t255 * t238) / 0.2e1 + t208 * (t255 * t236 - t265 * t238) / 0.2e1 + m(1) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + t198 * ((-t113 * t146 - t114 * t147 - t153 * t198) * t253 + ((-t116 * t232 + t118 * t233) * t147 + (-t115 * t232 + t117 * t233) * t146 + (-t154 * t232 + t155 * t233) * t198) * t250) / 0.2e1 + ((t168 * t174 + t169 * t175 + t180 * t185 + t181 * t186 + t291 * t309) * t223 + (t125 * t168 + t127 * t169 + t133 * t180 + t135 * t181 + t291 * t310) * t179 + (t168 * t124 + t169 * t126 + t180 * t132 + t181 * t134 + t311 * t291) * t178) * t178 / 0.2e1 + ((t170 * t174 + t171 * t175 + t182 * t185 + t183 * t186 + t288 * t309) * t223 + (t170 * t125 + t171 * t127 + t182 * t133 + t183 * t135 + t310 * t288) * t179 + (t124 * t170 + t126 * t171 + t132 * t182 + t134 * t183 + t288 * t311) * t178) * t179 / 0.2e1 + ((-t178 * t311 - t179 * t310 - t223 * t309) * t253 + ((-t174 * t235 + t175 * t237 - t185 * t249 + t186 * t252) * t223 + (-t125 * t235 + t127 * t237 - t133 * t249 + t135 * t252) * t179 + (-t124 * t235 + t126 * t237 - t132 * t249 + t134 * t252) * t178) * t250) * t223 / 0.2e1 + ((t159 * t253 + t161 * t250) * t209 + (t158 * t253 + t160 * t250) * t208 + (t253 * t213 + t250 * t216 + Icges(2,3) + Icges(3,3)) * t240) * t240 / 0.2e1 + ((-t194 * t236 + t196 * t238 - t214 * t251 + t217 * t254 + Icges(1,4)) * V_base(5) + (-t236 * t195 + t238 * t197 - t251 * t215 + t254 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t238 * t194 + t236 * t196 + t254 * t214 + t251 * t217 + Icges(1,2)) * V_base(5) + (t195 * t238 + t197 * t236 + t215 * t254 + t218 * t251 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t240 * (Icges(2,5) * t251 + Icges(3,5) * t236 + Icges(2,6) * t254 + Icges(3,6) * t238) + V_base(4) * t240 * (Icges(2,5) * t254 + Icges(3,5) * t238 - Icges(2,6) * t251 - Icges(3,6) * t236) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
