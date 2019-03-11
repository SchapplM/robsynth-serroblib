% Calculate kinetic energy for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:36
% EndTime: 2019-03-09 13:55:40
% DurationCPUTime: 4.88s
% Computational Cost: add. (1678->353), mult. (3174->517), div. (0->0), fcn. (3450->10), ass. (0->164)
t331 = Icges(3,4) - Icges(4,5);
t330 = Icges(3,1) + Icges(4,1);
t329 = Icges(3,2) + Icges(4,3);
t263 = cos(qJ(2));
t328 = t331 * t263;
t260 = sin(qJ(2));
t327 = t331 * t260;
t326 = Icges(4,4) + Icges(3,5);
t325 = Icges(3,6) - Icges(4,6);
t324 = t329 * t260 - t328;
t323 = t330 * t263 - t327;
t322 = Icges(4,2) + Icges(3,3);
t261 = sin(qJ(1));
t264 = cos(qJ(1));
t321 = t324 * t261 + t325 * t264;
t320 = -t325 * t261 + t324 * t264;
t319 = t323 * t261 - t326 * t264;
t318 = t326 * t261 + t323 * t264;
t317 = -t329 * t263 - t327;
t316 = t330 * t260 + t328;
t315 = -t325 * t260 + t326 * t263;
t245 = -qJD(2) * t264 + V_base(5);
t246 = qJD(2) * t261 + V_base(4);
t250 = V_base(6) + qJD(1);
t314 = (t317 * t260 + t316 * t263) * t250 + (t320 * t260 + t318 * t263) * t246 + (t321 * t260 + t319 * t263) * t245;
t313 = (t326 * t260 + t325 * t263) * t250 + (t322 * t261 + t315 * t264) * t246 + (t315 * t261 - t322 * t264) * t245;
t259 = sin(qJ(4));
t309 = cos(qJ(4));
t290 = t260 * t309;
t216 = -t263 * t259 + t290;
t308 = pkin(3) * t260;
t262 = cos(qJ(5));
t307 = pkin(5) * t262;
t305 = Icges(2,4) * t261;
t258 = sin(qJ(5));
t300 = t258 * t261;
t299 = t258 * t264;
t297 = t263 * t264;
t286 = pkin(2) * t263 + qJ(3) * t260;
t209 = t286 * t261;
t241 = pkin(1) * t261 - pkin(7) * t264;
t296 = -t209 - t241;
t295 = qJD(3) * t260;
t294 = V_base(5) * pkin(6) + V_base(1);
t219 = pkin(3) * t261 * t263 + pkin(8) * t264;
t291 = -t219 + t296;
t215 = t260 * t259 + t263 * t309;
t198 = qJD(5) * t215 + t250;
t236 = pkin(2) * t260 - qJ(3) * t263;
t289 = t245 * t236 + t264 * t295 + t294;
t288 = rSges(3,1) * t263 - rSges(3,2) * t260;
t287 = rSges(4,1) * t263 + rSges(4,3) * t260;
t279 = t245 * t308 + t289;
t213 = qJD(4) * t264 + t245;
t214 = -qJD(4) * t261 + t246;
t242 = pkin(1) * t264 + pkin(7) * t261;
t278 = -V_base(4) * pkin(6) + t250 * t242 + V_base(2);
t277 = V_base(4) * t241 - t242 * V_base(5) + V_base(3);
t204 = t216 * t261;
t163 = -qJD(5) * t204 + t213;
t206 = t259 * t297 - t264 * t290;
t164 = qJD(5) * t206 + t214;
t210 = t286 * t264;
t274 = t250 * t210 + t261 * t295 + t278;
t273 = -qJD(3) * t263 + t246 * t209 + t277;
t205 = t215 * t261;
t160 = t205 * pkin(4) - t204 * pkin(9);
t173 = t216 * pkin(4) + t215 * pkin(9);
t272 = t213 * t173 + (-t160 + t291) * t250 + t279;
t220 = pkin(3) * t297 - pkin(8) * t261;
t271 = t250 * t220 + (-t236 - t308) * t246 + t274;
t270 = t246 * t219 + (-t210 - t220) * t245 + t273;
t207 = t215 * t264;
t161 = t207 * pkin(4) + t206 * pkin(9);
t269 = t250 * t161 - t173 * t214 + t271;
t268 = t214 * t160 - t161 * t213 + t270;
t257 = qJ(5) + qJ(6);
t255 = Icges(2,4) * t264;
t254 = cos(t257);
t253 = sin(t257);
t240 = rSges(2,1) * t264 - rSges(2,2) * t261;
t239 = rSges(2,1) * t261 + rSges(2,2) * t264;
t238 = rSges(3,1) * t260 + rSges(3,2) * t263;
t237 = rSges(4,1) * t260 - rSges(4,3) * t263;
t235 = Icges(2,1) * t264 - t305;
t234 = Icges(2,1) * t261 + t255;
t231 = -Icges(2,2) * t261 + t255;
t230 = Icges(2,2) * t264 + t305;
t223 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t222 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t221 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t202 = rSges(3,3) * t261 + t264 * t288;
t201 = rSges(4,2) * t261 + t264 * t287;
t200 = -rSges(3,3) * t264 + t261 * t288;
t199 = -rSges(4,2) * t264 + t261 * t287;
t180 = V_base(5) * rSges(2,3) - t239 * t250 + t294;
t179 = t240 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t178 = t207 * t262 - t300;
t177 = -t207 * t258 - t261 * t262;
t176 = t205 * t262 + t299;
t175 = -t205 * t258 + t262 * t264;
t174 = t239 * V_base(4) - t240 * V_base(5) + V_base(3);
t172 = rSges(5,1) * t216 - rSges(5,2) * t215;
t171 = Icges(5,1) * t216 - Icges(5,4) * t215;
t170 = Icges(5,4) * t216 - Icges(5,2) * t215;
t169 = Icges(5,5) * t216 - Icges(5,6) * t215;
t168 = t207 * t254 - t253 * t261;
t167 = -t207 * t253 - t254 * t261;
t166 = t205 * t254 + t253 * t264;
t165 = -t205 * t253 + t254 * t264;
t162 = qJD(6) * t215 + t198;
t158 = rSges(5,1) * t207 - rSges(5,2) * t206 - rSges(5,3) * t261;
t157 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t264;
t156 = Icges(5,1) * t207 - Icges(5,4) * t206 - Icges(5,5) * t261;
t155 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t264;
t154 = Icges(5,4) * t207 - Icges(5,2) * t206 - Icges(5,6) * t261;
t153 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t264;
t152 = Icges(5,5) * t207 - Icges(5,6) * t206 - Icges(5,3) * t261;
t151 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t264;
t149 = rSges(6,3) * t215 + (rSges(6,1) * t262 - rSges(6,2) * t258) * t216;
t148 = Icges(6,5) * t215 + (Icges(6,1) * t262 - Icges(6,4) * t258) * t216;
t147 = Icges(6,6) * t215 + (Icges(6,4) * t262 - Icges(6,2) * t258) * t216;
t146 = Icges(6,3) * t215 + (Icges(6,5) * t262 - Icges(6,6) * t258) * t216;
t145 = qJD(6) * t206 + t164;
t144 = -qJD(6) * t204 + t163;
t143 = rSges(7,3) * t215 + (rSges(7,1) * t254 - rSges(7,2) * t253) * t216;
t141 = Icges(7,5) * t215 + (Icges(7,1) * t254 - Icges(7,4) * t253) * t216;
t140 = Icges(7,6) * t215 + (Icges(7,4) * t254 - Icges(7,2) * t253) * t216;
t139 = Icges(7,3) * t215 + (Icges(7,5) * t254 - Icges(7,6) * t253) * t216;
t138 = pkin(10) * t215 + t216 * t307;
t137 = t238 * t245 + (-t200 - t241) * t250 + t294;
t136 = t202 * t250 - t238 * t246 + t278;
t135 = rSges(6,1) * t178 + rSges(6,2) * t177 + rSges(6,3) * t206;
t134 = rSges(6,1) * t176 + rSges(6,2) * t175 - rSges(6,3) * t204;
t133 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t206;
t132 = Icges(6,1) * t176 + Icges(6,4) * t175 - Icges(6,5) * t204;
t131 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t206;
t130 = Icges(6,4) * t176 + Icges(6,2) * t175 - Icges(6,6) * t204;
t129 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t206;
t128 = Icges(6,5) * t176 + Icges(6,6) * t175 - Icges(6,3) * t204;
t127 = t200 * t246 - t202 * t245 + t277;
t126 = rSges(7,1) * t168 + rSges(7,2) * t167 + rSges(7,3) * t206;
t125 = rSges(7,1) * t166 + rSges(7,2) * t165 - rSges(7,3) * t204;
t124 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t206;
t123 = Icges(7,1) * t166 + Icges(7,4) * t165 - Icges(7,5) * t204;
t122 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t206;
t121 = Icges(7,4) * t166 + Icges(7,2) * t165 - Icges(7,6) * t204;
t120 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t206;
t119 = Icges(7,5) * t166 + Icges(7,6) * t165 - Icges(7,3) * t204;
t118 = -pkin(5) * t300 + pkin(10) * t206 + t207 * t307;
t117 = pkin(5) * t299 - pkin(10) * t204 + t205 * t307;
t116 = t237 * t245 + (-t199 + t296) * t250 + t289;
t115 = t201 * t250 + (-t236 - t237) * t246 + t274;
t114 = t199 * t246 + (-t201 - t210) * t245 + t273;
t113 = t172 * t213 + (-t157 + t291) * t250 + t279;
t112 = t158 * t250 - t172 * t214 + t271;
t111 = t157 * t214 - t158 * t213 + t270;
t110 = -t134 * t198 + t149 * t163 + t272;
t109 = t135 * t198 - t149 * t164 + t269;
t108 = t134 * t164 - t135 * t163 + t268;
t107 = -t117 * t198 - t125 * t162 + t138 * t163 + t143 * t144 + t272;
t106 = t118 * t198 + t126 * t162 - t138 * t164 - t143 * t145 + t269;
t105 = t117 * t164 - t118 * t163 + t125 * t145 - t126 * t144 + t268;
t1 = t198 * ((t128 * t163 + t129 * t164 + t146 * t198) * t215 + ((-t131 * t258 + t133 * t262) * t164 + (-t130 * t258 + t132 * t262) * t163 + (-t147 * t258 + t148 * t262) * t198) * t216) / 0.2e1 + t162 * ((t119 * t144 + t120 * t145 + t139 * t162) * t215 + ((-t122 * t253 + t124 * t254) * t145 + (-t121 * t253 + t123 * t254) * t144 + (-t140 * t253 + t141 * t254) * t162) * t216) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(2) * (t174 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + t144 * ((-t120 * t204 + t122 * t165 + t124 * t166) * t145 + (-t204 * t119 + t165 * t121 + t166 * t123) * t144 + (-t139 * t204 + t140 * t165 + t141 * t166) * t162) / 0.2e1 + t163 * ((-t129 * t204 + t131 * t175 + t133 * t176) * t164 + (-t204 * t128 + t175 * t130 + t176 * t132) * t163 + (-t146 * t204 + t147 * t175 + t148 * t176) * t198) / 0.2e1 + t145 * ((t206 * t120 + t167 * t122 + t168 * t124) * t145 + (t119 * t206 + t121 * t167 + t123 * t168) * t144 + (t139 * t206 + t140 * t167 + t141 * t168) * t162) / 0.2e1 + t164 * ((t206 * t129 + t177 * t131 + t178 * t133) * t164 + (t128 * t206 + t130 * t177 + t132 * t178) * t163 + (t146 * t206 + t147 * t177 + t148 * t178) * t198) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + t214 * ((-t261 * t152 - t206 * t154 + t207 * t156) * t214 + (-t151 * t261 - t153 * t206 + t155 * t207) * t213 + (-t169 * t261 - t170 * t206 + t171 * t207) * t250) / 0.2e1 + t213 * ((t152 * t264 + t154 * t204 + t156 * t205) * t214 + (t264 * t151 + t204 * t153 + t205 * t155) * t213 + (t169 * t264 + t170 * t204 + t171 * t205) * t250) / 0.2e1 + (t314 * t261 - t313 * t264) * t245 / 0.2e1 + (t313 * t261 + t314 * t264) * t246 / 0.2e1 + ((-t230 * t261 + t234 * t264 + Icges(1,4)) * V_base(5) + (-t261 * t231 + t264 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t264 * t230 + t261 * t234 + Icges(1,2)) * V_base(5) + (t231 * t264 + t235 * t261 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t154 * t215 + t156 * t216) * t214 + (-t153 * t215 + t155 * t216) * t213 + (t318 * t260 - t320 * t263) * t246 + (t319 * t260 - t321 * t263) * t245 + (-t215 * t170 + t216 * t171 + t316 * t260 - t317 * t263 + Icges(2,3)) * t250) * t250 / 0.2e1 + t250 * V_base(4) * (Icges(2,5) * t264 - Icges(2,6) * t261) + t250 * V_base(5) * (Icges(2,5) * t261 + Icges(2,6) * t264) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
