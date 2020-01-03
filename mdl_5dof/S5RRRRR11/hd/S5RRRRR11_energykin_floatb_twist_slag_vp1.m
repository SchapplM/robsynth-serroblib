% Calculate kinetic energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:43
% EndTime: 2019-12-31 22:38:46
% DurationCPUTime: 3.01s
% Computational Cost: add. (2253->360), mult. (4776->539), div. (0->0), fcn. (5730->12), ass. (0->158)
t294 = cos(qJ(3));
t257 = cos(pkin(5));
t293 = pkin(7) * t257;
t262 = cos(qJ(4));
t292 = pkin(4) * t262;
t261 = sin(qJ(1));
t290 = Icges(2,4) * t261;
t263 = cos(qJ(2));
t264 = cos(qJ(1));
t281 = t263 * t264;
t260 = sin(qJ(2));
t284 = t260 * t261;
t221 = -t257 * t281 + t284;
t258 = sin(qJ(4));
t289 = t221 * t258;
t282 = t261 * t263;
t283 = t260 * t264;
t223 = t257 * t282 + t283;
t288 = t223 * t258;
t256 = sin(pkin(5));
t287 = t256 * t261;
t286 = t256 * t263;
t285 = t256 * t264;
t280 = qJD(2) * t256;
t279 = V_base(5) * pkin(6) + V_base(1);
t276 = t258 * t286;
t275 = t256 * t294;
t233 = t261 * t280 + V_base(4);
t250 = V_base(6) + qJD(1);
t201 = qJD(3) * t223 + t233;
t234 = qJD(2) * t257 + t250;
t224 = -t257 * t284 + t281;
t259 = sin(qJ(3));
t204 = t224 * t259 - t261 * t275;
t166 = qJD(4) * t204 + t201;
t232 = -t264 * t280 + V_base(5);
t227 = pkin(1) * t261 - pkin(7) * t285;
t274 = -t227 * t250 + V_base(5) * t293 + t279;
t228 = pkin(1) * t264 + pkin(7) * t287;
t273 = V_base(4) * t227 - t228 * V_base(5) + V_base(3);
t200 = qJD(3) * t221 + t232;
t222 = t257 * t283 + t282;
t202 = t222 * t259 + t264 * t275;
t165 = qJD(4) * t202 + t200;
t217 = -qJD(3) * t286 + t234;
t272 = t250 * t228 + V_base(2) + (-pkin(6) - t293) * V_base(4);
t219 = t256 * t259 * t260 - t257 * t294;
t191 = qJD(4) * t219 + t217;
t192 = pkin(2) * t222 + pkin(8) * t221;
t226 = (pkin(2) * t260 - pkin(8) * t263) * t256;
t271 = -t192 * t234 + t232 * t226 + t274;
t193 = pkin(2) * t224 + pkin(8) * t223;
t270 = t233 * t192 - t193 * t232 + t273;
t269 = t234 * t193 - t226 * t233 + t272;
t203 = t222 * t294 - t259 * t285;
t163 = t203 * pkin(3) + t202 * pkin(9);
t220 = t257 * t259 + t260 * t275;
t190 = t220 * pkin(3) + t219 * pkin(9);
t268 = -t163 * t217 + t200 * t190 + t271;
t205 = t224 * t294 + t259 * t287;
t164 = t205 * pkin(3) + t204 * pkin(9);
t267 = t201 * t163 - t164 * t200 + t270;
t266 = t217 * t164 - t190 * t201 + t269;
t255 = qJ(4) + qJ(5);
t253 = Icges(2,4) * t264;
t252 = cos(t255);
t251 = sin(t255);
t242 = rSges(2,1) * t264 - rSges(2,2) * t261;
t241 = rSges(2,1) * t261 + rSges(2,2) * t264;
t240 = Icges(2,1) * t264 - t290;
t239 = Icges(2,1) * t261 + t253;
t238 = -Icges(2,2) * t261 + t253;
t237 = Icges(2,2) * t264 + t290;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = rSges(3,3) * t257 + (rSges(3,1) * t260 + rSges(3,2) * t263) * t256;
t212 = Icges(3,5) * t257 + (Icges(3,1) * t260 + Icges(3,4) * t263) * t256;
t211 = Icges(3,6) * t257 + (Icges(3,4) * t260 + Icges(3,2) * t263) * t256;
t210 = Icges(3,3) * t257 + (Icges(3,5) * t260 + Icges(3,6) * t263) * t256;
t209 = V_base(5) * rSges(2,3) - t241 * t250 + t279;
t208 = t242 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t206 = t241 * V_base(4) - t242 * V_base(5) + V_base(3);
t199 = t220 * t262 - t276;
t198 = -t220 * t258 - t262 * t286;
t195 = t220 * t252 - t251 * t286;
t194 = -t220 * t251 - t252 * t286;
t189 = rSges(3,1) * t224 - rSges(3,2) * t223 + rSges(3,3) * t287;
t188 = rSges(3,1) * t222 - rSges(3,2) * t221 - rSges(3,3) * t285;
t187 = Icges(3,1) * t224 - Icges(3,4) * t223 + Icges(3,5) * t287;
t186 = Icges(3,1) * t222 - Icges(3,4) * t221 - Icges(3,5) * t285;
t185 = Icges(3,4) * t224 - Icges(3,2) * t223 + Icges(3,6) * t287;
t184 = Icges(3,4) * t222 - Icges(3,2) * t221 - Icges(3,6) * t285;
t183 = Icges(3,5) * t224 - Icges(3,6) * t223 + Icges(3,3) * t287;
t182 = Icges(3,5) * t222 - Icges(3,6) * t221 - Icges(3,3) * t285;
t181 = rSges(4,1) * t220 - rSges(4,2) * t219 - rSges(4,3) * t286;
t180 = Icges(4,1) * t220 - Icges(4,4) * t219 - Icges(4,5) * t286;
t179 = Icges(4,4) * t220 - Icges(4,2) * t219 - Icges(4,6) * t286;
t178 = Icges(4,5) * t220 - Icges(4,6) * t219 - Icges(4,3) * t286;
t175 = t205 * t262 + t288;
t174 = -t205 * t258 + t223 * t262;
t173 = t203 * t262 + t289;
t172 = -t203 * t258 + t221 * t262;
t171 = qJD(5) * t219 + t191;
t170 = t205 * t252 + t223 * t251;
t169 = -t205 * t251 + t223 * t252;
t168 = t203 * t252 + t221 * t251;
t167 = -t203 * t251 + t221 * t252;
t161 = rSges(4,1) * t205 - rSges(4,2) * t204 + rSges(4,3) * t223;
t160 = rSges(4,1) * t203 - rSges(4,2) * t202 + rSges(4,3) * t221;
t158 = Icges(4,1) * t205 - Icges(4,4) * t204 + Icges(4,5) * t223;
t157 = Icges(4,1) * t203 - Icges(4,4) * t202 + Icges(4,5) * t221;
t156 = Icges(4,4) * t205 - Icges(4,2) * t204 + Icges(4,6) * t223;
t155 = Icges(4,4) * t203 - Icges(4,2) * t202 + Icges(4,6) * t221;
t154 = Icges(4,5) * t205 - Icges(4,6) * t204 + Icges(4,3) * t223;
t153 = Icges(4,5) * t203 - Icges(4,6) * t202 + Icges(4,3) * t221;
t152 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t219;
t151 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t219;
t150 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t219;
t149 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t219;
t148 = -pkin(4) * t276 + pkin(10) * t219 + t220 * t292;
t147 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t219;
t146 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t219;
t145 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t219;
t144 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t219;
t143 = qJD(5) * t204 + t166;
t142 = qJD(5) * t202 + t165;
t140 = -t188 * t234 + t213 * t232 + t274;
t139 = t189 * t234 - t213 * t233 + t272;
t138 = rSges(5,1) * t175 + rSges(5,2) * t174 + rSges(5,3) * t204;
t137 = rSges(5,1) * t173 + rSges(5,2) * t172 + rSges(5,3) * t202;
t136 = Icges(5,1) * t175 + Icges(5,4) * t174 + Icges(5,5) * t204;
t135 = Icges(5,1) * t173 + Icges(5,4) * t172 + Icges(5,5) * t202;
t134 = Icges(5,4) * t175 + Icges(5,2) * t174 + Icges(5,6) * t204;
t133 = Icges(5,4) * t173 + Icges(5,2) * t172 + Icges(5,6) * t202;
t132 = Icges(5,5) * t175 + Icges(5,6) * t174 + Icges(5,3) * t204;
t131 = Icges(5,5) * t173 + Icges(5,6) * t172 + Icges(5,3) * t202;
t130 = t188 * t233 - t189 * t232 + t273;
t129 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t204;
t128 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t202;
t127 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t204;
t126 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t202;
t125 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t204;
t124 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t202;
t123 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t204;
t122 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t202;
t121 = pkin(4) * t288 + pkin(10) * t204 + t205 * t292;
t120 = pkin(4) * t289 + pkin(10) * t202 + t203 * t292;
t119 = -t160 * t217 + t181 * t200 + t271;
t118 = t161 * t217 - t181 * t201 + t269;
t117 = t160 * t201 - t161 * t200 + t270;
t116 = -t137 * t191 + t152 * t165 + t268;
t115 = t138 * t191 - t152 * t166 + t266;
t114 = t137 * t166 - t138 * t165 + t267;
t113 = -t120 * t191 - t128 * t171 + t142 * t147 + t148 * t165 + t268;
t112 = t121 * t191 + t129 * t171 - t143 * t147 - t148 * t166 + t266;
t111 = t120 * t166 - t121 * t165 + t128 * t143 - t129 * t142 + t267;
t1 = m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(2) * (t206 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(3) * (t130 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + t233 * ((t183 * t287 - t223 * t185 + t224 * t187) * t233 + (t182 * t287 - t184 * t223 + t186 * t224) * t232 + (t210 * t287 - t211 * t223 + t212 * t224) * t234) / 0.2e1 + t232 * ((-t183 * t285 - t185 * t221 + t187 * t222) * t233 + (-t182 * t285 - t221 * t184 + t222 * t186) * t232 + (-t210 * t285 - t211 * t221 + t212 * t222) * t234) / 0.2e1 + t234 * ((t182 * t232 + t183 * t233 + t210 * t234) * t257 + ((t185 * t263 + t187 * t260) * t233 + (t184 * t263 + t186 * t260) * t232 + (t211 * t263 + t212 * t260) * t234) * t256) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t201 * ((t223 * t154 - t204 * t156 + t205 * t158) * t201 + (t153 * t223 - t155 * t204 + t157 * t205) * t200 + (t178 * t223 - t179 * t204 + t180 * t205) * t217) / 0.2e1 + t200 * ((t154 * t221 - t156 * t202 + t158 * t203) * t201 + (t221 * t153 - t202 * t155 + t203 * t157) * t200 + (t178 * t221 - t179 * t202 + t180 * t203) * t217) / 0.2e1 + t217 * ((-t154 * t286 - t156 * t219 + t158 * t220) * t201 + (-t153 * t286 - t155 * t219 + t157 * t220) * t200 + (-t178 * t286 - t219 * t179 + t220 * t180) * t217) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + t166 * ((t204 * t132 + t174 * t134 + t175 * t136) * t166 + (t131 * t204 + t133 * t174 + t135 * t175) * t165 + (t149 * t204 + t150 * t174 + t151 * t175) * t191) / 0.2e1 + t165 * ((t132 * t202 + t134 * t172 + t136 * t173) * t166 + (t202 * t131 + t172 * t133 + t173 * t135) * t165 + (t149 * t202 + t150 * t172 + t151 * t173) * t191) / 0.2e1 + t191 * ((t132 * t219 + t134 * t198 + t136 * t199) * t166 + (t131 * t219 + t133 * t198 + t135 * t199) * t165 + (t219 * t149 + t198 * t150 + t199 * t151) * t191) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + t143 * ((t204 * t123 + t169 * t125 + t170 * t127) * t143 + (t122 * t204 + t124 * t169 + t126 * t170) * t142 + (t144 * t204 + t145 * t169 + t146 * t170) * t171) / 0.2e1 + t142 * ((t123 * t202 + t125 * t167 + t127 * t168) * t143 + (t202 * t122 + t167 * t124 + t168 * t126) * t142 + (t144 * t202 + t145 * t167 + t146 * t168) * t171) / 0.2e1 + t171 * ((t123 * t219 + t125 * t194 + t127 * t195) * t143 + (t122 * t219 + t124 * t194 + t126 * t195) * t142 + (t219 * t144 + t194 * t145 + t195 * t146) * t171) / 0.2e1 + ((-t237 * t261 + t239 * t264 + Icges(1,4)) * V_base(5) + (-t261 * t238 + t264 * t240 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t264 * t237 + t261 * t239 + Icges(1,2)) * V_base(5) + (t238 * t264 + t240 * t261 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t261 + Icges(2,6) * t264) * V_base(5) + (Icges(2,5) * t264 - Icges(2,6) * t261) * V_base(4) + Icges(2,3) * t250 / 0.2e1) * t250;
T = t1;
