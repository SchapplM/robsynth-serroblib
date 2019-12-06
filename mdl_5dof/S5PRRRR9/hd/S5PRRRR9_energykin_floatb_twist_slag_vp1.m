% Calculate kinetic energy for
% S5PRRRR9
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:57
% EndTime: 2019-12-05 17:19:01
% DurationCPUTime: 3.21s
% Computational Cost: add. (2208->360), mult. (4776->537), div. (0->0), fcn. (5730->12), ass. (0->157)
t290 = cos(qJ(3));
t254 = cos(pkin(5));
t289 = pkin(6) * t254;
t258 = cos(qJ(4));
t288 = pkin(4) * t258;
t251 = sin(pkin(10));
t286 = Icges(2,4) * t251;
t253 = cos(pkin(10));
t257 = sin(qJ(2));
t259 = cos(qJ(2));
t278 = t254 * t259;
t213 = t251 * t257 - t253 * t278;
t255 = sin(qJ(4));
t285 = t213 * t255;
t215 = t251 * t278 + t253 * t257;
t284 = t215 * t255;
t252 = sin(pkin(5));
t283 = t251 * t252;
t282 = t252 * t253;
t256 = sin(qJ(3));
t281 = t252 * t256;
t280 = t252 * t259;
t279 = t254 * t257;
t277 = qJD(2) * t252;
t276 = V_base(5) * qJ(1) + V_base(1);
t272 = qJD(1) + V_base(3);
t271 = t255 * t280;
t270 = t252 * t290;
t229 = t251 * t277 + V_base(4);
t240 = qJD(2) * t254 + V_base(6);
t196 = qJD(3) * t215 + t229;
t216 = -t251 * t279 + t253 * t259;
t199 = t216 * t256 - t251 * t270;
t162 = qJD(4) * t199 + t196;
t228 = -t253 * t277 + V_base(5);
t195 = qJD(3) * t213 + t228;
t217 = -qJD(3) * t280 + t240;
t223 = pkin(1) * t251 - pkin(6) * t282;
t269 = -t223 * V_base(6) + V_base(5) * t289 + t276;
t224 = pkin(1) * t253 + pkin(6) * t283;
t268 = V_base(4) * t223 - t224 * V_base(5) + t272;
t214 = t251 * t259 + t253 * t279;
t197 = t214 * t256 + t253 * t270;
t161 = qJD(4) * t197 + t195;
t220 = -t254 * t290 + t257 * t281;
t189 = qJD(4) * t220 + t217;
t267 = V_base(6) * t224 + V_base(2) + (-qJ(1) - t289) * V_base(4);
t186 = pkin(2) * t214 + pkin(7) * t213;
t222 = (pkin(2) * t257 - pkin(7) * t259) * t252;
t266 = -t186 * t240 + t228 * t222 + t269;
t187 = pkin(2) * t216 + pkin(7) * t215;
t265 = t229 * t186 - t187 * t228 + t268;
t264 = t240 * t187 - t222 * t229 + t267;
t198 = t214 * t290 - t253 * t281;
t159 = t198 * pkin(3) + t197 * pkin(8);
t221 = t254 * t256 + t257 * t270;
t188 = t221 * pkin(3) + t220 * pkin(8);
t263 = -t159 * t217 + t195 * t188 + t266;
t200 = t216 * t290 + t251 * t281;
t160 = t200 * pkin(3) + t199 * pkin(8);
t262 = t196 * t159 - t160 * t195 + t265;
t261 = t217 * t160 - t188 * t196 + t264;
t250 = qJ(4) + qJ(5);
t248 = cos(t250);
t247 = sin(t250);
t246 = Icges(2,4) * t253;
t237 = rSges(2,1) * t253 - rSges(2,2) * t251;
t236 = rSges(2,1) * t251 + rSges(2,2) * t253;
t235 = Icges(2,1) * t253 - t286;
t234 = Icges(2,1) * t251 + t246;
t233 = -Icges(2,2) * t251 + t246;
t232 = Icges(2,2) * t253 + t286;
t227 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t226 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t225 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t209 = rSges(3,3) * t254 + (rSges(3,1) * t257 + rSges(3,2) * t259) * t252;
t208 = Icges(3,5) * t254 + (Icges(3,1) * t257 + Icges(3,4) * t259) * t252;
t207 = Icges(3,6) * t254 + (Icges(3,4) * t257 + Icges(3,2) * t259) * t252;
t206 = Icges(3,3) * t254 + (Icges(3,5) * t257 + Icges(3,6) * t259) * t252;
t205 = V_base(5) * rSges(2,3) - t236 * V_base(6) + t276;
t204 = t237 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t202 = t221 * t258 - t271;
t201 = -t221 * t255 - t258 * t280;
t194 = t236 * V_base(4) - t237 * V_base(5) + t272;
t191 = t221 * t248 - t247 * t280;
t190 = -t221 * t247 - t248 * t280;
t185 = rSges(4,1) * t221 - rSges(4,2) * t220 - rSges(4,3) * t280;
t184 = Icges(4,1) * t221 - Icges(4,4) * t220 - Icges(4,5) * t280;
t183 = Icges(4,4) * t221 - Icges(4,2) * t220 - Icges(4,6) * t280;
t182 = Icges(4,5) * t221 - Icges(4,6) * t220 - Icges(4,3) * t280;
t181 = rSges(3,1) * t216 - rSges(3,2) * t215 + rSges(3,3) * t283;
t180 = rSges(3,1) * t214 - rSges(3,2) * t213 - rSges(3,3) * t282;
t179 = Icges(3,1) * t216 - Icges(3,4) * t215 + Icges(3,5) * t283;
t178 = Icges(3,1) * t214 - Icges(3,4) * t213 - Icges(3,5) * t282;
t177 = Icges(3,4) * t216 - Icges(3,2) * t215 + Icges(3,6) * t283;
t176 = Icges(3,4) * t214 - Icges(3,2) * t213 - Icges(3,6) * t282;
t175 = Icges(3,5) * t216 - Icges(3,6) * t215 + Icges(3,3) * t283;
t174 = Icges(3,5) * t214 - Icges(3,6) * t213 - Icges(3,3) * t282;
t171 = qJD(5) * t220 + t189;
t170 = t200 * t258 + t284;
t169 = -t200 * t255 + t215 * t258;
t168 = t198 * t258 + t285;
t167 = -t198 * t255 + t213 * t258;
t166 = t200 * t248 + t215 * t247;
t165 = -t200 * t247 + t215 * t248;
t164 = t198 * t248 + t213 * t247;
t163 = -t198 * t247 + t213 * t248;
t157 = rSges(5,1) * t202 + rSges(5,2) * t201 + rSges(5,3) * t220;
t155 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t220;
t154 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t220;
t153 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t220;
t152 = rSges(4,1) * t200 - rSges(4,2) * t199 + rSges(4,3) * t215;
t151 = rSges(4,1) * t198 - rSges(4,2) * t197 + rSges(4,3) * t213;
t150 = Icges(4,1) * t200 - Icges(4,4) * t199 + Icges(4,5) * t215;
t149 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t213;
t148 = Icges(4,4) * t200 - Icges(4,2) * t199 + Icges(4,6) * t215;
t147 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t213;
t146 = Icges(4,5) * t200 - Icges(4,6) * t199 + Icges(4,3) * t215;
t145 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t213;
t144 = -pkin(4) * t271 + pkin(9) * t220 + t221 * t288;
t143 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t220;
t142 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t220;
t141 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t220;
t140 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t220;
t139 = qJD(5) * t199 + t162;
t138 = qJD(5) * t197 + t161;
t136 = -t180 * t240 + t209 * t228 + t269;
t135 = t181 * t240 - t209 * t229 + t267;
t134 = rSges(5,1) * t170 + rSges(5,2) * t169 + rSges(5,3) * t199;
t133 = rSges(5,1) * t168 + rSges(5,2) * t167 + rSges(5,3) * t197;
t132 = Icges(5,1) * t170 + Icges(5,4) * t169 + Icges(5,5) * t199;
t131 = Icges(5,1) * t168 + Icges(5,4) * t167 + Icges(5,5) * t197;
t130 = Icges(5,4) * t170 + Icges(5,2) * t169 + Icges(5,6) * t199;
t129 = Icges(5,4) * t168 + Icges(5,2) * t167 + Icges(5,6) * t197;
t128 = Icges(5,5) * t170 + Icges(5,6) * t169 + Icges(5,3) * t199;
t127 = Icges(5,5) * t168 + Icges(5,6) * t167 + Icges(5,3) * t197;
t126 = t180 * t229 - t181 * t228 + t268;
t125 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t199;
t124 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t197;
t123 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t199;
t122 = Icges(6,1) * t164 + Icges(6,4) * t163 + Icges(6,5) * t197;
t121 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t199;
t120 = Icges(6,4) * t164 + Icges(6,2) * t163 + Icges(6,6) * t197;
t119 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t199;
t118 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t197;
t117 = pkin(4) * t284 + pkin(9) * t199 + t200 * t288;
t116 = pkin(4) * t285 + pkin(9) * t197 + t198 * t288;
t115 = -t151 * t217 + t185 * t195 + t266;
t114 = t152 * t217 - t185 * t196 + t264;
t113 = t151 * t196 - t152 * t195 + t265;
t112 = -t133 * t189 + t157 * t161 + t263;
t111 = t134 * t189 - t157 * t162 + t261;
t110 = t133 * t162 - t134 * t161 + t262;
t109 = -t116 * t189 - t124 * t171 + t138 * t143 + t144 * t161 + t263;
t108 = t117 * t189 + t125 * t171 - t139 * t143 - t144 * t162 + t261;
t107 = t116 * t162 - t117 * t161 + t124 * t139 - t125 * t138 + t262;
t1 = m(1) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + m(2) * (t194 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(3) * (t126 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t229 * ((t175 * t283 - t177 * t215 + t179 * t216) * t229 + (t174 * t283 - t176 * t215 + t178 * t216) * t228 + (t206 * t283 - t207 * t215 + t208 * t216) * t240) / 0.2e1 + t228 * ((-t175 * t282 - t177 * t213 + t179 * t214) * t229 + (-t174 * t282 - t176 * t213 + t178 * t214) * t228 + (-t206 * t282 - t207 * t213 + t208 * t214) * t240) / 0.2e1 + t240 * ((t174 * t228 + t175 * t229 + t206 * t240) * t254 + ((t177 * t259 + t179 * t257) * t229 + (t176 * t259 + t178 * t257) * t228 + (t207 * t259 + t208 * t257) * t240) * t252) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t196 * ((t146 * t215 - t148 * t199 + t150 * t200) * t196 + (t145 * t215 - t147 * t199 + t149 * t200) * t195 + (t182 * t215 - t183 * t199 + t184 * t200) * t217) / 0.2e1 + t195 * ((t146 * t213 - t148 * t197 + t150 * t198) * t196 + (t145 * t213 - t147 * t197 + t149 * t198) * t195 + (t182 * t213 - t183 * t197 + t184 * t198) * t217) / 0.2e1 + t217 * ((-t146 * t280 - t148 * t220 + t150 * t221) * t196 + (-t145 * t280 - t147 * t220 + t149 * t221) * t195 + (-t182 * t280 - t183 * t220 + t221 * t184) * t217) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t162 * ((t199 * t128 + t169 * t130 + t170 * t132) * t162 + (t127 * t199 + t129 * t169 + t131 * t170) * t161 + (t153 * t199 + t154 * t169 + t155 * t170) * t189) / 0.2e1 + t161 * ((t128 * t197 + t130 * t167 + t132 * t168) * t162 + (t197 * t127 + t167 * t129 + t168 * t131) * t161 + (t153 * t197 + t154 * t167 + t155 * t168) * t189) / 0.2e1 + t189 * ((t128 * t220 + t130 * t201 + t132 * t202) * t162 + (t127 * t220 + t129 * t201 + t131 * t202) * t161 + (t153 * t220 + t154 * t201 + t155 * t202) * t189) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t139 * ((t199 * t119 + t165 * t121 + t166 * t123) * t139 + (t118 * t199 + t120 * t165 + t122 * t166) * t138 + (t140 * t199 + t141 * t165 + t142 * t166) * t171) / 0.2e1 + t138 * ((t119 * t197 + t121 * t163 + t123 * t164) * t139 + (t197 * t118 + t163 * t120 + t164 * t122) * t138 + (t140 * t197 + t141 * t163 + t142 * t164) * t171) / 0.2e1 + t171 * ((t119 * t220 + t121 * t190 + t123 * t191) * t139 + (t118 * t220 + t120 * t190 + t122 * t191) * t138 + (t140 * t220 + t141 * t190 + t142 * t191) * t171) / 0.2e1 + ((-t232 * t251 + t234 * t253 + Icges(1,4)) * V_base(5) + (-t233 * t251 + t235 * t253 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t232 * t253 + t234 * t251 + Icges(1,2)) * V_base(5) + (t233 * t253 + t235 * t251 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t251 + Icges(2,6) * t253 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t253 - Icges(2,6) * t251 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
