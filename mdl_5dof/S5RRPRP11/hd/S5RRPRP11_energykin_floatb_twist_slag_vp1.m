% Calculate kinetic energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:18
% EndTime: 2019-12-31 20:12:20
% DurationCPUTime: 2.83s
% Computational Cost: add. (833->234), mult. (1597->321), div. (0->0), fcn. (1501->6), ass. (0->121)
t282 = Icges(3,4) + Icges(4,6);
t281 = Icges(3,1) + Icges(4,2);
t280 = -Icges(3,2) - Icges(4,3);
t196 = cos(qJ(2));
t279 = t282 * t196;
t193 = sin(qJ(2));
t278 = t282 * t193;
t277 = Icges(4,4) - Icges(3,5);
t276 = Icges(4,5) - Icges(3,6);
t275 = t280 * t193 + t279;
t274 = t281 * t196 - t278;
t273 = Icges(4,1) + Icges(3,3);
t272 = Icges(5,1) + Icges(6,1);
t271 = Icges(5,4) - Icges(6,5);
t270 = Icges(6,4) + Icges(5,5);
t269 = Icges(5,2) + Icges(6,3);
t268 = Icges(6,6) - Icges(5,6);
t267 = Icges(5,3) + Icges(6,2);
t194 = sin(qJ(1));
t197 = cos(qJ(1));
t266 = t275 * t194 + t276 * t197;
t265 = -t276 * t194 + t275 * t197;
t264 = t274 * t194 + t277 * t197;
t263 = -t277 * t194 + t274 * t197;
t262 = t280 * t196 - t278;
t261 = t281 * t193 + t279;
t260 = t276 * t193 - t277 * t196;
t259 = rSges(6,1) + pkin(4);
t258 = rSges(6,3) + qJ(5);
t195 = cos(qJ(4));
t230 = t195 * t197;
t192 = sin(qJ(4));
t234 = t192 * t194;
t148 = -t193 * t230 + t234;
t232 = t194 * t195;
t233 = t192 * t197;
t149 = t193 * t233 + t232;
t229 = t196 * t197;
t257 = t268 * t148 + t270 * t149 + t267 * t229;
t150 = t193 * t232 + t233;
t151 = t193 * t234 - t230;
t231 = t194 * t196;
t256 = -t268 * t150 + t270 * t151 + t267 * t231;
t255 = -t269 * t148 + t271 * t149 - t268 * t229;
t254 = t269 * t150 + t271 * t151 - t268 * t231;
t253 = -t271 * t148 + t272 * t149 + t270 * t229;
t252 = t271 * t150 + t272 * t151 + t270 * t231;
t251 = (t271 * t192 + t269 * t195) * t196 + t268 * t193;
t250 = (-t270 * t192 + t268 * t195) * t196 + t267 * t193;
t249 = (-t272 * t192 - t271 * t195) * t196 + t270 * t193;
t184 = -qJD(2) * t197 + V_base(5);
t185 = qJD(2) * t194 + V_base(4);
t188 = V_base(6) + qJD(1);
t248 = (t262 * t193 + t261 * t196) * t188 + (-t265 * t193 + t263 * t196) * t185 + (-t266 * t193 + t264 * t196) * t184;
t247 = (-t277 * t193 - t276 * t196) * t188 + (t273 * t194 + t260 * t197) * t185 + (t260 * t194 - t273 * t197) * t184;
t240 = pkin(7) * t193;
t239 = Icges(2,4) * t194;
t228 = rSges(6,2) * t229 + t258 * t148 + t259 * t149;
t227 = rSges(6,2) * t231 - t258 * t150 + t259 * t151;
t226 = rSges(6,2) * t193 + (-t259 * t192 + t258 * t195) * t196;
t215 = pkin(2) * t196 + qJ(3) * t193;
t152 = t215 * t194;
t181 = t194 * pkin(1) - pkin(6) * t197;
t225 = -t152 - t181;
t224 = qJD(3) * t193;
t223 = qJD(3) * t196;
t222 = qJD(4) * t196;
t221 = V_base(5) * pkin(5) + V_base(1);
t176 = pkin(2) * t193 - qJ(3) * t196;
t218 = t184 * t176 + t197 * t224 + t221;
t217 = rSges(3,1) * t196 - rSges(3,2) * t193;
t216 = -rSges(4,2) * t196 + rSges(4,3) * t193;
t182 = pkin(1) * t197 + t194 * pkin(6);
t208 = -V_base(4) * pkin(5) + t188 * t182 + V_base(2);
t207 = V_base(4) * t181 - t182 * V_base(5) + V_base(3);
t206 = t185 * t152 + t207;
t154 = t215 * t197;
t203 = t188 * t154 + t194 * t224 + t208;
t159 = -pkin(3) * t197 + pkin(7) * t231;
t202 = t184 * t240 + (-t159 + t225) * t188 + t218;
t158 = t194 * pkin(3) + pkin(7) * t229;
t201 = t185 * t159 + (-t154 - t158) * t184 + t206;
t200 = t188 * t158 + (-t176 - t240) * t185 + t203;
t190 = Icges(2,4) * t197;
t180 = rSges(2,1) * t197 - t194 * rSges(2,2);
t179 = t194 * rSges(2,1) + rSges(2,2) * t197;
t178 = rSges(3,1) * t193 + rSges(3,2) * t196;
t177 = -rSges(4,2) * t193 - rSges(4,3) * t196;
t175 = qJD(4) * t193 + t188;
t174 = Icges(2,1) * t197 - t239;
t173 = Icges(2,1) * t194 + t190;
t171 = -Icges(2,2) * t194 + t190;
t170 = Icges(2,2) * t197 + t239;
t162 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t161 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t160 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t147 = t197 * t222 + t185;
t146 = t194 * t222 + t184;
t143 = -rSges(4,1) * t197 + t194 * t216;
t142 = t194 * rSges(4,1) + t197 * t216;
t141 = t194 * rSges(3,3) + t197 * t217;
t140 = rSges(5,3) * t193 + (-rSges(5,1) * t192 - rSges(5,2) * t195) * t196;
t138 = -rSges(3,3) * t197 + t194 * t217;
t116 = V_base(5) * rSges(2,3) - t179 * t188 + t221;
t115 = t180 * t188 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t114 = t179 * V_base(4) - t180 * V_base(5) + V_base(3);
t111 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t231;
t109 = t149 * rSges(5,1) - t148 * rSges(5,2) + rSges(5,3) * t229;
t95 = t178 * t184 + (-t138 - t181) * t188 + t221;
t94 = t141 * t188 - t178 * t185 + t208;
t93 = t138 * t185 - t141 * t184 + t207;
t92 = t177 * t184 + (-t143 + t225) * t188 + t218;
t91 = t142 * t188 + (-t176 - t177) * t185 + t203;
t90 = -t223 + t143 * t185 + (-t142 - t154) * t184 + t206;
t89 = -t111 * t175 + t140 * t146 + t202;
t88 = t109 * t175 - t140 * t147 + t200;
t87 = -t109 * t146 + t111 * t147 + t201 - t223;
t86 = qJD(5) * t148 + t146 * t226 - t175 * t227 + t202;
t85 = -qJD(5) * t150 - t147 * t226 + t175 * t228 + t200;
t84 = (qJD(5) * t195 - qJD(3)) * t196 + t227 * t147 - t228 * t146 + t201;
t1 = m(1) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(2) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(3) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + ((-t251 * t150 + t249 * t151 + t250 * t231) * t175 + (t255 * t150 + t253 * t151 + t257 * t231) * t147 + (t254 * t150 + t252 * t151 + t256 * t231) * t146) * t146 / 0.2e1 + ((t251 * t148 + t249 * t149 + t250 * t229) * t175 + (-t255 * t148 + t253 * t149 + t257 * t229) * t147 + (-t254 * t148 + t252 * t149 + t256 * t229) * t146) * t147 / 0.2e1 + (((-t249 * t192 + t251 * t195) * t175 + (-t253 * t192 - t255 * t195) * t147 + (-t252 * t192 - t254 * t195) * t146) * t196 + (t256 * t146 + t257 * t147 + t250 * t175) * t193) * t175 / 0.2e1 + (t248 * t194 - t247 * t197) * t184 / 0.2e1 + (t247 * t194 + t248 * t197) * t185 / 0.2e1 + ((-t194 * t170 + t173 * t197 + Icges(1,4)) * V_base(5) + (-t194 * t171 + t174 * t197 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t170 * t197 + t194 * t173 + Icges(1,2)) * V_base(5) + (t171 * t197 + t194 * t174 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t263 * t193 + t265 * t196) * t185 + (t264 * t193 + t266 * t196) * t184 + (t261 * t193 - t262 * t196 + Icges(2,3)) * t188) * t188 / 0.2e1 + V_base(4) * t188 * (Icges(2,5) * t197 - Icges(2,6) * t194) + t188 * V_base(5) * (Icges(2,5) * t194 + Icges(2,6) * t197) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
