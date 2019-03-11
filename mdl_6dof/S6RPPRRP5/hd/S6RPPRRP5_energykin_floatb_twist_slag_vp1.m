% Calculate kinetic energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:37
% EndTime: 2019-03-09 02:07:40
% DurationCPUTime: 2.41s
% Computational Cost: add. (865->257), mult. (1553->342), div. (0->0), fcn. (1391->6), ass. (0->130)
t289 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t288 = Icges(6,1) + Icges(7,1);
t287 = Icges(6,4) + Icges(7,4);
t286 = Icges(7,5) + Icges(6,5);
t285 = Icges(6,2) + Icges(7,2);
t284 = Icges(7,6) + Icges(6,6);
t283 = Icges(7,3) + Icges(6,3);
t282 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t281 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t280 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t279 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t204 = sin(qJ(1));
t278 = t289 * t204;
t207 = cos(qJ(1));
t277 = t289 * t207;
t203 = sin(qJ(4));
t205 = cos(qJ(5));
t240 = t205 * t207;
t202 = sin(qJ(5));
t243 = t204 * t202;
t141 = -t203 * t243 + t240;
t242 = t204 * t205;
t244 = t202 * t207;
t142 = t203 * t242 + t244;
t206 = cos(qJ(4));
t241 = t204 * t206;
t276 = t284 * t141 + t286 * t142 - t283 * t241;
t143 = -t203 * t244 - t242;
t144 = t203 * t240 - t243;
t239 = t206 * t207;
t275 = t284 * t143 + t286 * t144 - t283 * t239;
t274 = t285 * t141 + t287 * t142 - t284 * t241;
t273 = t285 * t143 + t287 * t144 - t284 * t239;
t272 = t287 * t141 + t288 * t142 - t286 * t241;
t271 = t287 * t143 + t288 * t144 - t286 * t239;
t270 = (-t284 * t202 + t286 * t205) * t206 + t283 * t203;
t269 = (-t285 * t202 + t287 * t205) * t206 + t284 * t203;
t268 = (-t287 * t202 + t288 * t205) * t206 + t286 * t203;
t267 = t282 * t207 - t278;
t266 = t282 * t204 + t277;
t265 = -t279 * t207 - t278;
t264 = t279 * t204 - t277;
t261 = -pkin(2) - pkin(6);
t256 = pkin(7) * t204;
t255 = pkin(7) * t207;
t254 = pkin(5) * t205;
t251 = Icges(5,4) * t203;
t250 = Icges(5,4) * t206;
t246 = qJ(3) * t204;
t245 = qJ(3) * t207;
t215 = -qJ(6) * t206 + t203 * t254;
t238 = rSges(7,1) * t142 + rSges(7,2) * t141 - rSges(7,3) * t241 + pkin(5) * t244 + t204 * t215;
t237 = t144 * rSges(7,1) + t143 * rSges(7,2) - rSges(7,3) * t239 - pkin(5) * t243 + t207 * t215;
t236 = (rSges(7,1) * t205 - rSges(7,2) * t202 + t254) * t206 + (qJ(6) + rSges(7,3)) * t203;
t235 = qJD(5) * t206;
t234 = qJD(6) * t206;
t174 = t204 * pkin(1) - qJ(2) * t207;
t233 = V_base(4) * t174 + V_base(3);
t232 = V_base(5) * pkin(6) + V_base(1);
t185 = qJD(4) * t207 + V_base(5);
t191 = V_base(6) + qJD(1);
t229 = V_base(4) * t246 + t233;
t228 = qJD(2) * t204 + t232;
t227 = -t174 - t246;
t179 = pkin(1) * t207 + t204 * qJ(2);
t226 = -t179 - t245;
t225 = pkin(4) * t203 - pkin(8) * t206;
t186 = -qJD(4) * t204 + V_base(4);
t224 = rSges(5,1) * t203 + rSges(5,2) * t206;
t223 = Icges(5,1) * t203 + t250;
t222 = Icges(5,2) * t206 + t251;
t221 = Icges(5,5) * t203 + Icges(5,6) * t206;
t220 = -qJD(2) * t207 + t191 * t179 + V_base(2);
t219 = V_base(5) * pkin(2) + qJD(3) * t207 + t228;
t218 = t227 - t255;
t217 = V_base(5) * pkin(3) + t219;
t216 = qJD(3) * t204 + t191 * t245 + t220;
t214 = (Icges(5,3) * t207 + t204 * t221) * t185 + (-Icges(5,3) * t204 + t207 * t221) * t186 + (Icges(5,5) * t206 - Icges(5,6) * t203) * t191;
t213 = V_base(4) * t255 + t229 + (t226 + t256) * V_base(5);
t212 = (-pkin(3) + t261) * V_base(4) + t216;
t146 = t225 * t204;
t147 = t225 * t207;
t211 = t186 * t146 - t185 * t147 + t213;
t184 = pkin(4) * t206 + pkin(8) * t203;
t210 = t191 * t147 - t186 * t184 + t212;
t209 = t185 * t184 + (-t146 + t218) * t191 + t217;
t128 = Icges(5,6) * t207 + t204 * t222;
t129 = -Icges(5,6) * t204 + t207 * t222;
t132 = Icges(5,5) * t207 + t204 * t223;
t133 = -Icges(5,5) * t204 + t207 * t223;
t163 = -Icges(5,2) * t203 + t250;
t170 = Icges(5,1) * t206 - t251;
t208 = (t129 * t206 + t133 * t203) * t186 + (t128 * t206 + t132 * t203) * t185 + (t163 * t206 + t170 * t203) * t191;
t182 = rSges(2,1) * t207 - t204 * rSges(2,2);
t181 = -rSges(3,2) * t207 + t204 * rSges(3,3);
t180 = -rSges(4,2) * t207 + t204 * rSges(4,3);
t178 = rSges(5,1) * t206 - rSges(5,2) * t203;
t177 = t204 * rSges(2,1) + rSges(2,2) * t207;
t176 = -t204 * rSges(3,2) - rSges(3,3) * t207;
t175 = t204 * rSges(4,2) + rSges(4,3) * t207;
t173 = qJD(5) * t203 + t191;
t151 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t150 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t149 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t140 = -t207 * t235 + t186;
t139 = -t204 * t235 + t185;
t137 = -t204 * rSges(5,3) + t207 * t224;
t136 = rSges(6,3) * t203 + (rSges(6,1) * t205 - rSges(6,2) * t202) * t206;
t134 = rSges(5,3) * t207 + t204 * t224;
t118 = V_base(5) * rSges(2,3) - t177 * t191 + t232;
t117 = t182 * t191 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t116 = t177 * V_base(4) - t182 * V_base(5) + V_base(3);
t115 = t144 * rSges(6,1) + t143 * rSges(6,2) - rSges(6,3) * t239;
t113 = rSges(6,1) * t142 + rSges(6,2) * t141 - rSges(6,3) * t241;
t97 = V_base(5) * rSges(3,1) + (-t174 - t176) * t191 + t228;
t96 = t191 * t181 + (-rSges(3,1) - pkin(6)) * V_base(4) + t220;
t95 = t176 * V_base(4) + (-t179 - t181) * V_base(5) + t233;
t94 = V_base(5) * rSges(4,1) + (-t180 + t227) * t191 + t219;
t93 = t191 * t175 + (-rSges(4,1) + t261) * V_base(4) + t216;
t92 = V_base(4) * t180 + (-t175 + t226) * V_base(5) + t229;
t91 = t185 * t178 + (-t134 + t218) * t191 + t217;
t90 = -t186 * t178 + (t137 - t256) * t191 + t212;
t89 = t186 * t134 - t185 * t137 + t213;
t88 = -t173 * t113 + t139 * t136 + t209;
t87 = t173 * t115 - t140 * t136 - t191 * t256 + t210;
t86 = t140 * t113 - t139 * t115 + t211;
t85 = t139 * t236 - t173 * t238 - t207 * t234 + t209;
t84 = (-pkin(7) * t191 - t234) * t204 + t237 * t173 - t236 * t140 + t210;
t83 = qJD(6) * t203 - t139 * t237 + t140 * t238 + t211;
t1 = m(1) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(3) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t185 * (t204 * t208 + t207 * t214) / 0.2e1 + t186 * (-t204 * t214 + t208 * t207) / 0.2e1 + ((t141 * t269 + t142 * t268 - t270 * t241) * t173 + (t141 * t273 + t271 * t142 - t275 * t241) * t140 + (t274 * t141 + t272 * t142 - t276 * t241) * t139) * t139 / 0.2e1 + ((t143 * t269 + t144 * t268 - t239 * t270) * t173 + (t273 * t143 + t271 * t144 - t275 * t239) * t140 + (t274 * t143 + t272 * t144 - t239 * t276) * t139) * t140 / 0.2e1 + (((-t202 * t269 + t205 * t268) * t173 + (-t202 * t273 + t205 * t271) * t140 + (-t202 * t274 + t205 * t272) * t139) * t206 + (t139 * t276 + t275 * t140 + t270 * t173) * t203) * t173 / 0.2e1 + ((-t129 * t203 + t133 * t206) * t186 + (-t128 * t203 + t132 * t206) * t185 + (-t163 * t203 + t170 * t206 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t191) * t191 / 0.2e1 + ((t204 * t265 + t207 * t266 + Icges(1,4)) * V_base(5) + (t264 * t204 + t267 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t266 * t204 - t265 * t207 + Icges(1,2)) * V_base(5) + (t204 * t267 - t207 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t191 * (t281 * t204 - t280 * t207) + V_base(4) * t191 * (t280 * t204 + t281 * t207) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
