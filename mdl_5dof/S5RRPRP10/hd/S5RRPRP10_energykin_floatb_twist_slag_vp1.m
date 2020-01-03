% Calculate kinetic energy for
% S5RRPRP10
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:20
% EndTime: 2019-12-31 20:09:23
% DurationCPUTime: 2.90s
% Computational Cost: add. (847->237), mult. (1600->329), div. (0->0), fcn. (1494->6), ass. (0->122)
t286 = Icges(3,4) + Icges(4,6);
t285 = Icges(3,1) + Icges(4,2);
t284 = -Icges(3,2) - Icges(4,3);
t198 = cos(qJ(2));
t283 = t286 * t198;
t195 = sin(qJ(2));
t282 = t286 * t195;
t281 = Icges(4,4) - Icges(3,5);
t280 = Icges(4,5) - Icges(3,6);
t279 = t284 * t195 + t283;
t278 = t285 * t198 - t282;
t277 = Icges(4,1) + Icges(3,3);
t276 = Icges(5,1) + Icges(6,1);
t275 = Icges(5,4) + Icges(6,4);
t274 = Icges(6,5) + Icges(5,5);
t273 = Icges(5,2) + Icges(6,2);
t272 = Icges(6,6) + Icges(5,6);
t271 = Icges(6,3) + Icges(5,3);
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t270 = t279 * t196 + t280 * t199;
t269 = -t280 * t196 + t279 * t199;
t268 = t278 * t196 + t281 * t199;
t267 = -t281 * t196 + t278 * t199;
t266 = t284 * t198 - t282;
t265 = t285 * t195 + t283;
t264 = t280 * t195 - t281 * t198;
t197 = cos(qJ(4));
t233 = t197 * t199;
t194 = sin(qJ(4));
t236 = t196 * t194;
t150 = t195 * t233 - t236;
t235 = t196 * t197;
t237 = t194 * t199;
t151 = t195 * t237 + t235;
t232 = t198 * t199;
t263 = t272 * t150 + t274 * t151 + t271 * t232;
t152 = t195 * t235 + t237;
t153 = t195 * t236 - t233;
t234 = t196 * t198;
t262 = t272 * t152 + t274 * t153 + t271 * t234;
t261 = t273 * t150 + t275 * t151 + t272 * t232;
t260 = t273 * t152 + t275 * t153 + t272 * t234;
t259 = t275 * t150 + t276 * t151 + t274 * t232;
t258 = t275 * t152 + t276 * t153 + t274 * t234;
t257 = (-t274 * t194 - t272 * t197) * t198 + t271 * t195;
t256 = (-t275 * t194 - t273 * t197) * t198 + t272 * t195;
t255 = (-t276 * t194 - t275 * t197) * t198 + t274 * t195;
t184 = -qJD(2) * t199 + V_base(5);
t185 = qJD(2) * t196 + V_base(4);
t189 = V_base(6) + qJD(1);
t254 = (t266 * t195 + t265 * t198) * t189 + (-t269 * t195 + t267 * t198) * t185 + (-t270 * t195 + t268 * t198) * t184;
t253 = (-t281 * t195 - t280 * t198) * t189 + (t277 * t196 + t264 * t199) * t185 + (t264 * t196 - t277 * t199) * t184;
t246 = pkin(4) * t194;
t245 = pkin(7) * t195;
t244 = pkin(4) * t197;
t242 = Icges(2,4) * t196;
t209 = qJ(5) * t198 + t195 * t246;
t231 = t151 * rSges(6,1) + t150 * rSges(6,2) + rSges(6,3) * t232 + t196 * t244 + t199 * t209;
t230 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t234 + t196 * t209 - t199 * t244;
t229 = (-rSges(6,1) * t194 - rSges(6,2) * t197 - t246) * t198 + (rSges(6,3) + qJ(5)) * t195;
t218 = pkin(2) * t198 + qJ(3) * t195;
t154 = t218 * t196;
t182 = t196 * pkin(1) - pkin(6) * t199;
t228 = -t154 - t182;
t227 = qJD(3) * t195;
t226 = qJD(4) * t198;
t225 = qJD(5) * t198;
t224 = V_base(5) * pkin(5) + V_base(1);
t177 = pkin(2) * t195 - qJ(3) * t198;
t221 = t184 * t177 + t199 * t227 + t224;
t220 = rSges(3,1) * t198 - rSges(3,2) * t195;
t219 = -rSges(4,2) * t198 + rSges(4,3) * t195;
t183 = pkin(1) * t199 + t196 * pkin(6);
t211 = -V_base(4) * pkin(5) + t189 * t183 + V_base(2);
t210 = V_base(4) * t182 - t183 * V_base(5) + V_base(3);
t155 = t218 * t199;
t206 = t189 * t155 + t196 * t227 + t211;
t205 = -qJD(3) * t198 + t185 * t154 + t210;
t160 = -pkin(3) * t199 + pkin(7) * t234;
t204 = t184 * t245 + (-t160 + t228) * t189 + t221;
t159 = t196 * pkin(3) + pkin(7) * t232;
t203 = t189 * t159 + (-t177 - t245) * t185 + t206;
t202 = t185 * t160 + (-t155 - t159) * t184 + t205;
t191 = Icges(2,4) * t199;
t181 = rSges(2,1) * t199 - t196 * rSges(2,2);
t180 = t196 * rSges(2,1) + rSges(2,2) * t199;
t179 = rSges(3,1) * t195 + rSges(3,2) * t198;
t178 = -rSges(4,2) * t195 - rSges(4,3) * t198;
t176 = qJD(4) * t195 + t189;
t175 = Icges(2,1) * t199 - t242;
t174 = Icges(2,1) * t196 + t191;
t172 = -Icges(2,2) * t196 + t191;
t171 = Icges(2,2) * t199 + t242;
t163 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t162 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t161 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t149 = t199 * t226 + t185;
t148 = t196 * t226 + t184;
t144 = -rSges(4,1) * t199 + t196 * t219;
t143 = t196 * rSges(4,1) + t199 * t219;
t142 = t196 * rSges(3,3) + t199 * t220;
t141 = rSges(5,3) * t195 + (-rSges(5,1) * t194 - rSges(5,2) * t197) * t198;
t139 = -rSges(3,3) * t199 + t196 * t220;
t117 = V_base(5) * rSges(2,3) - t180 * t189 + t224;
t116 = t181 * t189 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t115 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t112 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t234;
t110 = t151 * rSges(5,1) + t150 * rSges(5,2) + rSges(5,3) * t232;
t96 = t179 * t184 + (-t139 - t182) * t189 + t224;
t95 = t142 * t189 - t179 * t185 + t211;
t94 = t139 * t185 - t142 * t184 + t210;
t93 = t178 * t184 + (-t144 + t228) * t189 + t221;
t92 = t143 * t189 + (-t177 - t178) * t185 + t206;
t91 = t144 * t185 + (-t143 - t155) * t184 + t205;
t90 = -t112 * t176 + t141 * t148 + t204;
t89 = t110 * t176 - t141 * t149 + t203;
t88 = -t110 * t148 + t112 * t149 + t202;
t87 = t148 * t229 - t176 * t230 + t199 * t225 + t204;
t86 = -t149 * t229 + t176 * t231 + t196 * t225 + t203;
t85 = qJD(5) * t195 - t148 * t231 + t149 * t230 + t202;
t1 = m(1) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + ((t152 * t256 + t153 * t255 + t234 * t257) * t176 + (t261 * t152 + t259 * t153 + t234 * t263) * t149 + (t260 * t152 + t258 * t153 + t262 * t234) * t148) * t148 / 0.2e1 + ((t150 * t256 + t151 * t255 + t232 * t257) * t176 + (t261 * t150 + t259 * t151 + t263 * t232) * t149 + (t150 * t260 + t151 * t258 + t232 * t262) * t148) * t149 / 0.2e1 + (((-t194 * t255 - t197 * t256) * t176 + (-t194 * t259 - t197 * t261) * t149 + (-t194 * t258 - t197 * t260) * t148) * t198 + (t262 * t148 + t149 * t263 + t257 * t176) * t195) * t176 / 0.2e1 + (t254 * t196 - t253 * t199) * t184 / 0.2e1 + (t253 * t196 + t254 * t199) * t185 / 0.2e1 + ((-t196 * t171 + t174 * t199 + Icges(1,4)) * V_base(5) + (-t196 * t172 + t175 * t199 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t171 * t199 + t196 * t174 + Icges(1,2)) * V_base(5) + (t172 * t199 + t196 * t175 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t267 * t195 + t269 * t198) * t185 + (t268 * t195 + t270 * t198) * t184 + (t265 * t195 - t266 * t198 + Icges(2,3)) * t189) * t189 / 0.2e1 + t189 * V_base(4) * (Icges(2,5) * t199 - Icges(2,6) * t196) + V_base(5) * t189 * (Icges(2,5) * t196 + Icges(2,6) * t199) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
