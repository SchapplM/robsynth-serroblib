% Calculate kinetic energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:55
% EndTime: 2019-12-05 16:11:57
% DurationCPUTime: 2.41s
% Computational Cost: add. (1153->277), mult. (2526->382), div. (0->0), fcn. (2745->8), ass. (0->131)
t287 = Icges(5,1) + Icges(6,1);
t286 = -Icges(5,4) + Icges(6,5);
t285 = Icges(6,4) + Icges(5,5);
t284 = Icges(5,2) + Icges(6,3);
t283 = Icges(6,2) + Icges(5,3);
t282 = -Icges(5,6) + Icges(6,6);
t281 = rSges(6,1) + pkin(4);
t280 = rSges(6,3) + qJ(5);
t223 = sin(pkin(7));
t224 = cos(pkin(7));
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t228 = cos(qJ(2));
t254 = t227 * t228;
t184 = t223 * t254 - t224 * t225;
t222 = sin(pkin(8));
t226 = sin(qJ(2));
t262 = cos(pkin(8));
t244 = t226 * t262;
t158 = t184 * t222 - t223 * t244;
t258 = t223 * t226;
t159 = t184 * t262 + t222 * t258;
t255 = t225 * t228;
t183 = t223 * t255 + t224 * t227;
t277 = t284 * t158 + t286 * t159 + t282 * t183;
t188 = t223 * t225 + t224 * t254;
t160 = t188 * t222 - t224 * t244;
t257 = t224 * t226;
t161 = t188 * t262 + t222 * t257;
t187 = -t223 * t227 + t224 * t255;
t276 = t284 * t160 + t286 * t161 + t282 * t187;
t275 = t282 * t158 + t285 * t159 + t283 * t183;
t274 = t282 * t160 + t285 * t161 + t283 * t187;
t273 = t286 * t158 + t287 * t159 + t285 * t183;
t272 = t286 * t160 + t287 * t161 + t285 * t187;
t185 = t226 * t227 * t222 + t228 * t262;
t186 = -t228 * t222 + t227 * t244;
t256 = t225 * t226;
t271 = t284 * t185 + t286 * t186 + t282 * t256;
t270 = t282 * t185 + t285 * t186 + t283 * t256;
t269 = t286 * t185 + t287 * t186 + t285 * t256;
t139 = Icges(4,4) * t184 - Icges(4,2) * t183 + Icges(4,6) * t258;
t268 = -t139 + t275;
t140 = Icges(4,4) * t188 - Icges(4,2) * t187 + Icges(4,6) * t257;
t267 = -t140 + t274;
t174 = -Icges(4,6) * t228 + (Icges(4,4) * t227 - Icges(4,2) * t225) * t226;
t266 = -t174 + t270;
t261 = Icges(2,4) * t223;
t260 = Icges(3,4) * t226;
t259 = Icges(3,4) * t228;
t253 = rSges(6,2) * t183 + t280 * t158 + t159 * t281;
t252 = rSges(6,2) * t187 + t280 * t160 + t161 * t281;
t251 = rSges(6,2) * t256 + t280 * t185 + t186 * t281;
t250 = qJD(3) * t226;
t249 = V_base(5) * qJ(1) + V_base(1);
t245 = qJD(1) + V_base(3);
t215 = qJD(2) * t223 + V_base(4);
t243 = pkin(2) * t228 + pkin(6) * t226;
t214 = -qJD(2) * t224 + V_base(5);
t242 = rSges(3,1) * t228 - rSges(3,2) * t226;
t241 = Icges(3,1) * t228 - t260;
t240 = -Icges(3,2) * t226 + t259;
t239 = Icges(3,5) * t228 - Icges(3,6) * t226;
t206 = pkin(1) * t224 + pkin(5) * t223;
t238 = -V_base(4) * qJ(1) + V_base(6) * t206 + V_base(2);
t205 = pkin(1) * t223 - pkin(5) * t224;
t237 = V_base(4) * t205 - t206 * V_base(5) + t245;
t189 = t243 * t223;
t213 = t226 * pkin(2) - pkin(6) * t228;
t236 = t214 * t213 + (-t189 - t205) * V_base(6) + t249;
t235 = (-Icges(3,3) * t224 + t223 * t239) * t214 + (Icges(3,3) * t223 + t224 * t239) * t215 + (Icges(3,5) * t226 + Icges(3,6) * t228) * V_base(6);
t190 = t243 * t224;
t234 = V_base(6) * t190 - t213 * t215 + t238;
t181 = t223 * t250 + t214;
t191 = (pkin(3) * t227 + qJ(4) * t225) * t226;
t233 = qJD(4) * t187 + t181 * t191 + t236;
t232 = t215 * t189 - t190 * t214 + t237;
t155 = pkin(3) * t188 + qJ(4) * t187;
t216 = -qJD(3) * t228 + V_base(6);
t231 = qJD(4) * t183 + t216 * t155 + t234;
t153 = pkin(3) * t184 + qJ(4) * t183;
t182 = t224 * t250 + t215;
t230 = qJD(4) * t256 + t182 * t153 + t232;
t167 = -Icges(3,6) * t224 + t223 * t240;
t168 = Icges(3,6) * t223 + t224 * t240;
t169 = -Icges(3,5) * t224 + t223 * t241;
t170 = Icges(3,5) * t223 + t224 * t241;
t208 = Icges(3,2) * t228 + t260;
t209 = Icges(3,1) * t226 + t259;
t229 = (-t168 * t226 + t170 * t228) * t215 + (-t167 * t226 + t169 * t228) * t214 + (-t208 * t226 + t209 * t228) * V_base(6);
t220 = Icges(2,4) * t224;
t210 = t226 * rSges(3,1) + rSges(3,2) * t228;
t204 = rSges(2,1) * t224 - rSges(2,2) * t223;
t203 = rSges(2,1) * t223 + rSges(2,2) * t224;
t202 = Icges(2,1) * t224 - t261;
t201 = Icges(2,1) * t223 + t220;
t200 = -Icges(2,2) * t223 + t220;
t199 = Icges(2,2) * t224 + t261;
t196 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t195 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t194 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t176 = -rSges(4,3) * t228 + (rSges(4,1) * t227 - rSges(4,2) * t225) * t226;
t175 = -Icges(4,5) * t228 + (Icges(4,1) * t227 - Icges(4,4) * t225) * t226;
t173 = -Icges(4,3) * t228 + (Icges(4,5) * t227 - Icges(4,6) * t225) * t226;
t172 = t223 * rSges(3,3) + t224 * t242;
t171 = -t224 * rSges(3,3) + t223 * t242;
t163 = V_base(5) * rSges(2,3) - t203 * V_base(6) + t249;
t162 = t204 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t157 = t203 * V_base(4) - t204 * V_base(5) + t245;
t152 = rSges(5,1) * t186 - rSges(5,2) * t185 + rSges(5,3) * t256;
t150 = rSges(4,1) * t188 - rSges(4,2) * t187 + rSges(4,3) * t257;
t149 = rSges(4,1) * t184 - rSges(4,2) * t183 + rSges(4,3) * t258;
t142 = Icges(4,1) * t188 - Icges(4,4) * t187 + Icges(4,5) * t257;
t141 = Icges(4,1) * t184 - Icges(4,4) * t183 + Icges(4,5) * t258;
t138 = Icges(4,5) * t188 - Icges(4,6) * t187 + Icges(4,3) * t257;
t137 = Icges(4,5) * t184 - Icges(4,6) * t183 + Icges(4,3) * t258;
t134 = t210 * t214 + (-t171 - t205) * V_base(6) + t249;
t133 = t172 * V_base(6) - t210 * t215 + t238;
t130 = rSges(5,1) * t161 - rSges(5,2) * t160 + rSges(5,3) * t187;
t128 = rSges(5,1) * t159 - rSges(5,2) * t158 + rSges(5,3) * t183;
t114 = t171 * t215 - t172 * t214 + t237;
t113 = -t149 * t216 + t176 * t181 + t236;
t112 = t150 * t216 - t176 * t182 + t234;
t111 = t149 * t182 - t150 * t181 + t232;
t110 = t152 * t181 + (-t128 - t153) * t216 + t233;
t109 = t130 * t216 + (-t152 - t191) * t182 + t231;
t108 = t128 * t182 + (-t130 - t155) * t181 + t230;
t107 = qJD(5) * t160 + t251 * t181 + (-t153 - t253) * t216 + t233;
t106 = qJD(5) * t158 + t252 * t216 + (-t191 - t251) * t182 + t231;
t105 = qJD(5) * t185 + t253 * t182 + (-t155 - t252) * t181 + t230;
t1 = m(1) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(2) * (t157 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t114 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + t215 * (t235 * t223 + t229 * t224) / 0.2e1 + t214 * (t229 * t223 - t235 * t224) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + ((-t199 * t223 + t201 * t224 + Icges(1,4)) * V_base(5) + (-t200 * t223 + t202 * t224 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t199 * t224 + t201 * t223 + Icges(1,2)) * V_base(5) + (t200 * t224 + t202 * t223 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t158 * t271 + t159 * t269 + t173 * t258 + t175 * t184 + t183 * t266) * t216 + (t138 * t258 + t142 * t184 + t158 * t276 + t159 * t272 + t183 * t267) * t182 + (t137 * t258 + t141 * t184 + t277 * t158 + t273 * t159 + t268 * t183) * t181) * t181 / 0.2e1 + ((t160 * t271 + t161 * t269 + t173 * t257 + t175 * t188 + t187 * t266) * t216 + (t138 * t257 + t142 * t188 + t276 * t160 + t272 * t161 + t267 * t187) * t182 + (t137 * t257 + t141 * t188 + t160 * t277 + t161 * t273 + t187 * t268) * t181) * t182 / 0.2e1 + ((-t173 * t228 + (-t174 * t225 + t175 * t227) * t226 + t270 * t256 + t269 * t186 + t271 * t185) * t216 + (-t138 * t228 + (-t140 * t225 + t142 * t227) * t226 + t274 * t256 + t272 * t186 + t276 * t185) * t182 + (-t137 * t228 + (-t139 * t225 + t141 * t227) * t226 + t275 * t256 + t273 * t186 + t277 * t185) * t181) * t216 / 0.2e1 + ((t168 * t228 + t226 * t170) * t215 + (t167 * t228 + t226 * t169) * t214 + (t208 * t228 + t226 * t209 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t224 - Icges(2,6) * t223 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t223 + Icges(2,6) * t224 + Icges(1,6));
T = t1;
