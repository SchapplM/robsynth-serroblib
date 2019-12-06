% Calculate kinetic energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:42
% EndTime: 2019-12-05 15:04:44
% DurationCPUTime: 2.82s
% Computational Cost: add. (1423->328), mult. (2196->463), div. (0->0), fcn. (2278->10), ass. (0->149)
t285 = -Icges(5,3) - Icges(4,3);
t229 = sin(pkin(7));
t231 = cos(pkin(7));
t284 = Icges(2,5) * t231 - Icges(2,6) * t229 + Icges(1,5);
t283 = Icges(2,5) * t229 + Icges(2,6) * t231 + Icges(1,6);
t227 = qJ(3) + pkin(9);
t223 = sin(t227);
t224 = cos(t227);
t230 = cos(pkin(8));
t266 = t229 * t230;
t172 = t223 * t266 + t224 * t231;
t173 = -t223 * t231 + t224 * t266;
t236 = cos(qJ(3));
t261 = t231 * t236;
t234 = sin(qJ(3));
t265 = t229 * t234;
t186 = -t230 * t265 - t261;
t262 = t231 * t234;
t264 = t229 * t236;
t187 = t230 * t264 - t262;
t228 = sin(pkin(8));
t270 = t228 * t229;
t282 = Icges(4,5) * t187 + Icges(5,5) * t173 + Icges(4,6) * t186 - Icges(5,6) * t172 - t285 * t270;
t263 = t230 * t231;
t174 = t223 * t263 - t229 * t224;
t175 = t223 * t229 + t224 * t263;
t188 = -t230 * t262 + t264;
t189 = t230 * t261 + t265;
t269 = t228 * t231;
t281 = Icges(4,5) * t189 + Icges(5,5) * t175 + Icges(4,6) * t188 - Icges(5,6) * t174 - t285 * t269;
t280 = t285 * t230 + (Icges(4,5) * t236 + Icges(5,5) * t224 - Icges(4,6) * t234 - Icges(5,6) * t223) * t228;
t276 = pkin(3) * t236;
t274 = Icges(2,4) * t229;
t273 = Icges(3,4) * t228;
t272 = Icges(3,4) * t230;
t271 = t223 * t228;
t233 = sin(qJ(5));
t268 = t228 * t233;
t235 = cos(qJ(5));
t267 = t228 * t235;
t260 = qJD(3) * t228;
t259 = qJD(4) * t228;
t258 = V_base(5) * qJ(1) + V_base(1);
t254 = qJD(1) + V_base(3);
t199 = t231 * t260 + V_base(4);
t198 = t229 * t260 + V_base(5);
t253 = qJD(2) * t229 + t258;
t210 = pkin(1) * t229 - qJ(2) * t231;
t252 = V_base(4) * t210 + t254;
t251 = pkin(2) * t230 + pkin(5) * t228;
t216 = -qJD(3) * t230 + V_base(6);
t250 = rSges(3,1) * t230 - rSges(3,2) * t228;
t249 = Icges(3,1) * t230 - t273;
t248 = -Icges(3,2) * t228 + t272;
t247 = Icges(3,5) * t230 - Icges(3,6) * t228;
t212 = pkin(1) * t231 + qJ(2) * t229;
t246 = -qJD(2) * t231 + V_base(6) * t212 + V_base(2);
t245 = qJ(4) * t228 + t230 * t276;
t190 = t251 * t229;
t214 = pkin(2) * t228 - pkin(5) * t230;
t244 = V_base(5) * t214 + (-t190 - t210) * V_base(6) + t253;
t191 = t251 * t231;
t243 = V_base(4) * t190 + (-t191 - t212) * V_base(5) + t252;
t242 = (-Icges(3,3) * t231 + t229 * t247) * V_base(5) + (Icges(3,3) * t229 + t231 * t247) * V_base(4) + (Icges(3,5) * t228 + Icges(3,6) * t230) * V_base(6);
t241 = V_base(6) * t191 + (-qJ(1) - t214) * V_base(4) + t246;
t157 = -qJ(4) * t230 + t228 * t276;
t240 = t198 * t157 + t231 * t259 + t244;
t146 = pkin(3) * t265 + t231 * t245;
t239 = t216 * t146 + t229 * t259 + t241;
t145 = -pkin(3) * t262 + t229 * t245;
t238 = -qJD(4) * t230 + t199 * t145 + t243;
t166 = -Icges(3,6) * t231 + t229 * t248;
t167 = Icges(3,6) * t229 + t231 * t248;
t168 = -Icges(3,5) * t231 + t229 * t249;
t169 = Icges(3,5) * t229 + t231 * t249;
t203 = Icges(3,2) * t230 + t273;
t206 = Icges(3,1) * t228 + t272;
t237 = (-t167 * t228 + t169 * t230) * V_base(4) + (-t166 * t228 + t168 * t230) * V_base(5) + (-t203 * t228 + t206 * t230) * V_base(6);
t225 = Icges(2,4) * t231;
t213 = rSges(2,1) * t231 - rSges(2,2) * t229;
t211 = rSges(2,1) * t229 + rSges(2,2) * t231;
t209 = rSges(3,1) * t228 + rSges(3,2) * t230;
t208 = Icges(2,1) * t231 - t274;
t207 = Icges(2,1) * t229 + t225;
t205 = -Icges(2,2) * t229 + t225;
t204 = Icges(2,2) * t231 + t274;
t197 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t196 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t195 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t185 = t224 * t267 - t230 * t233;
t184 = -t224 * t268 - t230 * t235;
t183 = qJD(5) * t271 + t216;
t180 = (pkin(4) * t224 + pkin(6) * t223) * t228;
t179 = -t230 * rSges(4,3) + (rSges(4,1) * t236 - rSges(4,2) * t234) * t228;
t178 = -Icges(4,5) * t230 + (Icges(4,1) * t236 - Icges(4,4) * t234) * t228;
t177 = -Icges(4,6) * t230 + (Icges(4,4) * t236 - Icges(4,2) * t234) * t228;
t171 = rSges(3,3) * t229 + t231 * t250;
t170 = -rSges(3,3) * t231 + t229 * t250;
t163 = -rSges(5,3) * t230 + (rSges(5,1) * t224 - rSges(5,2) * t223) * t228;
t162 = V_base(5) * rSges(2,3) - t211 * V_base(6) + t258;
t161 = t213 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t160 = -Icges(5,5) * t230 + (Icges(5,1) * t224 - Icges(5,4) * t223) * t228;
t159 = -Icges(5,6) * t230 + (Icges(5,4) * t224 - Icges(5,2) * t223) * t228;
t156 = t211 * V_base(4) - t213 * V_base(5) + t254;
t155 = t175 * t235 + t231 * t268;
t154 = -t175 * t233 + t231 * t267;
t153 = t173 * t235 + t229 * t268;
t152 = -t173 * t233 + t229 * t267;
t151 = qJD(5) * t174 + t199;
t150 = qJD(5) * t172 + t198;
t148 = rSges(4,1) * t189 + rSges(4,2) * t188 + rSges(4,3) * t269;
t147 = rSges(4,1) * t187 + rSges(4,2) * t186 + rSges(4,3) * t270;
t144 = Icges(4,1) * t189 + Icges(4,4) * t188 + Icges(4,5) * t269;
t143 = Icges(4,1) * t187 + Icges(4,4) * t186 + Icges(4,5) * t270;
t142 = Icges(4,4) * t189 + Icges(4,2) * t188 + Icges(4,6) * t269;
t141 = Icges(4,4) * t187 + Icges(4,2) * t186 + Icges(4,6) * t270;
t138 = pkin(4) * t175 + pkin(6) * t174;
t137 = pkin(4) * t173 + pkin(6) * t172;
t136 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t271;
t135 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t271;
t134 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t271;
t133 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t271;
t132 = rSges(5,1) * t175 - rSges(5,2) * t174 + rSges(5,3) * t269;
t131 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t270;
t129 = Icges(5,1) * t175 - Icges(5,4) * t174 + Icges(5,5) * t269;
t128 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t270;
t127 = Icges(5,4) * t175 - Icges(5,2) * t174 + Icges(5,6) * t269;
t126 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t270;
t122 = t209 * V_base(5) + (-t170 - t210) * V_base(6) + t253;
t121 = t171 * V_base(6) + (-qJ(1) - t209) * V_base(4) + t246;
t120 = t170 * V_base(4) + (-t171 - t212) * V_base(5) + t252;
t119 = rSges(6,1) * t155 + rSges(6,2) * t154 + rSges(6,3) * t174;
t118 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t172;
t117 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t174;
t116 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t172;
t115 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t174;
t114 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t172;
t113 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t174;
t112 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t172;
t111 = -t147 * t216 + t179 * t198 + t244;
t110 = t148 * t216 - t179 * t199 + t241;
t109 = t147 * t199 - t148 * t198 + t243;
t108 = t163 * t198 + (-t131 - t145) * t216 + t240;
t107 = t132 * t216 + (-t157 - t163) * t199 + t239;
t106 = t131 * t199 + (-t132 - t146) * t198 + t238;
t105 = -t118 * t183 + t136 * t150 + t180 * t198 + (-t137 - t145) * t216 + t240;
t104 = t119 * t183 - t136 * t151 + t138 * t216 + (-t157 - t180) * t199 + t239;
t103 = t118 * t151 - t119 * t150 + t137 * t199 + (-t138 - t146) * t198 + t238;
t1 = m(1) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + m(2) * (t156 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + t151 * ((t174 * t113 + t154 * t115 + t155 * t117) * t151 + (t112 * t174 + t114 * t154 + t116 * t155) * t150 + (t133 * t174 + t134 * t154 + t135 * t155) * t183) / 0.2e1 + t150 * ((t113 * t172 + t115 * t152 + t117 * t153) * t151 + (t172 * t112 + t152 * t114 + t153 * t116) * t150 + (t133 * t172 + t134 * t152 + t135 * t153) * t183) / 0.2e1 + t183 * ((t113 * t271 + t115 * t184 + t117 * t185) * t151 + (t112 * t271 + t114 * t184 + t116 * t185) * t150 + (t133 * t271 + t134 * t184 + t135 * t185) * t183) / 0.2e1 + ((-t159 * t172 + t160 * t173 + t177 * t186 + t178 * t187 + t270 * t280) * t216 + (-t127 * t172 + t129 * t173 + t142 * t186 + t144 * t187 + t270 * t281) * t199 + (-t126 * t172 + t128 * t173 + t141 * t186 + t143 * t187 + t270 * t282) * t198) * t198 / 0.2e1 + ((-t159 * t174 + t160 * t175 + t177 * t188 + t178 * t189 + t269 * t280) * t216 + (-t127 * t174 + t129 * t175 + t142 * t188 + t144 * t189 + t269 * t281) * t199 + (-t126 * t174 + t128 * t175 + t141 * t188 + t143 * t189 + t269 * t282) * t198) * t199 / 0.2e1 + ((-t198 * t282 - t199 * t281 - t216 * t280) * t230 + ((-t159 * t223 + t160 * t224 - t177 * t234 + t178 * t236) * t216 + (-t127 * t223 + t129 * t224 - t142 * t234 + t144 * t236) * t199 + (-t126 * t223 + t128 * t224 - t141 * t234 + t143 * t236) * t198) * t228) * t216 / 0.2e1 + (t242 * t229 + t237 * t231 + t284 * V_base(6) + (-t204 * t229 + t207 * t231 + Icges(1,4)) * V_base(5) + (-t205 * t229 + t208 * t231 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t237 * t229 - t242 * t231 + t283 * V_base(6) + (t204 * t231 + t207 * t229 + Icges(1,2)) * V_base(5) + (t205 * t231 + t208 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t203 * t230 + t206 * t228 + Icges(1,3) + Icges(2,3)) * V_base(6) + (t166 * t230 + t168 * t228 + t283) * V_base(5) + (t167 * t230 + t169 * t228 + t284) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
