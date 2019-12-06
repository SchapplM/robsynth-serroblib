% Calculate kinetic energy for
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:46
% EndTime: 2019-12-05 14:57:48
% DurationCPUTime: 2.23s
% Computational Cost: add. (1396->336), mult. (2151->473), div. (0->0), fcn. (2233->10), ass. (0->152)
t226 = sin(pkin(7));
t229 = cos(pkin(7));
t275 = Icges(2,5) * t229 - Icges(2,6) * t226 + Icges(1,5);
t274 = Icges(2,5) * t226 + Icges(2,6) * t229 + Icges(1,6);
t227 = cos(pkin(9));
t273 = pkin(3) * t227;
t272 = Icges(2,4) * t226;
t225 = sin(pkin(8));
t271 = Icges(3,4) * t225;
t228 = cos(pkin(8));
t270 = Icges(3,4) * t228;
t223 = pkin(9) + qJ(4);
t219 = sin(t223);
t269 = t219 * t225;
t224 = sin(pkin(9));
t268 = t224 * t226;
t267 = t224 * t229;
t266 = t225 * t226;
t265 = t225 * t229;
t231 = sin(qJ(5));
t264 = t225 * t231;
t232 = cos(qJ(5));
t263 = t225 * t232;
t262 = t226 * t228;
t261 = t228 * t229;
t205 = pkin(2) * t225 - qJ(3) * t228;
t260 = -qJ(1) - t205;
t246 = pkin(2) * t228 + qJ(3) * t225;
t186 = t246 * t226;
t207 = pkin(1) * t226 - qJ(2) * t229;
t258 = -t186 - t207;
t187 = t246 * t229;
t209 = pkin(1) * t229 + qJ(2) * t226;
t257 = -t187 - t209;
t256 = qJD(3) * t225;
t255 = qJD(4) * t225;
t254 = V_base(5) * qJ(1) + V_base(1);
t250 = qJD(1) + V_base(3);
t195 = t229 * t255 + V_base(4);
t194 = t226 * t255 + V_base(5);
t249 = qJD(2) * t226 + t254;
t248 = V_base(4) * t207 + t250;
t212 = -qJD(4) * t228 + V_base(6);
t247 = rSges(3,1) * t228 - rSges(3,2) * t225;
t245 = Icges(3,1) * t228 - t271;
t244 = -Icges(3,2) * t225 + t270;
t243 = Icges(3,5) * t228 - Icges(3,6) * t225;
t242 = -qJD(2) * t229 + V_base(6) * t209 + V_base(2);
t241 = V_base(5) * t205 + t229 * t256 + t249;
t240 = V_base(6) * t187 + t226 * t256 + t242;
t239 = -qJD(3) * t228 + V_base(4) * t186 + t248;
t238 = pkin(5) * t225 + t228 * t273;
t237 = (-Icges(3,3) * t229 + t226 * t243) * V_base(5) + (Icges(3,3) * t226 + t229 * t243) * V_base(4) + (Icges(3,5) * t225 + Icges(3,6) * t228) * V_base(6);
t139 = -pkin(3) * t267 + t226 * t238;
t153 = -pkin(5) * t228 + t225 * t273;
t236 = V_base(5) * t153 + (-t139 + t258) * V_base(6) + t241;
t140 = pkin(3) * t268 + t229 * t238;
t235 = V_base(6) * t140 + (-t153 + t260) * V_base(4) + t240;
t234 = V_base(4) * t139 + (-t140 + t257) * V_base(5) + t239;
t164 = -Icges(3,6) * t229 + t226 * t244;
t165 = Icges(3,6) * t226 + t229 * t244;
t167 = -Icges(3,5) * t229 + t226 * t245;
t168 = Icges(3,5) * t226 + t229 * t245;
t199 = Icges(3,2) * t228 + t271;
t202 = Icges(3,1) * t225 + t270;
t233 = (-t165 * t225 + t168 * t228) * V_base(4) + (-t164 * t225 + t167 * t228) * V_base(5) + (-t199 * t225 + t202 * t228) * V_base(6);
t221 = Icges(2,4) * t229;
t220 = cos(t223);
t210 = rSges(2,1) * t229 - rSges(2,2) * t226;
t208 = rSges(2,1) * t226 + rSges(2,2) * t229;
t206 = rSges(3,1) * t225 + rSges(3,2) * t228;
t204 = Icges(2,1) * t229 - t272;
t203 = Icges(2,1) * t226 + t221;
t201 = -Icges(2,2) * t226 + t221;
t200 = Icges(2,2) * t229 + t272;
t193 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t192 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t191 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t185 = t227 * t261 + t268;
t184 = -t224 * t261 + t226 * t227;
t183 = t227 * t262 - t267;
t182 = -t224 * t262 - t227 * t229;
t181 = t220 * t263 - t228 * t231;
t180 = -t220 * t264 - t228 * t232;
t179 = qJD(5) * t269 + t212;
t176 = (pkin(4) * t220 + pkin(6) * t219) * t225;
t175 = t219 * t226 + t220 * t261;
t174 = t219 * t261 - t226 * t220;
t173 = -t219 * t229 + t220 * t262;
t172 = t219 * t262 + t220 * t229;
t171 = rSges(3,3) * t226 + t229 * t247;
t170 = -rSges(3,3) * t229 + t226 * t247;
t169 = -rSges(4,3) * t228 + (rSges(4,1) * t227 - rSges(4,2) * t224) * t225;
t166 = -Icges(4,5) * t228 + (Icges(4,1) * t227 - Icges(4,4) * t224) * t225;
t163 = -Icges(4,6) * t228 + (Icges(4,4) * t227 - Icges(4,2) * t224) * t225;
t160 = -Icges(4,3) * t228 + (Icges(4,5) * t227 - Icges(4,6) * t224) * t225;
t159 = -rSges(5,3) * t228 + (rSges(5,1) * t220 - rSges(5,2) * t219) * t225;
t158 = V_base(5) * rSges(2,3) - t208 * V_base(6) + t254;
t157 = t210 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t156 = -Icges(5,5) * t228 + (Icges(5,1) * t220 - Icges(5,4) * t219) * t225;
t155 = -Icges(5,6) * t228 + (Icges(5,4) * t220 - Icges(5,2) * t219) * t225;
t154 = -Icges(5,3) * t228 + (Icges(5,5) * t220 - Icges(5,6) * t219) * t225;
t151 = t208 * V_base(4) - t210 * V_base(5) + t250;
t150 = t175 * t232 + t229 * t264;
t149 = -t175 * t231 + t229 * t263;
t148 = t173 * t232 + t226 * t264;
t147 = -t173 * t231 + t226 * t263;
t146 = qJD(5) * t174 + t195;
t145 = qJD(5) * t172 + t194;
t144 = pkin(4) * t175 + pkin(6) * t174;
t143 = pkin(4) * t173 + pkin(6) * t172;
t142 = rSges(4,1) * t185 + rSges(4,2) * t184 + rSges(4,3) * t265;
t141 = rSges(4,1) * t183 + rSges(4,2) * t182 + rSges(4,3) * t266;
t138 = Icges(4,1) * t185 + Icges(4,4) * t184 + Icges(4,5) * t265;
t137 = Icges(4,1) * t183 + Icges(4,4) * t182 + Icges(4,5) * t266;
t136 = Icges(4,4) * t185 + Icges(4,2) * t184 + Icges(4,6) * t265;
t135 = Icges(4,4) * t183 + Icges(4,2) * t182 + Icges(4,6) * t266;
t134 = Icges(4,5) * t185 + Icges(4,6) * t184 + Icges(4,3) * t265;
t133 = Icges(4,5) * t183 + Icges(4,6) * t182 + Icges(4,3) * t266;
t130 = rSges(6,1) * t181 + rSges(6,2) * t180 + rSges(6,3) * t269;
t129 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t269;
t128 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t269;
t127 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t269;
t126 = rSges(5,1) * t175 - rSges(5,2) * t174 + rSges(5,3) * t265;
t125 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t266;
t124 = Icges(5,1) * t175 - Icges(5,4) * t174 + Icges(5,5) * t265;
t123 = Icges(5,1) * t173 - Icges(5,4) * t172 + Icges(5,5) * t266;
t122 = Icges(5,4) * t175 - Icges(5,2) * t174 + Icges(5,6) * t265;
t121 = Icges(5,4) * t173 - Icges(5,2) * t172 + Icges(5,6) * t266;
t120 = Icges(5,5) * t175 - Icges(5,6) * t174 + Icges(5,3) * t265;
t119 = Icges(5,5) * t173 - Icges(5,6) * t172 + Icges(5,3) * t266;
t118 = t206 * V_base(5) + (-t170 - t207) * V_base(6) + t249;
t117 = t171 * V_base(6) + (-qJ(1) - t206) * V_base(4) + t242;
t116 = t170 * V_base(4) + (-t171 - t209) * V_base(5) + t248;
t115 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t174;
t114 = rSges(6,1) * t148 + rSges(6,2) * t147 + rSges(6,3) * t172;
t113 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t174;
t112 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t172;
t111 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t174;
t110 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t172;
t109 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t174;
t108 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t172;
t107 = t169 * V_base(5) + (-t141 + t258) * V_base(6) + t241;
t106 = t142 * V_base(6) + (-t169 + t260) * V_base(4) + t240;
t105 = t141 * V_base(4) + (-t142 + t257) * V_base(5) + t239;
t104 = -t125 * t212 + t159 * t194 + t236;
t103 = t126 * t212 - t159 * t195 + t235;
t102 = t125 * t195 - t126 * t194 + t234;
t101 = -t114 * t179 + t130 * t145 - t143 * t212 + t176 * t194 + t236;
t100 = t115 * t179 - t130 * t146 + t144 * t212 - t176 * t195 + t235;
t99 = t114 * t146 - t115 * t145 + t143 * t195 - t144 * t194 + t234;
t1 = m(1) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(2) * (t151 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t195 * ((t120 * t265 - t174 * t122 + t175 * t124) * t195 + (t119 * t265 - t121 * t174 + t123 * t175) * t194 + (t154 * t265 - t155 * t174 + t156 * t175) * t212) / 0.2e1 + t194 * ((t120 * t266 - t122 * t172 + t124 * t173) * t195 + (t119 * t266 - t172 * t121 + t173 * t123) * t194 + (t154 * t266 - t155 * t172 + t156 * t173) * t212) / 0.2e1 + t212 * ((-t119 * t194 - t120 * t195 - t154 * t212) * t228 + ((-t122 * t219 + t124 * t220) * t195 + (-t121 * t219 + t123 * t220) * t194 + (-t155 * t219 + t156 * t220) * t212) * t225) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + t146 * ((t174 * t109 + t149 * t111 + t150 * t113) * t146 + (t108 * t174 + t110 * t149 + t112 * t150) * t145 + (t127 * t174 + t128 * t149 + t129 * t150) * t179) / 0.2e1 + t145 * ((t109 * t172 + t111 * t147 + t113 * t148) * t146 + (t172 * t108 + t147 * t110 + t148 * t112) * t145 + (t127 * t172 + t128 * t147 + t129 * t148) * t179) / 0.2e1 + t179 * ((t109 * t269 + t111 * t180 + t113 * t181) * t146 + (t108 * t269 + t110 * t180 + t112 * t181) * t145 + (t127 * t269 + t180 * t128 + t181 * t129) * t179) / 0.2e1 + (t237 * t226 + t233 * t229 + (t160 * t265 + t163 * t184 + t166 * t185 + t275) * V_base(6) + (t133 * t265 + t135 * t184 + t137 * t185 - t200 * t226 + t203 * t229 + Icges(1,4)) * V_base(5) + (t134 * t265 + t184 * t136 + t185 * t138 - t226 * t201 + t229 * t204 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t233 * t226 - t237 * t229 + (t160 * t266 + t163 * t182 + t166 * t183 + t274) * V_base(6) + (t133 * t266 + t182 * t135 + t183 * t137 + t229 * t200 + t226 * t203 + Icges(1,2)) * V_base(5) + (t134 * t266 + t136 * t182 + t138 * t183 + t201 * t229 + t204 * t226 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(1,3) + Icges(2,3) + (t199 - t160) * t228 + (-t224 * t163 + t227 * t166 + t202) * t225) * V_base(6) + ((t164 - t133) * t228 + (-t135 * t224 + t137 * t227 + t167) * t225 + t274) * V_base(5) + ((t165 - t134) * t228 + (-t136 * t224 + t138 * t227 + t168) * t225 + t275) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
