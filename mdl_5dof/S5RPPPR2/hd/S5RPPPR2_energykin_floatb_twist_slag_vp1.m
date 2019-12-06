% Calculate kinetic energy for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:54
% EndTime: 2019-12-05 17:30:56
% DurationCPUTime: 2.69s
% Computational Cost: add. (1293->336), mult. (2923->444), div. (0->0), fcn. (3349->10), ass. (0->151)
t232 = sin(pkin(8));
t235 = cos(pkin(7));
t239 = cos(qJ(1));
t234 = cos(pkin(8));
t237 = sin(qJ(1));
t268 = t237 * t234;
t196 = t232 * t239 - t235 * t268;
t231 = sin(pkin(9));
t233 = sin(pkin(7));
t278 = cos(pkin(9));
t257 = t233 * t278;
t168 = t196 * t231 + t237 * t257;
t272 = t233 * t237;
t169 = t196 * t278 - t231 * t272;
t269 = t237 * t232;
t195 = t234 * t239 + t235 * t269;
t130 = Icges(5,5) * t169 - Icges(5,6) * t168 - Icges(5,3) * t195;
t153 = Icges(4,4) * t196 + Icges(4,2) * t195 - Icges(4,6) * t272;
t285 = t130 - t153;
t270 = t235 * t239;
t198 = t234 * t270 + t269;
t170 = t198 * t231 - t239 * t257;
t271 = t233 * t239;
t171 = t198 * t278 + t231 * t271;
t197 = t232 * t270 - t268;
t131 = Icges(5,5) * t171 - Icges(5,6) * t170 + Icges(5,3) * t197;
t154 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t271;
t284 = t131 - t154;
t193 = t233 * t234 * t231 + t235 * t278;
t194 = -t235 * t231 + t234 * t257;
t273 = t232 * t233;
t147 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t273;
t177 = -Icges(4,6) * t235 + (Icges(4,4) * t234 - Icges(4,2) * t232) * t233;
t283 = t147 - t177;
t274 = Icges(3,4) * t235;
t251 = -Icges(3,2) * t233 + t274;
t183 = Icges(3,6) * t239 - t237 * t251;
t275 = Icges(3,4) * t233;
t252 = Icges(3,1) * t235 - t275;
t185 = Icges(3,5) * t239 - t237 * t252;
t276 = Icges(2,4) * t239;
t282 = -Icges(2,1) * t237 - t183 * t233 + t185 * t235 - t276;
t184 = Icges(3,6) * t237 + t239 * t251;
t186 = Icges(3,5) * t237 + t239 * t252;
t277 = Icges(2,4) * t237;
t281 = Icges(2,1) * t239 - t184 * t233 + t186 * t235 - t277;
t211 = pkin(2) * t233 - qJ(3) * t235;
t279 = -pkin(5) - t211;
t253 = pkin(2) * t235 + qJ(3) * t233;
t200 = t253 * t237;
t219 = -t237 * pkin(1) + qJ(2) * t239;
t267 = t200 - t219;
t201 = t253 * t239;
t221 = pkin(1) * t239 + t237 * qJ(2);
t266 = -t201 - t221;
t265 = qJD(3) * t233;
t264 = V_base(5) * t221 + V_base(1);
t263 = V_base(6) * pkin(5) + V_base(2);
t199 = (pkin(3) * t234 + qJ(4) * t232) * t233;
t260 = -t199 + t279;
t162 = pkin(3) * t196 - qJ(4) * t195;
t259 = -t162 + t267;
t163 = pkin(3) * t198 + qJ(4) * t197;
t258 = -t163 + t266;
t227 = V_base(4) + qJD(1);
t256 = qJD(2) * t237 + t227 * t219 + V_base(3);
t255 = qJD(2) * t239 + t263;
t254 = rSges(3,1) * t235 - rSges(3,2) * t233;
t250 = Icges(3,5) * t235 - Icges(3,6) * t233;
t209 = Icges(3,2) * t235 + t275;
t210 = Icges(3,1) * t233 + t274;
t247 = t209 * t233 - t210 * t235;
t246 = -t227 * t200 + t239 * t265 + t256;
t245 = -qJD(3) * t235 + V_base(5) * t201 + t264;
t244 = qJD(4) * t197 + t227 * t162 + t246;
t243 = V_base(6) * t211 - t237 * t265 + t255;
t242 = qJD(4) * t273 + V_base(5) * t163 + t245;
t241 = (Icges(3,3) * t239 - t237 * t250) * V_base(5) + (Icges(3,3) * t237 + t239 * t250) * V_base(6) + (Icges(3,5) * t233 + Icges(3,6) * t235) * t227;
t240 = -qJD(4) * t195 + V_base(6) * t199 + t243;
t238 = cos(qJ(5));
t236 = sin(qJ(5));
t222 = rSges(2,1) * t239 - t237 * rSges(2,2);
t220 = -t237 * rSges(2,1) - rSges(2,2) * t239;
t216 = -Icges(2,2) * t237 + t276;
t215 = -Icges(2,2) * t239 - t277;
t214 = Icges(2,5) * t239 - Icges(2,6) * t237;
t213 = -Icges(2,5) * t237 - Icges(2,6) * t239;
t212 = rSges(3,1) * t233 + rSges(3,2) * t235;
t207 = -V_base(5) * rSges(1,1) + rSges(1,2) * V_base(4) + V_base(3);
t206 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t205 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t188 = t237 * rSges(3,3) + t239 * t254;
t187 = rSges(3,3) * t239 - t237 * t254;
t179 = -rSges(4,3) * t235 + (rSges(4,1) * t234 - rSges(4,2) * t232) * t233;
t178 = -Icges(4,5) * t235 + (Icges(4,1) * t234 - Icges(4,4) * t232) * t233;
t176 = -Icges(4,3) * t235 + (Icges(4,5) * t234 - Icges(4,6) * t232) * t233;
t175 = qJD(5) * t193 + t227;
t174 = V_base(6) * rSges(2,3) - t222 * t227 + t263;
t173 = t220 * t227 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t172 = -t220 * V_base(6) + t222 * V_base(5) + V_base(1);
t167 = t194 * t238 + t236 * t273;
t166 = -t194 * t236 + t238 * t273;
t165 = qJD(5) * t168 + V_base(5);
t164 = qJD(5) * t170 + V_base(6);
t161 = pkin(4) * t194 + pkin(6) * t193;
t158 = t198 * rSges(4,1) - t197 * rSges(4,2) + rSges(4,3) * t271;
t157 = rSges(4,1) * t196 + rSges(4,2) * t195 - rSges(4,3) * t272;
t156 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t271;
t155 = Icges(4,1) * t196 + Icges(4,4) * t195 - Icges(4,5) * t272;
t152 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t271;
t151 = Icges(4,5) * t196 + Icges(4,6) * t195 - Icges(4,3) * t272;
t150 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t273;
t149 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t273;
t148 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t273;
t146 = t171 * t238 + t197 * t236;
t145 = -t171 * t236 + t197 * t238;
t144 = t169 * t238 - t195 * t236;
t143 = -t169 * t236 - t195 * t238;
t142 = pkin(4) * t171 + pkin(6) * t170;
t141 = pkin(4) * t169 + pkin(6) * t168;
t140 = t212 * V_base(6) + (-t188 - t221) * t227 + t255;
t139 = t187 * t227 + (-pkin(5) - t212) * V_base(5) + t256;
t138 = t188 * V_base(5) + (-t187 - t219) * V_base(6) + t264;
t137 = rSges(5,1) * t171 - rSges(5,2) * t170 + rSges(5,3) * t197;
t136 = rSges(5,1) * t169 - rSges(5,2) * t168 - rSges(5,3) * t195;
t135 = Icges(5,1) * t171 - Icges(5,4) * t170 + Icges(5,5) * t197;
t134 = Icges(5,1) * t169 - Icges(5,4) * t168 - Icges(5,5) * t195;
t133 = Icges(5,4) * t171 - Icges(5,2) * t170 + Icges(5,6) * t197;
t132 = Icges(5,4) * t169 - Icges(5,2) * t168 - Icges(5,6) * t195;
t129 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t193;
t128 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t193;
t127 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t193;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t193;
t125 = t179 * V_base(6) + (-t158 + t266) * t227 + t243;
t124 = t157 * t227 + (-t179 + t279) * V_base(5) + t246;
t123 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t170;
t122 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t168;
t121 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t170;
t120 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t168;
t119 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t170;
t118 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t168;
t117 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t170;
t116 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t168;
t115 = t158 * V_base(5) + (-t157 + t267) * V_base(6) + t245;
t114 = t150 * V_base(6) + (-t137 + t258) * t227 + t240;
t113 = t136 * t227 + (-t150 + t260) * V_base(5) + t244;
t112 = t137 * V_base(5) + (-t136 + t259) * V_base(6) + t242;
t111 = -t123 * t175 + t129 * t164 + t161 * V_base(6) + (-t142 + t258) * t227 + t240;
t110 = t122 * t175 - t129 * t165 + t141 * t227 + (-t161 + t260) * V_base(5) + t244;
t109 = -t122 * t164 + t123 * t165 + t142 * V_base(5) + (-t141 + t259) * V_base(6) + t242;
t1 = m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(3) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t175 * ((t126 * t193 + t127 * t166 + t128 * t167) * t175 + (t116 * t193 + t118 * t166 + t120 * t167) * t165 + (t117 * t193 + t119 * t166 + t121 * t167) * t164) / 0.2e1 + t165 * ((t126 * t168 + t127 * t143 + t128 * t144) * t175 + (t116 * t168 + t118 * t143 + t120 * t144) * t165 + (t117 * t168 + t119 * t143 + t121 * t144) * t164) / 0.2e1 + t164 * ((t126 * t170 + t127 * t145 + t128 * t146) * t175 + (t116 * t170 + t118 * t145 + t120 * t146) * t165 + (t117 * t170 + t119 * t145 + t121 * t146) * t164) / 0.2e1 + ((t131 * t273 - t133 * t193 + t135 * t194 + t214 + (t184 - t152) * t235 + (-t154 * t232 + t156 * t234 + t186) * t233) * V_base(6) + (t130 * t273 - t132 * t193 + t134 * t194 + t213 + (t183 - t151) * t235 + (-t153 * t232 + t155 * t234 + t185) * t233) * V_base(5) + (t147 * t273 - t148 * t193 + t149 * t194 + Icges(2,3) + (t209 - t176) * t235 + (-t177 * t232 + t178 * t234 + t210) * t233) * t227) * t227 / 0.2e1 + (t241 * t239 + (-t148 * t168 + t149 * t169 - t176 * t272 + t178 * t196 - t195 * t283 + t247 * t237 + t213) * t227 + (-t133 * t168 + t135 * t169 - t152 * t272 + t156 * t196 - t195 * t284 - t216 * t239 - t237 * t281 + Icges(1,6)) * V_base(6) + (-t132 * t168 + t134 * t169 - t151 * t272 + t155 * t196 - t285 * t195 - t215 * t239 - t282 * t237 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t241 * t237 + (-t148 * t170 + t149 * t171 + t176 * t271 + t198 * t178 + t197 * t283 - t247 * t239 + t214) * t227 + (-t133 * t170 + t135 * t171 + t152 * t271 + t198 * t156 + t284 * t197 - t237 * t216 + t281 * t239 + Icges(1,3)) * V_base(6) + (-t132 * t170 + t134 * t171 + t151 * t271 + t198 * t155 + t197 * t285 - t237 * t215 + t282 * t239 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
