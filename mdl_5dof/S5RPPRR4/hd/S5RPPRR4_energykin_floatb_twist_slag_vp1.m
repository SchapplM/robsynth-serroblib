% Calculate kinetic energy for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:14
% DurationCPUTime: 2.80s
% Computational Cost: add. (1471->359), mult. (1861->481), div. (0->0), fcn. (1811->10), ass. (0->163)
t225 = sin(pkin(9));
t271 = pkin(3) * t225;
t226 = sin(pkin(8));
t228 = cos(pkin(8));
t191 = pkin(2) * t226 - qJ(3) * t228;
t270 = -pkin(5) - t191;
t229 = qJ(3) + pkin(6);
t227 = cos(pkin(9));
t211 = t227 * pkin(3) + pkin(2);
t230 = sin(qJ(1));
t269 = Icges(2,4) * t230;
t268 = Icges(3,4) * t226;
t267 = Icges(3,4) * t228;
t266 = t226 * t230;
t231 = cos(qJ(1));
t265 = t226 * t231;
t264 = t228 * t231;
t224 = pkin(9) + qJ(4);
t214 = qJ(5) + t224;
t209 = sin(t214);
t263 = t230 * t209;
t210 = cos(t214);
t262 = t230 * t210;
t212 = sin(t224);
t261 = t230 * t212;
t213 = cos(t224);
t260 = t230 * t213;
t259 = t230 * t225;
t258 = t230 * t227;
t185 = pkin(2) * t228 + qJ(3) * t226 + pkin(1);
t220 = t230 * pkin(1);
t167 = t185 * t230 - t220;
t200 = -qJ(2) * t231 + t220;
t257 = -t167 - t200;
t221 = t231 * pkin(1);
t168 = t185 * t231 - t221;
t217 = t230 * qJ(2);
t202 = t221 + t217;
t256 = -t168 - t202;
t171 = t211 * t228 + t226 * t229 + pkin(1);
t255 = t171 - t185;
t208 = qJ(2) + t271;
t254 = pkin(4) * t212 - t208 + t271;
t253 = qJD(3) * t226;
t252 = qJD(4) * t226;
t251 = qJD(5) * t226;
t250 = V_base(4) * t200 + V_base(3);
t249 = V_base(5) * pkin(5) + V_base(1);
t187 = t231 * t252 + V_base(4);
t186 = t230 * t252 + V_base(5);
t215 = V_base(6) + qJD(1);
t246 = qJD(2) * t230 + t249;
t245 = rSges(3,1) * t228 - rSges(3,2) * t226;
t244 = Icges(3,1) * t228 - t268;
t243 = -Icges(3,2) * t226 + t267;
t242 = Icges(3,5) * t228 - Icges(3,6) * t226;
t241 = -qJD(2) * t231 + t215 * t202 + V_base(2);
t240 = V_base(5) * t191 + t231 * t253 + t246;
t180 = pkin(4) * t213 + t211;
t223 = -pkin(7) - t229;
t239 = t180 * t228 - t223 * t226 - t171;
t238 = -qJD(3) * t228 + V_base(4) * t167 + t250;
t237 = t215 * t168 + t230 * t253 + t241;
t236 = (-Icges(3,3) * t231 + t230 * t242) * V_base(5) + (Icges(3,3) * t230 + t231 * t242) * V_base(4) + (Icges(3,5) * t226 + Icges(3,6) * t228) * t215;
t121 = (qJ(2) - t208) * t231 + t255 * t230;
t138 = (qJ(3) - t229) * t228 + (-pkin(2) + t211) * t226;
t235 = V_base(5) * t138 + (-t121 + t257) * t215 + t240;
t122 = t208 * t230 + t231 * t255 - t217;
t234 = V_base(4) * t121 + (-t122 + t256) * V_base(5) + t238;
t233 = t215 * t122 + (-t138 + t270) * V_base(4) + t237;
t157 = -Icges(3,6) * t231 + t230 * t243;
t158 = Icges(3,6) * t230 + t231 * t243;
t159 = -Icges(3,5) * t231 + t230 * t244;
t160 = Icges(3,5) * t230 + t231 * t244;
t189 = Icges(3,2) * t228 + t268;
t190 = Icges(3,1) * t226 + t267;
t232 = (-t158 * t226 + t160 * t228) * V_base(4) + (-t157 * t226 + t159 * t228) * V_base(5) + (-t189 * t226 + t190 * t228) * t215;
t219 = Icges(2,4) * t231;
t203 = rSges(2,1) * t231 - t230 * rSges(2,2);
t201 = t230 * rSges(2,1) + rSges(2,2) * t231;
t199 = Icges(2,1) * t231 - t269;
t198 = Icges(2,1) * t230 + t219;
t197 = -Icges(2,2) * t230 + t219;
t196 = Icges(2,2) * t231 + t269;
t195 = Icges(2,5) * t231 - Icges(2,6) * t230;
t194 = Icges(2,5) * t230 + Icges(2,6) * t231;
t193 = -qJD(4) * t228 + t215;
t192 = rSges(3,1) * t226 + rSges(3,2) * t228;
t184 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t183 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t182 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t176 = (-qJD(4) - qJD(5)) * t228 + t215;
t175 = t227 * t264 + t259;
t174 = -t225 * t264 + t258;
t173 = -t225 * t231 + t228 * t258;
t172 = -t227 * t231 - t228 * t259;
t170 = t231 * t251 + t187;
t169 = t230 * t251 + t186;
t166 = t230 * rSges(3,3) + t231 * t245;
t165 = -rSges(3,3) * t231 + t230 * t245;
t164 = t213 * t264 + t261;
t163 = -t212 * t264 + t260;
t162 = -t212 * t231 + t228 * t260;
t161 = -t213 * t231 - t228 * t261;
t154 = -rSges(4,3) * t228 + (rSges(4,1) * t227 - rSges(4,2) * t225) * t226;
t153 = -Icges(4,5) * t228 + (Icges(4,1) * t227 - Icges(4,4) * t225) * t226;
t152 = -Icges(4,6) * t228 + (Icges(4,4) * t227 - Icges(4,2) * t225) * t226;
t151 = -Icges(4,3) * t228 + (Icges(4,5) * t227 - Icges(4,6) * t225) * t226;
t149 = t210 * t264 + t263;
t148 = -t209 * t264 + t262;
t147 = -t209 * t231 + t228 * t262;
t146 = -t210 * t231 - t228 * t263;
t144 = -rSges(5,3) * t228 + (rSges(5,1) * t213 - rSges(5,2) * t212) * t226;
t143 = -Icges(5,5) * t228 + (Icges(5,1) * t213 - Icges(5,4) * t212) * t226;
t142 = -Icges(5,6) * t228 + (Icges(5,4) * t213 - Icges(5,2) * t212) * t226;
t141 = -Icges(5,3) * t228 + (Icges(5,5) * t213 - Icges(5,6) * t212) * t226;
t140 = V_base(5) * rSges(2,3) - t201 * t215 + t249;
t139 = t203 * t215 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t137 = -rSges(6,3) * t228 + (rSges(6,1) * t210 - rSges(6,2) * t209) * t226;
t136 = -Icges(6,5) * t228 + (Icges(6,1) * t210 - Icges(6,4) * t209) * t226;
t135 = -Icges(6,6) * t228 + (Icges(6,4) * t210 - Icges(6,2) * t209) * t226;
t134 = -Icges(6,3) * t228 + (Icges(6,5) * t210 - Icges(6,6) * t209) * t226;
t133 = t201 * V_base(4) - t203 * V_base(5) + V_base(3);
t131 = (t223 + t229) * t228 + (t180 - t211) * t226;
t130 = t175 * rSges(4,1) + t174 * rSges(4,2) + rSges(4,3) * t265;
t129 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t266;
t128 = Icges(4,1) * t175 + Icges(4,4) * t174 + Icges(4,5) * t265;
t127 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t266;
t126 = Icges(4,4) * t175 + Icges(4,2) * t174 + Icges(4,6) * t265;
t125 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t266;
t124 = Icges(4,5) * t175 + Icges(4,6) * t174 + Icges(4,3) * t265;
t123 = Icges(4,5) * t173 + Icges(4,6) * t172 + Icges(4,3) * t266;
t120 = t164 * rSges(5,1) + t163 * rSges(5,2) + rSges(5,3) * t265;
t119 = rSges(5,1) * t162 + rSges(5,2) * t161 + rSges(5,3) * t266;
t118 = Icges(5,1) * t164 + Icges(5,4) * t163 + Icges(5,5) * t265;
t117 = Icges(5,1) * t162 + Icges(5,4) * t161 + Icges(5,5) * t266;
t116 = Icges(5,4) * t164 + Icges(5,2) * t163 + Icges(5,6) * t265;
t115 = Icges(5,4) * t162 + Icges(5,2) * t161 + Icges(5,6) * t266;
t114 = Icges(5,5) * t164 + Icges(5,6) * t163 + Icges(5,3) * t265;
t113 = Icges(5,5) * t162 + Icges(5,6) * t161 + Icges(5,3) * t266;
t110 = t149 * rSges(6,1) + t148 * rSges(6,2) + rSges(6,3) * t265;
t109 = rSges(6,1) * t147 + rSges(6,2) * t146 + rSges(6,3) * t266;
t108 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t265;
t107 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t266;
t106 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t265;
t105 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t266;
t104 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t265;
t103 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t266;
t102 = t192 * V_base(5) + (-t165 - t200) * t215 + t246;
t101 = t215 * t166 + (-pkin(5) - t192) * V_base(4) + t241;
t100 = t165 * V_base(4) + (-t166 - t202) * V_base(5) + t250;
t99 = t230 * t254 + t231 * t239 + t202;
t98 = t220 + (-qJ(2) - t254) * t231 + t239 * t230;
t97 = t154 * V_base(5) + (-t129 + t257) * t215 + t240;
t96 = t215 * t130 + (-t154 + t270) * V_base(4) + t237;
t95 = t129 * V_base(4) + (-t130 + t256) * V_base(5) + t238;
t94 = -t119 * t193 + t144 * t186 + t235;
t93 = t193 * t120 - t187 * t144 + t233;
t92 = t119 * t187 - t120 * t186 + t234;
t91 = -t109 * t176 + t131 * t186 + t137 * t169 - t193 * t98 + t235;
t90 = t176 * t110 - t187 * t131 - t170 * t137 + t193 * t99 + t233;
t89 = t109 * t170 - t110 * t169 - t186 * t99 + t187 * t98 + t234;
t1 = m(1) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t187 * ((t114 * t265 + t163 * t116 + t164 * t118) * t187 + (t113 * t265 + t163 * t115 + t164 * t117) * t186 + (t141 * t265 + t163 * t142 + t164 * t143) * t193) / 0.2e1 + t186 * ((t114 * t266 + t116 * t161 + t118 * t162) * t187 + (t113 * t266 + t115 * t161 + t117 * t162) * t186 + (t141 * t266 + t142 * t161 + t143 * t162) * t193) / 0.2e1 + t193 * ((-t113 * t186 - t114 * t187 - t141 * t193) * t228 + ((-t116 * t212 + t118 * t213) * t187 + (-t115 * t212 + t117 * t213) * t186 + (-t142 * t212 + t143 * t213) * t193) * t226) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t170 * ((t104 * t265 + t148 * t106 + t149 * t108) * t170 + (t103 * t265 + t148 * t105 + t149 * t107) * t169 + (t134 * t265 + t148 * t135 + t149 * t136) * t176) / 0.2e1 + t169 * ((t104 * t266 + t106 * t146 + t108 * t147) * t170 + (t103 * t266 + t105 * t146 + t107 * t147) * t169 + (t134 * t266 + t135 * t146 + t136 * t147) * t176) / 0.2e1 + t176 * ((-t103 * t169 - t104 * t170 - t134 * t176) * t228 + ((-t106 * t209 + t108 * t210) * t170 + (-t105 * t209 + t107 * t210) * t169 + (-t135 * t209 + t136 * t210) * t176) * t226) / 0.2e1 + ((t194 + (t157 - t123) * t228 + (-t125 * t225 + t127 * t227 + t159) * t226) * V_base(5) + (t195 + (t158 - t124) * t228 + (-t126 * t225 + t128 * t227 + t160) * t226) * V_base(4) + (Icges(2,3) + (t189 - t151) * t228 + (-t152 * t225 + t153 * t227 + t190) * t226) * t215) * t215 / 0.2e1 + (t230 * t236 + t231 * t232 + (t151 * t265 + t174 * t152 + t175 * t153 + t195) * t215 + (t123 * t265 + t174 * t125 + t175 * t127 - t230 * t196 + t198 * t231 + Icges(1,4)) * V_base(5) + (t124 * t265 + t174 * t126 + t175 * t128 - t230 * t197 + t199 * t231 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t230 * t232 - t231 * t236 + (t151 * t266 + t152 * t172 + t153 * t173 + t194) * t215 + (t123 * t266 + t125 * t172 + t127 * t173 + t196 * t231 + t230 * t198 + Icges(1,2)) * V_base(5) + (t124 * t266 + t126 * t172 + t128 * t173 + t197 * t231 + t230 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
