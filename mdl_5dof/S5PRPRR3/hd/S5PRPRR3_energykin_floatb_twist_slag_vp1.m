% Calculate kinetic energy for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:33
% EndTime: 2019-12-05 15:46:35
% DurationCPUTime: 2.68s
% Computational Cost: add. (1442->304), mult. (1640->442), div. (0->0), fcn. (1534->10), ass. (0->149)
t276 = Icges(3,3) + Icges(4,3);
t209 = qJ(2) + pkin(9);
t203 = sin(t209);
t204 = cos(t209);
t215 = sin(qJ(2));
t217 = cos(qJ(2));
t275 = Icges(3,5) * t217 + Icges(4,5) * t204 - Icges(3,6) * t215 - Icges(4,6) * t203;
t211 = sin(pkin(8));
t212 = cos(pkin(8));
t259 = Icges(4,4) * t204;
t233 = -Icges(4,2) * t203 + t259;
t138 = -Icges(4,6) * t212 + t211 * t233;
t139 = Icges(4,6) * t211 + t212 * t233;
t260 = Icges(4,4) * t203;
t235 = Icges(4,1) * t204 - t260;
t140 = -Icges(4,5) * t212 + t211 * t235;
t141 = Icges(4,5) * t211 + t212 * t235;
t261 = Icges(3,4) * t217;
t234 = -Icges(3,2) * t215 + t261;
t154 = -Icges(3,6) * t212 + t211 * t234;
t155 = Icges(3,6) * t211 + t212 * t234;
t262 = Icges(3,4) * t215;
t236 = Icges(3,1) * t217 - t262;
t156 = -Icges(3,5) * t212 + t211 * t236;
t157 = Icges(3,5) * t211 + t212 * t236;
t170 = Icges(4,2) * t204 + t260;
t171 = Icges(4,1) * t203 + t259;
t194 = Icges(3,2) * t217 + t262;
t195 = Icges(3,1) * t215 + t261;
t197 = -qJD(2) * t212 + V_base(5);
t198 = qJD(2) * t211 + V_base(4);
t272 = (-t170 * t203 + t171 * t204 - t194 * t215 + t195 * t217) * V_base(6) + (-t139 * t203 + t141 * t204 - t155 * t215 + t157 * t217) * t198 + (-t138 * t203 + t140 * t204 - t154 * t215 + t156 * t217) * t197;
t271 = (Icges(3,5) * t215 + Icges(4,5) * t203 + Icges(3,6) * t217 + Icges(4,6) * t204) * V_base(6) + (t276 * t211 + t275 * t212) * t198 + (t275 * t211 - t276 * t212) * t197;
t268 = pkin(2) * t215;
t267 = pkin(2) * t217;
t216 = cos(qJ(4));
t266 = pkin(4) * t216;
t263 = Icges(2,4) * t211;
t258 = t203 * t211;
t257 = t203 * t212;
t210 = qJ(4) + qJ(5);
t206 = sin(t210);
t256 = t206 * t211;
t255 = t206 * t212;
t207 = cos(t210);
t254 = t207 * t211;
t253 = t207 * t212;
t214 = sin(qJ(4));
t252 = t211 * t214;
t251 = t211 * t216;
t250 = t212 * t214;
t249 = t212 * t216;
t134 = -qJ(3) * t212 + t211 * t267;
t191 = pkin(1) * t211 - pkin(5) * t212;
t248 = -t134 - t191;
t247 = qJD(4) * t203;
t246 = qJD(5) * t203;
t245 = V_base(5) * qJ(1) + V_base(1);
t241 = qJD(1) + V_base(3);
t163 = t212 * t247 + t198;
t240 = qJD(3) * t211 + t197 * t268 + t245;
t239 = pkin(3) * t204 + pkin(6) * t203;
t238 = rSges(3,1) * t217 - rSges(3,2) * t215;
t237 = rSges(4,1) * t204 - rSges(4,2) * t203;
t162 = t211 * t247 + t197;
t192 = pkin(1) * t212 + pkin(5) * t211;
t230 = -V_base(4) * qJ(1) + V_base(6) * t192 + V_base(2);
t229 = V_base(4) * t191 - t192 * V_base(5) + t241;
t228 = pkin(7) * t203 + t204 * t266;
t227 = t198 * t134 + t229;
t135 = qJ(3) * t211 + t212 * t267;
t224 = -qJD(3) * t212 + V_base(6) * t135 + t230;
t160 = t239 * t211;
t173 = t203 * pkin(3) - t204 * pkin(6);
t223 = t197 * t173 + (-t160 + t248) * V_base(6) + t240;
t161 = t239 * t212;
t222 = t198 * t160 + (-t135 - t161) * t197 + t227;
t221 = V_base(6) * t161 + (-t173 - t268) * t198 + t224;
t205 = Icges(2,4) * t212;
t196 = rSges(3,1) * t215 + rSges(3,2) * t217;
t188 = rSges(2,1) * t212 - rSges(2,2) * t211;
t187 = rSges(2,1) * t211 + rSges(2,2) * t212;
t186 = -qJD(4) * t204 + V_base(6);
t185 = Icges(2,1) * t212 - t263;
t184 = Icges(2,1) * t211 + t205;
t183 = -Icges(2,2) * t211 + t205;
t182 = Icges(2,2) * t212 + t263;
t179 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t178 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t177 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = rSges(4,1) * t203 + rSges(4,2) * t204;
t168 = V_base(6) + (-qJD(4) - qJD(5)) * t204;
t167 = t204 * t249 + t252;
t166 = -t204 * t250 + t251;
t165 = t204 * t251 - t250;
t164 = -t204 * t252 - t249;
t159 = rSges(3,3) * t211 + t212 * t238;
t158 = -rSges(3,3) * t212 + t211 * t238;
t151 = t204 * t253 + t256;
t150 = -t204 * t255 + t254;
t149 = t204 * t254 - t255;
t148 = -t204 * t256 - t253;
t145 = rSges(4,3) * t211 + t212 * t237;
t144 = -rSges(4,3) * t212 + t211 * t237;
t143 = V_base(5) * rSges(2,3) - t187 * V_base(6) + t245;
t142 = t188 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t133 = -rSges(5,3) * t204 + (rSges(5,1) * t216 - rSges(5,2) * t214) * t203;
t132 = -Icges(5,5) * t204 + (Icges(5,1) * t216 - Icges(5,4) * t214) * t203;
t131 = -Icges(5,6) * t204 + (Icges(5,4) * t216 - Icges(5,2) * t214) * t203;
t130 = -Icges(5,3) * t204 + (Icges(5,5) * t216 - Icges(5,6) * t214) * t203;
t129 = t212 * t246 + t163;
t128 = t211 * t246 + t162;
t125 = -rSges(6,3) * t204 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t203;
t124 = -Icges(6,5) * t204 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t203;
t123 = -Icges(6,6) * t204 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t203;
t122 = -Icges(6,3) * t204 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t203;
t121 = t187 * V_base(4) - t188 * V_base(5) + t241;
t120 = -pkin(7) * t204 + t203 * t266;
t118 = rSges(5,1) * t167 + rSges(5,2) * t166 + rSges(5,3) * t257;
t117 = rSges(5,1) * t165 + rSges(5,2) * t164 + rSges(5,3) * t258;
t116 = Icges(5,1) * t167 + Icges(5,4) * t166 + Icges(5,5) * t257;
t115 = Icges(5,1) * t165 + Icges(5,4) * t164 + Icges(5,5) * t258;
t114 = Icges(5,4) * t167 + Icges(5,2) * t166 + Icges(5,6) * t257;
t113 = Icges(5,4) * t165 + Icges(5,2) * t164 + Icges(5,6) * t258;
t112 = Icges(5,5) * t167 + Icges(5,6) * t166 + Icges(5,3) * t257;
t111 = Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t258;
t110 = pkin(4) * t252 + t212 * t228;
t109 = -pkin(4) * t250 + t211 * t228;
t108 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t257;
t107 = rSges(6,1) * t149 + rSges(6,2) * t148 + rSges(6,3) * t258;
t106 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t257;
t105 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t258;
t104 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t257;
t103 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t258;
t102 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t257;
t101 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t258;
t100 = t196 * t197 + (-t158 - t191) * V_base(6) + t245;
t99 = t159 * V_base(6) - t196 * t198 + t230;
t98 = t158 * t198 - t159 * t197 + t229;
t97 = t172 * t197 + (-t144 + t248) * V_base(6) + t240;
t96 = t145 * V_base(6) + (-t172 - t268) * t198 + t224;
t95 = t144 * t198 + (-t135 - t145) * t197 + t227;
t94 = -t117 * t186 + t133 * t162 + t223;
t93 = t118 * t186 - t133 * t163 + t221;
t92 = t117 * t163 - t118 * t162 + t222;
t91 = -t107 * t168 - t109 * t186 + t120 * t162 + t125 * t128 + t223;
t90 = t108 * t168 + t110 * t186 - t120 * t163 - t125 * t129 + t221;
t89 = t107 * t129 - t108 * t128 + t109 * t163 - t110 * t162 + t222;
t1 = m(1) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(2) * (t121 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t163 * ((t112 * t257 + t166 * t114 + t167 * t116) * t163 + (t111 * t257 + t113 * t166 + t115 * t167) * t162 + (t130 * t257 + t131 * t166 + t132 * t167) * t186) / 0.2e1 + t162 * ((t112 * t258 + t114 * t164 + t116 * t165) * t163 + (t111 * t258 + t164 * t113 + t165 * t115) * t162 + (t130 * t258 + t131 * t164 + t132 * t165) * t186) / 0.2e1 + t186 * ((-t111 * t162 - t112 * t163 - t130 * t186) * t204 + ((-t114 * t214 + t116 * t216) * t163 + (-t113 * t214 + t115 * t216) * t162 + (-t131 * t214 + t132 * t216) * t186) * t203) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t129 * ((t102 * t257 + t150 * t104 + t151 * t106) * t129 + (t101 * t257 + t103 * t150 + t105 * t151) * t128 + (t122 * t257 + t123 * t150 + t124 * t151) * t168) / 0.2e1 + t128 * ((t102 * t258 + t104 * t148 + t106 * t149) * t129 + (t101 * t258 + t148 * t103 + t149 * t105) * t128 + (t122 * t258 + t123 * t148 + t124 * t149) * t168) / 0.2e1 + t168 * ((-t101 * t128 - t102 * t129 - t122 * t168) * t204 + ((-t104 * t206 + t106 * t207) * t129 + (-t103 * t206 + t105 * t207) * t128 + (-t123 * t206 + t124 * t207) * t168) * t203) / 0.2e1 + (t272 * t211 - t271 * t212) * t197 / 0.2e1 + (t271 * t211 + t272 * t212) * t198 / 0.2e1 + ((-t182 * t211 + t184 * t212 + Icges(1,4)) * V_base(5) + (-t211 * t183 + t212 * t185 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t212 * t182 + t211 * t184 + Icges(1,2)) * V_base(5) + (t183 * t212 + t185 * t211 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t139 * t204 + t141 * t203 + t155 * t217 + t157 * t215) * t198 + (t138 * t204 + t140 * t203 + t154 * t217 + t156 * t215) * t197 + (t204 * t170 + t203 * t171 + t217 * t194 + t215 * t195 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t212 - Icges(2,6) * t211 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t211 + Icges(2,6) * t212 + Icges(1,6));
T = t1;
