% Calculate kinetic energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:17
% DurationCPUTime: 2.31s
% Computational Cost: add. (1050->272), mult. (1884->354), div. (0->0), fcn. (2056->8), ass. (0->134)
t278 = Icges(2,4) - Icges(3,5);
t277 = Icges(4,4) + Icges(5,6);
t276 = Icges(2,1) + Icges(3,1);
t275 = Icges(4,1) + Icges(5,2);
t274 = Icges(3,4) + Icges(2,5);
t273 = -Icges(5,4) + Icges(4,5);
t272 = Icges(5,5) - Icges(4,6);
t271 = Icges(2,2) + Icges(3,3);
t270 = Icges(4,2) + Icges(5,3);
t269 = Icges(2,6) - Icges(3,6);
t247 = sin(pkin(9));
t248 = cos(pkin(9));
t251 = sin(qJ(1));
t252 = cos(qJ(1));
t162 = t247 * t252 - t248 * t251;
t268 = t277 * t162;
t267 = t278 * t251;
t266 = t278 * t252;
t161 = -t247 * t251 - t248 * t252;
t265 = t277 * t161;
t264 = t162 * t270 + t265;
t263 = -t161 * t270 + t268;
t262 = t161 * t275 + t268;
t261 = t162 * t275 - t265;
t260 = -t252 * t271 - t267;
t259 = t251 * t271 - t266;
t258 = t251 * t276 + t266;
t257 = t252 * t276 - t267;
t206 = cos(qJ(5));
t204 = sin(qJ(5));
t245 = Icges(6,4) * t204;
t221 = Icges(6,2) * t206 + t245;
t107 = -Icges(6,6) * t161 + t162 * t221;
t108 = -Icges(6,6) * t162 - t161 * t221;
t244 = Icges(6,4) * t206;
t222 = Icges(6,1) * t204 + t244;
t109 = -Icges(6,5) * t161 + t162 * t222;
t110 = -Icges(6,5) * t162 - t161 * t222;
t151 = -qJD(5) * t162 + V_base(5);
t152 = -qJD(5) * t161 + V_base(4);
t172 = Icges(6,2) * t204 - t244;
t177 = -Icges(6,1) * t206 + t245;
t196 = V_base(6) + qJD(1);
t254 = (t107 * t206 + t109 * t204) * t152 + (t108 * t206 + t110 * t204) * t151 + (t172 * t206 + t177 * t204) * t196;
t250 = pkin(7) * t161;
t249 = pkin(7) * t162;
t242 = t161 * t206;
t241 = t162 * t206;
t203 = sin(qJ(6));
t240 = t203 * t204;
t205 = cos(qJ(6));
t239 = t204 * t205;
t238 = qJD(6) * t206;
t183 = pkin(1) * t251 - qJ(2) * t252;
t237 = V_base(4) * t183 + V_base(3);
t236 = V_base(5) * pkin(6) + V_base(1);
t233 = t252 * pkin(2);
t232 = t251 * pkin(2);
t229 = qJD(2) * t251 + t236;
t228 = -t183 - t232;
t187 = pkin(1) * t252 + qJ(2) * t251;
t227 = -t187 - t233;
t226 = qJD(4) * t162 + t229;
t225 = pkin(5) * t204 - pkin(8) * t206;
t224 = V_base(4) * t232 - qJD(3) + t237;
t223 = rSges(6,1) * t204 + rSges(6,2) * t206;
t220 = Icges(6,5) * t204 + Icges(6,6) * t206;
t139 = -pkin(3) * t162 - qJ(4) * t161;
t216 = -t139 + t228;
t141 = -pkin(3) * t161 + qJ(4) * t162;
t215 = -t141 + t227;
t214 = V_base(4) * t139 + t224;
t213 = -qJD(2) * t252 + t196 * t187 + V_base(2);
t212 = -(-Icges(6,3) * t161 + t162 * t220) * t152 - (-Icges(6,3) * t162 - t161 * t220) * t151 - (-Icges(6,5) * t206 + Icges(6,6) * t204) * t196;
t211 = V_base(4) * qJ(3) + t196 * t233 + t213;
t210 = -qJD(4) * t161 + t196 * t141 + t211;
t209 = t196 * t249 + (-pkin(4) - qJ(3)) * V_base(5) + t226;
t208 = t210 + (pkin(4) - pkin(6)) * V_base(4);
t207 = -t249 * V_base(4) + t214 + (t215 + t250) * V_base(5);
t190 = -pkin(5) * t206 - t204 * pkin(8);
t189 = rSges(2,1) * t252 - rSges(2,2) * t251;
t188 = rSges(3,1) * t252 + rSges(3,3) * t251;
t186 = -rSges(6,1) * t206 + t204 * rSges(6,2);
t185 = rSges(2,1) * t251 + rSges(2,2) * t252;
t184 = rSges(3,1) * t251 - rSges(3,3) * t252;
t182 = -qJD(6) * t204 + t196;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = -t204 * rSges(7,3) + (-rSges(7,1) * t205 + rSges(7,2) * t203) * t206;
t148 = -Icges(7,5) * t204 + (-Icges(7,1) * t205 + Icges(7,4) * t203) * t206;
t147 = -Icges(7,6) * t204 + (-Icges(7,4) * t205 + Icges(7,2) * t203) * t206;
t146 = -Icges(7,3) * t204 + (-Icges(7,5) * t205 + Icges(7,6) * t203) * t206;
t145 = V_base(5) * rSges(2,3) - t185 * t196 + t236;
t144 = t189 * t196 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t143 = t185 * V_base(4) - t189 * V_base(5) + V_base(3);
t142 = -rSges(4,1) * t161 - rSges(4,2) * t162;
t140 = -rSges(4,1) * t162 + rSges(4,2) * t161;
t138 = rSges(5,2) * t162 - rSges(5,3) * t161;
t137 = rSges(5,2) * t161 + rSges(5,3) * t162;
t123 = t225 * t161;
t122 = t225 * t162;
t120 = -t161 * t239 - t162 * t203;
t119 = t161 * t240 - t162 * t205;
t118 = -t161 * t203 + t162 * t239;
t117 = -t161 * t205 - t162 * t240;
t116 = -t162 * t238 + t152;
t115 = t161 * t238 + t151;
t114 = V_base(5) * rSges(3,2) + (-t183 - t184) * t196 + t229;
t113 = t196 * t188 + (-pkin(6) - rSges(3,2)) * V_base(4) + t213;
t112 = -t162 * rSges(6,3) - t161 * t223;
t111 = -t161 * rSges(6,3) + t162 * t223;
t104 = t184 * V_base(4) + (-t187 - t188) * V_base(5) + t237;
t103 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t140 + t228) * t196 + t229;
t102 = t196 * t142 + (rSges(4,3) - pkin(6)) * V_base(4) + t211;
t101 = t120 * rSges(7,1) + t119 * rSges(7,2) + rSges(7,3) * t242;
t100 = t118 * rSges(7,1) + t117 * rSges(7,2) - rSges(7,3) * t241;
t99 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t242;
t98 = Icges(7,1) * t118 + Icges(7,4) * t117 - Icges(7,5) * t241;
t97 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t242;
t96 = Icges(7,4) * t118 + Icges(7,2) * t117 - Icges(7,6) * t241;
t95 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t242;
t94 = Icges(7,5) * t118 + Icges(7,6) * t117 - Icges(7,3) * t241;
t93 = V_base(4) * t140 + (-t142 + t227) * V_base(5) + t224;
t92 = (-qJ(3) - rSges(5,1)) * V_base(5) + (-t138 + t216) * t196 + t226;
t91 = t196 * t137 + (rSges(5,1) - pkin(6)) * V_base(4) + t210;
t90 = V_base(4) * t138 + (-t137 + t215) * V_base(5) + t214;
t89 = t151 * t186 + (-t112 + t216) * t196 + t209;
t88 = -t152 * t186 + (t111 - t250) * t196 + t208;
t87 = -t151 * t111 + t152 * t112 + t207;
t86 = -t182 * t101 + t115 * t150 + t151 * t190 + (t123 + t216) * t196 + t209;
t85 = t182 * t100 - t116 * t150 - t152 * t190 + (t122 - t250) * t196 + t208;
t84 = -t115 * t100 + t116 * t101 - t151 * t122 - t152 * t123 + t207;
t1 = t152 * (t212 * t161 + t254 * t162) / 0.2e1 + t151 * (-t254 * t161 + t212 * t162) / 0.2e1 + m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(3) * (t104 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t182 * ((-t95 * t115 - t94 * t116 - t146 * t182) * t204 + ((t203 * t96 - t205 * t98) * t116 + (t203 * t97 - t205 * t99) * t115 + (t147 * t203 - t148 * t205) * t182) * t206) / 0.2e1 + t116 * ((t117 * t96 + t118 * t98 - t94 * t241) * t116 + (t117 * t97 + t118 * t99 - t241 * t95) * t115 + (t117 * t147 + t118 * t148 - t146 * t241) * t182) / 0.2e1 + t115 * ((t119 * t96 + t120 * t98 + t242 * t94) * t116 + (t119 * t97 + t120 * t99 + t95 * t242) * t115 + (t119 * t147 + t120 * t148 + t146 * t242) * t182) / 0.2e1 + ((t204 * t107 - t109 * t206) * t152 + (t204 * t108 - t110 * t206) * t151 + (t204 * t172 - t177 * t206 + Icges(5,1) + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t196) * t196 / 0.2e1 + ((t161 * t261 + t162 * t263 + t251 * t260 + t252 * t258 + Icges(1,4)) * V_base(5) + (t262 * t161 + t264 * t162 + t259 * t251 + t257 * t252 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t263 * t161 + t261 * t162 + t258 * t251 - t260 * t252 + Icges(1,2)) * V_base(5) + (-t161 * t264 + t262 * t162 + t257 * t251 - t259 * t252 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t196 * (t161 * t272 + t162 * t273 + t251 * t274 + t252 * t269) + V_base(4) * t196 * (t161 * t273 - t162 * t272 - t269 * t251 + t274 * t252) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
