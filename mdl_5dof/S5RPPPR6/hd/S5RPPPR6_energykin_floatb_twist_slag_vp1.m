% Calculate kinetic energy for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:33
% EndTime: 2019-12-31 17:47:36
% DurationCPUTime: 3.23s
% Computational Cost: add. (870->291), mult. (1817->385), div. (0->0), fcn. (1843->8), ass. (0->140)
t278 = Icges(3,4) + Icges(4,6);
t277 = Icges(3,1) + Icges(4,2);
t276 = Icges(3,2) + Icges(4,3);
t207 = cos(pkin(7));
t275 = t278 * t207;
t205 = sin(pkin(7));
t274 = t278 * t205;
t273 = Icges(4,4) - Icges(3,5);
t272 = Icges(4,5) - Icges(3,6);
t271 = t276 * t205 - t275;
t270 = t277 * t207 - t274;
t269 = Icges(4,1) + Icges(3,3);
t209 = sin(qJ(1));
t211 = cos(qJ(1));
t268 = t271 * t209 - t272 * t211;
t267 = t272 * t209 + t271 * t211;
t266 = t270 * t209 + t273 * t211;
t265 = -t273 * t209 + t270 * t211;
t264 = -t276 * t207 - t274;
t263 = t277 * t205 + t275;
t262 = t272 * t205 - t273 * t207;
t200 = V_base(6) + qJD(1);
t261 = (t268 * t205 + t266 * t207) * V_base(5) + (t267 * t205 + t265 * t207) * V_base(4) + (t264 * t205 + t263 * t207) * t200;
t260 = (t262 * t209 - t269 * t211) * V_base(5) + (t269 * t209 + t262 * t211) * V_base(4) + (-t273 * t205 - t272 * t207) * t200;
t180 = pkin(2) * t205 - qJ(3) * t207;
t258 = -pkin(5) - t180;
t257 = Icges(2,4) * t209;
t252 = qJ(4) * t205;
t204 = sin(pkin(8));
t251 = t204 * t209;
t250 = t204 * t211;
t206 = cos(pkin(8));
t249 = t206 * t207;
t248 = t206 * t211;
t247 = t207 * t209;
t210 = cos(qJ(5));
t246 = t207 * t210;
t245 = t207 * t211;
t244 = t209 * t206;
t229 = pkin(2) * t207 + qJ(3) * t205;
t163 = t229 * t209;
t189 = t209 * pkin(1) - qJ(2) * t211;
t243 = -t163 - t189;
t164 = t229 * t211;
t191 = pkin(1) * t211 + t209 * qJ(2);
t242 = -t164 - t191;
t241 = qJD(3) * t205;
t240 = qJD(4) * t207;
t239 = V_base(4) * t189 + V_base(3);
t238 = V_base(5) * pkin(5) + V_base(1);
t170 = -pkin(3) * t211 + qJ(4) * t247;
t235 = -t170 + t243;
t169 = t209 * pkin(3) + qJ(4) * t245;
t234 = -t169 + t242;
t233 = qJD(2) * t209 + t238;
t232 = -t252 + t258;
t231 = rSges(3,1) * t207 - rSges(3,2) * t205;
t230 = -rSges(4,2) * t207 + rSges(4,3) * t205;
t222 = -qJD(2) * t211 + t200 * t191 + V_base(2);
t221 = V_base(5) * t180 + t211 * t241 + t233;
t220 = -qJD(3) * t207 + V_base(4) * t163 + t239;
t219 = t200 * t164 + t209 * t241 + t222;
t218 = t211 * t240 + V_base(5) * t252 + t221;
t217 = qJD(4) * t205 + V_base(4) * t170 + t220;
t216 = t200 * t169 + t209 * t240 + t219;
t208 = sin(qJ(5));
t202 = Icges(2,4) * t211;
t192 = rSges(2,1) * t211 - t209 * rSges(2,2);
t190 = t209 * rSges(2,1) + rSges(2,2) * t211;
t188 = Icges(2,1) * t211 - t257;
t187 = Icges(2,1) * t209 + t202;
t186 = -Icges(2,2) * t209 + t202;
t185 = Icges(2,2) * t211 + t257;
t184 = Icges(2,5) * t211 - Icges(2,6) * t209;
t183 = Icges(2,5) * t209 + Icges(2,6) * t211;
t182 = rSges(3,1) * t205 + rSges(3,2) * t207;
t181 = -rSges(4,2) * t205 - rSges(4,3) * t207;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t167 = qJD(5) * t249 + t200;
t161 = (-pkin(4) * t204 + pkin(6) * t206) * t207;
t160 = t205 * t251 - t248;
t159 = t205 * t244 + t250;
t158 = t205 * t250 + t244;
t157 = -t205 * t248 + t251;
t156 = -t204 * t246 + t205 * t208;
t155 = t204 * t207 * t208 + t205 * t210;
t152 = -rSges(4,1) * t211 + t209 * t230;
t151 = t209 * rSges(4,1) + t211 * t230;
t150 = t209 * rSges(3,3) + t211 * t231;
t149 = -rSges(3,3) * t211 + t209 * t231;
t135 = rSges(5,3) * t205 + (-rSges(5,1) * t204 - rSges(5,2) * t206) * t207;
t134 = Icges(5,5) * t205 + (-Icges(5,1) * t204 - Icges(5,4) * t206) * t207;
t133 = Icges(5,6) * t205 + (-Icges(5,4) * t204 - Icges(5,2) * t206) * t207;
t132 = Icges(5,3) * t205 + (-Icges(5,5) * t204 - Icges(5,6) * t206) * t207;
t131 = qJD(5) * t157 + V_base(4);
t130 = -qJD(5) * t159 + V_base(5);
t129 = V_base(5) * rSges(2,3) - t190 * t200 + t238;
t128 = t192 * t200 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t127 = t190 * V_base(4) - t192 * V_base(5) + V_base(3);
t126 = t160 * t210 + t208 * t247;
t125 = -t160 * t208 + t209 * t246;
t124 = t158 * t210 + t208 * t245;
t123 = -t158 * t208 + t210 * t245;
t122 = pkin(4) * t160 - pkin(6) * t159;
t121 = pkin(4) * t158 + pkin(6) * t157;
t120 = rSges(5,1) * t160 + rSges(5,2) * t159 + rSges(5,3) * t247;
t119 = t158 * rSges(5,1) - t157 * rSges(5,2) + rSges(5,3) * t245;
t118 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t249;
t117 = Icges(5,1) * t160 + Icges(5,4) * t159 + Icges(5,5) * t247;
t116 = Icges(5,1) * t158 - Icges(5,4) * t157 + Icges(5,5) * t245;
t115 = Icges(5,4) * t160 + Icges(5,2) * t159 + Icges(5,6) * t247;
t114 = Icges(5,4) * t158 - Icges(5,2) * t157 + Icges(5,6) * t245;
t113 = Icges(5,5) * t160 + Icges(5,6) * t159 + Icges(5,3) * t247;
t112 = Icges(5,5) * t158 - Icges(5,6) * t157 + Icges(5,3) * t245;
t111 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t249;
t110 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t249;
t109 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t249;
t108 = t182 * V_base(5) + (-t149 - t189) * t200 + t233;
t107 = t200 * t150 + (-pkin(5) - t182) * V_base(4) + t222;
t106 = t149 * V_base(4) + (-t150 - t191) * V_base(5) + t239;
t105 = rSges(6,1) * t126 + rSges(6,2) * t125 - rSges(6,3) * t159;
t104 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t157;
t103 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t159;
t102 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t157;
t101 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t159;
t100 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t157;
t99 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t159;
t98 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t157;
t97 = t181 * V_base(5) + (-t152 + t243) * t200 + t221;
t96 = t200 * t151 + (-t181 + t258) * V_base(4) + t219;
t95 = t152 * V_base(4) + (-t151 + t242) * V_base(5) + t220;
t94 = t135 * V_base(5) + (-t120 + t235) * t200 + t218;
t93 = t200 * t119 + (-t135 + t232) * V_base(4) + t216;
t92 = t120 * V_base(4) + (-t119 + t234) * V_base(5) + t217;
t91 = -t105 * t167 + t118 * t130 + t161 * V_base(5) + (-t122 + t235) * t200 + t218;
t90 = t167 * t104 - t131 * t118 + t200 * t121 + (-t161 + t232) * V_base(4) + t216;
t89 = -t104 * t130 + t105 * t131 + t122 * V_base(4) + (-t121 + t234) * V_base(5) + t217;
t1 = m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t131 * ((t100 * t123 + t102 * t124 + t157 * t98) * t131 + (t101 * t123 + t103 * t124 + t157 * t99) * t130 + (t109 * t157 + t110 * t123 + t111 * t124) * t167) / 0.2e1 + t130 * ((t100 * t125 + t102 * t126 - t159 * t98) * t131 + (t101 * t125 + t103 * t126 - t159 * t99) * t130 + (-t109 * t159 + t110 * t125 + t111 * t126) * t167) / 0.2e1 + t167 * ((t100 * t155 + t102 * t156 + t249 * t98) * t131 + (t101 * t155 + t103 * t156 + t249 * t99) * t130 + (t109 * t249 + t110 * t155 + t111 * t156) * t167) / 0.2e1 + ((t183 + (-t115 * t206 - t117 * t204 - t268) * t207 + (t113 + t266) * t205) * V_base(5) + (t184 + (-t114 * t206 - t116 * t204 - t267) * t207 + (t112 + t265) * t205) * V_base(4) + (Icges(2,3) + (-t133 * t206 - t134 * t204 - t264) * t207 + (t132 + t263) * t205) * t200) * t200 / 0.2e1 + (t261 * t211 + t260 * t209 + (t132 * t245 - t157 * t133 + t158 * t134 + t184) * t200 + (t113 * t245 - t157 * t115 + t158 * t117 - t209 * t185 + t187 * t211 + Icges(1,4)) * V_base(5) + (t112 * t245 - t157 * t114 + t158 * t116 - t209 * t186 + t188 * t211 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (-t260 * t211 + t261 * t209 + (t132 * t247 + t133 * t159 + t134 * t160 + t183) * t200 + (t113 * t247 + t115 * t159 + t117 * t160 + t185 * t211 + t209 * t187 + Icges(1,2)) * V_base(5) + (t112 * t247 + t114 * t159 + t116 * t160 + t186 * t211 + t209 * t188 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
