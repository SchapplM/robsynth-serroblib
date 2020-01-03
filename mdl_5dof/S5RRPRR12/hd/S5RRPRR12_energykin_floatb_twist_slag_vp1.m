% Calculate kinetic energy for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:57
% EndTime: 2019-12-31 20:29:00
% DurationCPUTime: 3.26s
% Computational Cost: add. (983->274), mult. (1985->396), div. (0->0), fcn. (2057->8), ass. (0->130)
t272 = Icges(3,4) - Icges(4,5);
t271 = Icges(3,1) + Icges(4,1);
t270 = Icges(3,2) + Icges(4,3);
t209 = sin(qJ(2));
t269 = t272 * t209;
t212 = cos(qJ(2));
t268 = t272 * t212;
t267 = Icges(4,4) + Icges(3,5);
t266 = Icges(3,6) - Icges(4,6);
t265 = t270 * t209 - t268;
t264 = t271 * t212 - t269;
t263 = Icges(4,2) + Icges(3,3);
t210 = sin(qJ(1));
t213 = cos(qJ(1));
t262 = t265 * t210 + t266 * t213;
t261 = -t266 * t210 + t265 * t213;
t260 = t264 * t210 - t267 * t213;
t259 = t267 * t210 + t264 * t213;
t258 = -t270 * t212 - t269;
t257 = t271 * t209 + t268;
t256 = -t266 * t209 + t267 * t212;
t198 = -qJD(2) * t213 + V_base(5);
t199 = qJD(2) * t210 + V_base(4);
t202 = V_base(6) + qJD(1);
t255 = (t258 * t209 + t257 * t212) * t202 + (t261 * t209 + t259 * t212) * t199 + (t262 * t209 + t260 * t212) * t198;
t254 = (t267 * t209 + t266 * t212) * t202 + (t263 * t210 + t256 * t213) * t199 + (t256 * t210 - t263 * t213) * t198;
t208 = sin(qJ(4));
t250 = cos(qJ(4));
t235 = t209 * t250;
t169 = -t212 * t208 + t235;
t249 = pkin(3) * t209;
t248 = Icges(2,4) * t210;
t242 = t212 * t213;
t231 = pkin(2) * t212 + qJ(3) * t209;
t163 = t231 * t210;
t194 = t210 * pkin(1) - pkin(6) * t213;
t241 = -t163 - t194;
t240 = qJD(3) * t209;
t239 = V_base(5) * pkin(5) + V_base(1);
t172 = t210 * t212 * pkin(3) + pkin(7) * t213;
t236 = -t172 + t241;
t189 = pkin(2) * t209 - qJ(3) * t212;
t234 = t198 * t189 + t213 * t240 + t239;
t233 = rSges(3,1) * t212 - rSges(3,2) * t209;
t232 = rSges(4,1) * t212 + rSges(4,3) * t209;
t224 = t198 * t249 + t234;
t166 = qJD(4) * t213 + t198;
t167 = -qJD(4) * t210 + t199;
t195 = pkin(1) * t213 + t210 * pkin(6);
t223 = -V_base(4) * pkin(5) + t202 * t195 + V_base(2);
t168 = t209 * t208 + t212 * t250;
t222 = V_base(4) * t194 - t195 * V_base(5) + V_base(3);
t164 = t231 * t213;
t219 = t202 * t164 + t210 * t240 + t223;
t218 = -qJD(3) * t212 + t199 * t163 + t222;
t173 = pkin(3) * t242 - t210 * pkin(7);
t217 = t202 * t173 + (-t189 - t249) * t199 + t219;
t216 = t199 * t172 + (-t164 - t173) * t198 + t218;
t211 = cos(qJ(5));
t207 = sin(qJ(5));
t205 = Icges(2,4) * t213;
t193 = rSges(2,1) * t213 - t210 * rSges(2,2);
t192 = t210 * rSges(2,1) + rSges(2,2) * t213;
t191 = rSges(3,1) * t209 + rSges(3,2) * t212;
t190 = rSges(4,1) * t209 - rSges(4,3) * t212;
t188 = Icges(2,1) * t213 - t248;
t187 = Icges(2,1) * t210 + t205;
t184 = -Icges(2,2) * t210 + t205;
t183 = Icges(2,2) * t213 + t248;
t176 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t175 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t174 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = t168 * t213;
t160 = t208 * t242 - t213 * t235;
t159 = t168 * t210;
t158 = t169 * t210;
t156 = t210 * rSges(3,3) + t213 * t233;
t155 = t210 * rSges(4,2) + t213 * t232;
t154 = -rSges(3,3) * t213 + t210 * t233;
t153 = -rSges(4,2) * t213 + t210 * t232;
t152 = qJD(5) * t168 + t202;
t136 = V_base(5) * rSges(2,3) - t192 * t202 + t239;
t135 = t193 * t202 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t134 = t161 * t211 - t207 * t210;
t133 = -t161 * t207 - t210 * t211;
t132 = t159 * t211 + t207 * t213;
t131 = -t159 * t207 + t211 * t213;
t130 = t192 * V_base(4) - t193 * V_base(5) + V_base(3);
t129 = pkin(4) * t169 + pkin(8) * t168;
t128 = rSges(5,1) * t169 - rSges(5,2) * t168;
t127 = Icges(5,1) * t169 - Icges(5,4) * t168;
t126 = Icges(5,4) * t169 - Icges(5,2) * t168;
t125 = Icges(5,5) * t169 - Icges(5,6) * t168;
t124 = qJD(5) * t160 + t167;
t123 = -qJD(5) * t158 + t166;
t122 = pkin(4) * t161 + pkin(8) * t160;
t121 = pkin(4) * t159 - pkin(8) * t158;
t120 = rSges(5,1) * t161 - rSges(5,2) * t160 - rSges(5,3) * t210;
t119 = t159 * rSges(5,1) + t158 * rSges(5,2) + rSges(5,3) * t213;
t118 = Icges(5,1) * t161 - Icges(5,4) * t160 - Icges(5,5) * t210;
t117 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t213;
t116 = Icges(5,4) * t161 - Icges(5,2) * t160 - Icges(5,6) * t210;
t115 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t213;
t114 = Icges(5,5) * t161 - Icges(5,6) * t160 - Icges(5,3) * t210;
t113 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t213;
t112 = rSges(6,3) * t168 + (rSges(6,1) * t211 - rSges(6,2) * t207) * t169;
t111 = Icges(6,5) * t168 + (Icges(6,1) * t211 - Icges(6,4) * t207) * t169;
t110 = Icges(6,6) * t168 + (Icges(6,4) * t211 - Icges(6,2) * t207) * t169;
t109 = Icges(6,3) * t168 + (Icges(6,5) * t211 - Icges(6,6) * t207) * t169;
t108 = t191 * t198 + (-t154 - t194) * t202 + t239;
t107 = t156 * t202 - t191 * t199 + t223;
t106 = rSges(6,1) * t134 + rSges(6,2) * t133 + rSges(6,3) * t160;
t105 = rSges(6,1) * t132 + rSges(6,2) * t131 - rSges(6,3) * t158;
t104 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t160;
t103 = Icges(6,1) * t132 + Icges(6,4) * t131 - Icges(6,5) * t158;
t102 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t160;
t101 = Icges(6,4) * t132 + Icges(6,2) * t131 - Icges(6,6) * t158;
t100 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t160;
t99 = Icges(6,5) * t132 + Icges(6,6) * t131 - Icges(6,3) * t158;
t98 = t154 * t199 - t156 * t198 + t222;
t97 = t190 * t198 + (-t153 + t241) * t202 + t234;
t96 = t155 * t202 + (-t189 - t190) * t199 + t219;
t95 = t153 * t199 + (-t155 - t164) * t198 + t218;
t94 = t128 * t166 + (-t119 + t236) * t202 + t224;
t93 = t120 * t202 - t128 * t167 + t217;
t92 = t119 * t167 - t120 * t166 + t216;
t91 = -t105 * t152 + t112 * t123 + t129 * t166 + (-t121 + t236) * t202 + t224;
t90 = t106 * t152 - t112 * t124 + t122 * t202 - t129 * t167 + t217;
t89 = t105 * t124 - t106 * t123 + t121 * t167 - t122 * t166 + t216;
t1 = m(1) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(2) * (t130 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(3) * (t107 ^ 2 + t108 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t167 * ((-t114 * t210 - t116 * t160 + t118 * t161) * t167 + (-t113 * t210 - t115 * t160 + t117 * t161) * t166 + (-t125 * t210 - t126 * t160 + t127 * t161) * t202) / 0.2e1 + t166 * ((t114 * t213 + t158 * t116 + t159 * t118) * t167 + (t113 * t213 + t158 * t115 + t159 * t117) * t166 + (t125 * t213 + t158 * t126 + t159 * t127) * t202) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t124 * ((t100 * t160 + t102 * t133 + t104 * t134) * t124 + (t101 * t133 + t103 * t134 + t160 * t99) * t123 + (t109 * t160 + t110 * t133 + t111 * t134) * t152) / 0.2e1 + t123 * ((-t100 * t158 + t102 * t131 + t104 * t132) * t124 + (t101 * t131 + t103 * t132 - t158 * t99) * t123 + (-t109 * t158 + t110 * t131 + t111 * t132) * t152) / 0.2e1 + t152 * ((t100 * t124 + t109 * t152 + t99 * t123) * t168 + ((-t102 * t207 + t104 * t211) * t124 + (-t101 * t207 + t103 * t211) * t123 + (-t110 * t207 + t111 * t211) * t152) * t169) / 0.2e1 + (t255 * t210 - t254 * t213) * t198 / 0.2e1 + (t254 * t210 + t255 * t213) * t199 / 0.2e1 + ((-t210 * t183 + t187 * t213 + Icges(1,4)) * V_base(5) + (-t210 * t184 + t188 * t213 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t183 * t213 + t210 * t187 + Icges(1,2)) * V_base(5) + (t184 * t213 + t210 * t188 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t116 * t168 + t118 * t169) * t167 + (-t115 * t168 + t117 * t169) * t166 + (t259 * t209 - t261 * t212) * t199 + (t260 * t209 - t262 * t212) * t198 + (-t126 * t168 + t127 * t169 + t257 * t209 - t258 * t212 + Icges(2,3)) * t202) * t202 / 0.2e1 + t202 * V_base(4) * (Icges(2,5) * t213 - Icges(2,6) * t210) + V_base(5) * t202 * (Icges(2,5) * t210 + Icges(2,6) * t213) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
