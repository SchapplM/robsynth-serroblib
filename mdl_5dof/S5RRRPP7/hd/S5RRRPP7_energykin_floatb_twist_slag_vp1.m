% Calculate kinetic energy for
% S5RRRPP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:53
% EndTime: 2019-12-31 21:03:55
% DurationCPUTime: 1.99s
% Computational Cost: add. (953->229), mult. (1883->327), div. (0->0), fcn. (1857->6), ass. (0->110)
t246 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t245 = -Icges(4,4) + Icges(6,4) + Icges(5,5);
t244 = Icges(6,5) - Icges(5,4) - Icges(4,5);
t243 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t242 = -Icges(5,6) + Icges(6,6) + Icges(4,6);
t241 = Icges(6,3) + Icges(4,3) + Icges(5,2);
t240 = rSges(6,1) + pkin(4);
t239 = rSges(6,3) + qJ(5);
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t195 = cos(qJ(1));
t192 = sin(qJ(1));
t194 = cos(qJ(2));
t220 = t192 * t194;
t153 = t190 * t220 + t193 * t195;
t154 = -t190 * t195 + t193 * t220;
t191 = sin(qJ(2));
t222 = t191 * t192;
t238 = -t153 * t242 - t154 * t244 + t222 * t241;
t219 = t194 * t195;
t155 = t190 * t219 - t192 * t193;
t156 = t190 * t192 + t193 * t219;
t221 = t191 * t195;
t237 = -t155 * t242 - t156 * t244 + t221 * t241;
t236 = t153 * t243 + t154 * t245 - t222 * t242;
t235 = t155 * t243 + t156 * t245 - t221 * t242;
t234 = t153 * t245 + t154 * t246 - t222 * t244;
t233 = t155 * t245 + t156 * t246 - t221 * t244;
t232 = t241 * t194 + (t190 * t242 + t193 * t244) * t191;
t231 = t242 * t194 + (t190 * t243 + t193 * t245) * t191;
t230 = t244 * t194 + (t190 * t245 + t193 * t246) * t191;
t225 = Icges(2,4) * t192;
t224 = Icges(3,4) * t191;
t223 = Icges(3,4) * t194;
t218 = rSges(6,2) * t153 + t154 * t240 - t222 * t239;
t217 = t155 * rSges(6,2) + t156 * t240 - t221 * t239;
t216 = t239 * t194 + (rSges(6,2) * t190 + t193 * t240) * t191;
t215 = qJD(3) * t191;
t214 = qJD(5) * t191;
t213 = V_base(5) * pkin(5) + V_base(1);
t183 = qJD(2) * t192 + V_base(4);
t186 = V_base(6) + qJD(1);
t210 = pkin(2) * t194 + pkin(7) * t191;
t182 = -qJD(2) * t195 + V_base(5);
t209 = rSges(3,1) * t194 - rSges(3,2) * t191;
t208 = Icges(3,1) * t194 - t224;
t207 = -Icges(3,2) * t191 + t223;
t206 = Icges(3,5) * t194 - Icges(3,6) * t191;
t181 = pkin(1) * t195 + pkin(6) * t192;
t205 = -V_base(4) * pkin(5) + t186 * t181 + V_base(2);
t180 = pkin(1) * t192 - pkin(6) * t195;
t204 = V_base(4) * t180 - t181 * V_base(5) + V_base(3);
t159 = t210 * t192;
t179 = pkin(2) * t191 - pkin(7) * t194;
t203 = t182 * t179 + (-t159 - t180) * t186 + t213;
t202 = (-Icges(3,3) * t195 + t192 * t206) * t182 + (Icges(3,3) * t192 + t195 * t206) * t183 + (Icges(3,5) * t191 + Icges(3,6) * t194) * t186;
t160 = t210 * t195;
t201 = t186 * t160 - t179 * t183 + t205;
t151 = t192 * t215 + t182;
t157 = (pkin(3) * t193 + qJ(4) * t190) * t191;
t200 = qJD(4) * t155 + t151 * t157 + t203;
t199 = t183 * t159 - t160 * t182 + t204;
t119 = pkin(3) * t156 + qJ(4) * t155;
t175 = -qJD(3) * t194 + t186;
t198 = qJD(4) * t153 + t175 * t119 + t201;
t118 = pkin(3) * t154 + qJ(4) * t153;
t152 = t195 * t215 + t183;
t197 = qJD(4) * t191 * t190 + t152 * t118 + t199;
t136 = -Icges(3,6) * t195 + t192 * t207;
t137 = Icges(3,6) * t192 + t195 * t207;
t141 = -Icges(3,5) * t195 + t192 * t208;
t142 = Icges(3,5) * t192 + t195 * t208;
t169 = Icges(3,2) * t194 + t224;
t172 = Icges(3,1) * t191 + t223;
t196 = (-t137 * t191 + t142 * t194) * t183 + (-t136 * t191 + t141 * t194) * t182 + (-t169 * t191 + t172 * t194) * t186;
t188 = Icges(2,4) * t195;
t178 = rSges(2,1) * t195 - rSges(2,2) * t192;
t177 = rSges(2,1) * t192 + rSges(2,2) * t195;
t176 = rSges(3,1) * t191 + rSges(3,2) * t194;
t174 = Icges(2,1) * t195 - t225;
t173 = Icges(2,1) * t192 + t188;
t171 = -Icges(2,2) * t192 + t188;
t170 = Icges(2,2) * t195 + t225;
t165 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t164 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t163 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t147 = t192 * rSges(3,3) + t195 * t209;
t146 = -rSges(3,3) * t195 + t192 * t209;
t145 = -rSges(4,3) * t194 + (rSges(4,1) * t193 - rSges(4,2) * t190) * t191;
t144 = -rSges(5,2) * t194 + (rSges(5,1) * t193 + rSges(5,3) * t190) * t191;
t123 = V_base(5) * rSges(2,3) - t177 * t186 + t213;
t122 = t178 * t186 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t121 = t177 * V_base(4) - t178 * V_base(5) + V_base(3);
t117 = rSges(4,1) * t156 - rSges(4,2) * t155 + rSges(4,3) * t221;
t116 = rSges(5,1) * t156 + rSges(5,2) * t221 + rSges(5,3) * t155;
t114 = rSges(4,1) * t154 - rSges(4,2) * t153 + rSges(4,3) * t222;
t113 = rSges(5,1) * t154 + rSges(5,2) * t222 + rSges(5,3) * t153;
t91 = t176 * t182 + (-t146 - t180) * t186 + t213;
t90 = t147 * t186 - t176 * t183 + t205;
t89 = t146 * t183 - t147 * t182 + t204;
t88 = -t114 * t175 + t145 * t151 + t203;
t87 = t117 * t175 - t145 * t152 + t201;
t86 = t114 * t152 - t117 * t151 + t199;
t85 = t144 * t151 + (-t113 - t118) * t175 + t200;
t84 = t116 * t175 + (-t144 - t157) * t152 + t198;
t83 = t113 * t152 + (-t116 - t119) * t151 + t197;
t82 = -t195 * t214 + t216 * t151 + (-t118 - t218) * t175 + t200;
t81 = -t192 * t214 + t217 * t175 + (-t157 - t216) * t152 + t198;
t80 = qJD(5) * t194 + t218 * t152 + (-t119 - t217) * t151 + t197;
t1 = m(1) * (t163 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(2) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t183 * (t202 * t192 + t196 * t195) / 0.2e1 + t182 * (t196 * t192 - t202 * t195) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + ((t137 * t194 + t142 * t191) * t183 + (t136 * t194 + t141 * t191) * t182 + (t169 * t194 + t172 * t191 + Icges(2,3)) * t186) * t186 / 0.2e1 + ((-t170 * t192 + t173 * t195 + Icges(1,4)) * V_base(5) + (-t192 * t171 + t174 * t195 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t170 * t195 + t192 * t173 + Icges(1,2)) * V_base(5) + (t171 * t195 + t174 * t192 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t153 * t231 + t154 * t230 - t222 * t232) * t175 + (t153 * t235 + t154 * t233 + t222 * t237) * t152 + (t236 * t153 + t234 * t154 + t238 * t222) * t151) * t151 / 0.2e1 + ((t155 * t231 + t156 * t230 - t221 * t232) * t175 + (t235 * t155 + t233 * t156 + t237 * t221) * t152 + (t155 * t236 + t156 * t234 + t221 * t238) * t151) * t152 / 0.2e1 + ((-t151 * t238 - t152 * t237 + t175 * t232) * t194 + ((t190 * t231 + t193 * t230) * t175 + (t190 * t235 + t193 * t233) * t152 + (t190 * t236 + t193 * t234) * t151) * t191) * t175 / 0.2e1 + V_base(4) * t186 * (Icges(2,5) * t195 - Icges(2,6) * t192) + V_base(5) * t186 * (Icges(2,5) * t192 + Icges(2,6) * t195) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
