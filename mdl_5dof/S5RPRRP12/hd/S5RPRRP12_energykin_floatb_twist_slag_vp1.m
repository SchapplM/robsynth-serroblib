% Calculate kinetic energy for
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:12
% EndTime: 2019-12-31 18:56:14
% DurationCPUTime: 2.28s
% Computational Cost: add. (767->235), mult. (1423->319), div. (0->0), fcn. (1317->6), ass. (0->118)
t260 = Icges(2,4) + Icges(3,6);
t259 = Icges(2,1) + Icges(3,2);
t258 = Icges(5,1) + Icges(6,1);
t257 = -Icges(3,4) + Icges(2,5);
t256 = Icges(5,4) + Icges(6,4);
t255 = Icges(3,5) - Icges(2,6);
t254 = Icges(6,5) + Icges(5,5);
t253 = Icges(2,2) + Icges(3,3);
t252 = Icges(5,2) + Icges(6,2);
t251 = Icges(6,6) + Icges(5,6);
t250 = Icges(6,3) + Icges(5,3);
t183 = cos(qJ(1));
t249 = t260 * t183;
t180 = sin(qJ(1));
t248 = t260 * t180;
t179 = sin(qJ(3));
t181 = cos(qJ(4));
t212 = t181 * t183;
t178 = sin(qJ(4));
t215 = t180 * t178;
t130 = -t179 * t215 + t212;
t214 = t180 * t181;
t216 = t178 * t183;
t131 = t179 * t214 + t216;
t182 = cos(qJ(3));
t213 = t180 * t182;
t247 = t251 * t130 + t254 * t131 - t250 * t213;
t132 = t179 * t216 + t214;
t133 = -t179 * t212 + t215;
t211 = t182 * t183;
t246 = t251 * t132 + t254 * t133 + t250 * t211;
t245 = t252 * t130 + t256 * t131 - t251 * t213;
t244 = t252 * t132 + t256 * t133 + t251 * t211;
t243 = t256 * t130 + t258 * t131 - t254 * t213;
t242 = t256 * t132 + t258 * t133 + t254 * t211;
t241 = (-t251 * t178 + t254 * t181) * t182 + t250 * t179;
t240 = (-t252 * t178 + t256 * t181) * t182 + t251 * t179;
t239 = (-t256 * t178 + t258 * t181) * t182 + t254 * t179;
t238 = -t253 * t183 - t248;
t237 = t253 * t180 - t249;
t236 = t259 * t180 + t249;
t235 = t259 * t183 - t248;
t224 = pkin(4) * t181;
t232 = -qJ(5) * t182 + t179 * t224;
t220 = Icges(4,4) * t179;
t194 = Icges(4,2) * t182 + t220;
t117 = Icges(4,6) * t183 + t180 * t194;
t118 = Icges(4,6) * t180 - t183 * t194;
t219 = Icges(4,4) * t182;
t195 = Icges(4,1) * t179 + t219;
t121 = Icges(4,5) * t183 + t180 * t195;
t122 = Icges(4,5) * t180 - t183 * t195;
t148 = -Icges(4,2) * t179 + t219;
t153 = Icges(4,1) * t182 - t220;
t166 = qJD(3) * t180 + V_base(5);
t167 = qJD(3) * t183 + V_base(4);
t170 = V_base(6) + qJD(1);
t231 = (t117 * t182 + t121 * t179) * t167 + (t118 * t182 + t122 * t179) * t166 + (t148 * t182 + t153 * t179) * t170;
t226 = pkin(6) * t180;
t225 = pkin(6) * t183;
t222 = rSges(6,1) * t131 + rSges(6,2) * t130 - rSges(6,3) * t213 + pkin(4) * t216 + t180 * t232;
t210 = t133 * rSges(6,1) + t132 * rSges(6,2) + rSges(6,3) * t211 + pkin(4) * t215 - t183 * t232;
t209 = (rSges(6,1) * t181 - rSges(6,2) * t178 + t224) * t182 + (qJ(5) + rSges(6,3)) * t179;
t208 = qJD(2) * t183;
t207 = qJD(4) * t182;
t206 = qJD(5) * t182;
t161 = pkin(1) * t183 + t180 * qJ(2);
t205 = t170 * t161 + V_base(2);
t157 = t180 * pkin(1) - qJ(2) * t183;
t204 = V_base(4) * t157 + V_base(3);
t203 = V_base(5) * pkin(5) + V_base(1);
t200 = -t157 - t226;
t199 = qJD(2) * t180 + t203;
t198 = V_base(5) * pkin(2) + t199;
t197 = pkin(3) * t179 - pkin(7) * t182;
t196 = rSges(4,1) * t179 + rSges(4,2) * t182;
t193 = Icges(4,5) * t179 + Icges(4,6) * t182;
t189 = (Icges(4,3) * t183 + t180 * t193) * t167 + (Icges(4,3) * t180 - t183 * t193) * t166 + (Icges(4,5) * t182 - Icges(4,6) * t179) * t170;
t188 = t170 * t225 + (-pkin(2) - pkin(5)) * V_base(4) + t205;
t187 = V_base(4) * t226 + (-t161 - t225) * V_base(5) + t204;
t135 = t197 * t180;
t164 = pkin(3) * t182 + pkin(7) * t179;
t186 = t170 * t135 - t167 * t164 + t188;
t136 = t197 * t183;
t185 = t166 * t164 + (t136 + t200) * t170 + t198;
t184 = -t166 * t135 - t167 * t136 + t187;
t163 = rSges(2,1) * t183 - t180 * rSges(2,2);
t162 = -rSges(3,2) * t183 + t180 * rSges(3,3);
t160 = rSges(4,1) * t182 - rSges(4,2) * t179;
t159 = t180 * rSges(2,1) + rSges(2,2) * t183;
t158 = -t180 * rSges(3,2) - rSges(3,3) * t183;
t156 = qJD(4) * t179 + t170;
t140 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t139 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t138 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = -t180 * t207 + t167;
t128 = t183 * t207 + t166;
t126 = t180 * rSges(4,3) - t183 * t196;
t125 = rSges(5,3) * t179 + (rSges(5,1) * t181 - rSges(5,2) * t178) * t182;
t123 = rSges(4,3) * t183 + t180 * t196;
t107 = V_base(5) * rSges(2,3) - t159 * t170 + t203;
t106 = t163 * t170 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t105 = t159 * V_base(4) - t163 * V_base(5) + V_base(3);
t104 = t133 * rSges(5,1) + t132 * rSges(5,2) + rSges(5,3) * t211;
t102 = rSges(5,1) * t131 + rSges(5,2) * t130 - rSges(5,3) * t213;
t86 = V_base(5) * rSges(3,1) + (-t157 - t158) * t170 + t199;
t85 = -t208 + t170 * t162 + (-rSges(3,1) - pkin(5)) * V_base(4) + t205;
t84 = t158 * V_base(4) + (-t161 - t162) * V_base(5) + t204;
t83 = t160 * t166 + (-t126 + t200) * t170 + t198;
t82 = t170 * t123 - t167 * t160 + t188 - t208;
t81 = -t166 * t123 + t167 * t126 + t187;
t80 = -t104 * t156 + t125 * t128 + t185;
t79 = t156 * t102 - t129 * t125 + t186 - t208;
t78 = -t128 * t102 + t129 * t104 + t184;
t77 = t128 * t209 - t156 * t210 - t180 * t206 + t185;
t76 = (-qJD(2) + t206) * t183 + t222 * t156 - t209 * t129 + t186;
t75 = qJD(5) * t179 - t128 * t222 + t129 * t210 + t184;
t1 = m(1) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t167 * (t231 * t180 + t189 * t183) / 0.2e1 + t166 * (t189 * t180 - t231 * t183) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + ((t132 * t240 + t133 * t239 + t211 * t241) * t156 + (t245 * t132 + t243 * t133 + t211 * t247) * t129 + (t244 * t132 + t242 * t133 + t246 * t211) * t128) * t128 / 0.2e1 + ((t130 * t240 + t131 * t239 - t213 * t241) * t156 + (t245 * t130 + t243 * t131 - t247 * t213) * t129 + (t130 * t244 + t131 * t242 - t213 * t246) * t128) * t129 / 0.2e1 + (((-t178 * t240 + t181 * t239) * t156 + (-t178 * t245 + t181 * t243) * t129 + (-t178 * t244 + t181 * t242) * t128) * t182 + (t246 * t128 + t129 * t247 + t241 * t156) * t179) * t156 / 0.2e1 + ((-t117 * t179 + t121 * t182) * t167 + (-t118 * t179 + t122 * t182) * t166 + (-t148 * t179 + t153 * t182 + Icges(3,1) + Icges(2,3)) * t170) * t170 / 0.2e1 + ((t180 * t238 + t183 * t236 + Icges(1,4)) * V_base(5) + (t237 * t180 + t235 * t183 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t236 * t180 - t238 * t183 + Icges(1,2)) * V_base(5) + (t180 * t235 - t183 * t237 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t170 * (t257 * t180 - t255 * t183) + V_base(4) * t170 * (t255 * t180 + t257 * t183) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
