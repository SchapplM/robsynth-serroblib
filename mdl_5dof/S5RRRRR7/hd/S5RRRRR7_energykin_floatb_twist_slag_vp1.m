% Calculate kinetic energy for
% S5RRRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:19
% EndTime: 2019-12-31 22:21:20
% DurationCPUTime: 1.86s
% Computational Cost: add. (1452->283), mult. (1432->429), div. (0->0), fcn. (1270->10), ass. (0->149)
t199 = sin(qJ(2));
t257 = pkin(2) * t199;
t197 = qJ(2) + qJ(3);
t190 = sin(t197);
t256 = pkin(3) * t190;
t202 = cos(qJ(2));
t255 = t202 * pkin(2);
t200 = sin(qJ(1));
t253 = Icges(2,4) * t200;
t252 = Icges(3,4) * t199;
t251 = Icges(3,4) * t202;
t250 = Icges(4,4) * t190;
t191 = cos(t197);
t249 = Icges(4,4) * t191;
t193 = qJ(4) + t197;
t183 = sin(t193);
t248 = Icges(5,4) * t183;
t184 = cos(t193);
t247 = Icges(5,4) * t184;
t246 = t183 * t200;
t203 = cos(qJ(1));
t245 = t183 * t203;
t198 = sin(qJ(5));
t244 = t198 * t200;
t243 = t198 * t203;
t201 = cos(qJ(5));
t242 = t200 * t201;
t241 = t201 * t203;
t121 = -pkin(7) * t203 + t255 * t200;
t179 = t200 * pkin(1) - t203 * pkin(6);
t240 = -t121 - t179;
t239 = pkin(3) * t191;
t237 = qJD(5) * t183;
t236 = -qJD(2) - qJD(3);
t235 = V_base(5) * pkin(5) + V_base(1);
t100 = -pkin(8) * t203 + t239 * t200;
t232 = -t100 + t240;
t182 = qJD(2) * t200 + V_base(4);
t186 = V_base(6) + qJD(1);
t181 = -qJD(2) * t203 + V_base(5);
t231 = t181 * t257 + t235;
t159 = qJD(3) * t200 + t182;
t158 = t236 * t203 + V_base(5);
t230 = t158 * t256 + t231;
t229 = pkin(4) * t184 + pkin(9) * t183;
t228 = rSges(3,1) * t202 - rSges(3,2) * t199;
t227 = rSges(4,1) * t191 - rSges(4,2) * t190;
t226 = rSges(5,1) * t184 - rSges(5,2) * t183;
t147 = qJD(4) * t200 + t159;
t225 = Icges(3,1) * t202 - t252;
t224 = Icges(4,1) * t191 - t250;
t223 = Icges(5,1) * t184 - t248;
t222 = -Icges(3,2) * t199 + t251;
t221 = -Icges(4,2) * t190 + t249;
t220 = -Icges(5,2) * t183 + t247;
t219 = Icges(3,5) * t202 - Icges(3,6) * t199;
t218 = Icges(4,5) * t191 - Icges(4,6) * t190;
t217 = Icges(5,5) * t184 - Icges(5,6) * t183;
t180 = t203 * pkin(1) + t200 * pkin(6);
t216 = -V_base(4) * pkin(5) + t186 * t180 + V_base(2);
t215 = V_base(4) * t179 - t180 * V_base(5) + V_base(3);
t146 = V_base(5) + (-qJD(4) + t236) * t203;
t214 = (-Icges(5,3) * t203 + t217 * t200) * t146 + (Icges(5,3) * t200 + t217 * t203) * t147 + (Icges(5,5) * t183 + Icges(5,6) * t184) * t186;
t213 = (-Icges(4,3) * t203 + t218 * t200) * t158 + (Icges(4,3) * t200 + t218 * t203) * t159 + (Icges(4,5) * t190 + Icges(4,6) * t191) * t186;
t212 = (-Icges(3,3) * t203 + t219 * t200) * t181 + (Icges(3,3) * t200 + t219 * t203) * t182 + (Icges(3,5) * t199 + Icges(3,6) * t202) * t186;
t122 = pkin(7) * t200 + t255 * t203;
t211 = t182 * t121 - t122 * t181 + t215;
t210 = t186 * t122 - t182 * t257 + t216;
t101 = pkin(8) * t200 + t239 * t203;
t209 = t159 * t100 - t101 * t158 + t211;
t208 = t186 * t101 - t159 * t256 + t210;
t115 = -Icges(5,6) * t203 + t220 * t200;
t116 = Icges(5,6) * t200 + t220 * t203;
t117 = -Icges(5,5) * t203 + t223 * t200;
t118 = Icges(5,5) * t200 + t223 * t203;
t149 = Icges(5,2) * t184 + t248;
t150 = Icges(5,1) * t183 + t247;
t207 = (-t116 * t183 + t118 * t184) * t147 + (-t115 * t183 + t117 * t184) * t146 + (-t149 * t183 + t150 * t184) * t186;
t125 = -Icges(4,6) * t203 + t221 * t200;
t126 = Icges(4,6) * t200 + t221 * t203;
t127 = -Icges(4,5) * t203 + t224 * t200;
t128 = Icges(4,5) * t200 + t224 * t203;
t155 = Icges(4,2) * t191 + t250;
t156 = Icges(4,1) * t190 + t249;
t206 = (-t126 * t190 + t128 * t191) * t159 + (-t125 * t190 + t127 * t191) * t158 + (-t155 * t190 + t156 * t191) * t186;
t135 = -Icges(3,6) * t203 + t222 * t200;
t136 = Icges(3,6) * t200 + t222 * t203;
t137 = -Icges(3,5) * t203 + t225 * t200;
t138 = Icges(3,5) * t200 + t225 * t203;
t170 = Icges(3,2) * t202 + t252;
t173 = Icges(3,1) * t199 + t251;
t205 = (-t136 * t199 + t138 * t202) * t182 + (-t135 * t199 + t137 * t202) * t181 + (-t170 * t199 + t173 * t202) * t186;
t192 = Icges(2,4) * t203;
t178 = rSges(2,1) * t203 - rSges(2,2) * t200;
t177 = rSges(2,1) * t200 + rSges(2,2) * t203;
t176 = rSges(3,1) * t199 + rSges(3,2) * t202;
t175 = Icges(2,1) * t203 - t253;
t174 = Icges(2,1) * t200 + t192;
t172 = -Icges(2,2) * t200 + t192;
t171 = Icges(2,2) * t203 + t253;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t160 = -qJD(5) * t184 + t186;
t157 = rSges(4,1) * t190 + rSges(4,2) * t191;
t152 = pkin(4) * t183 - pkin(9) * t184;
t151 = rSges(5,1) * t183 + rSges(5,2) * t184;
t145 = t184 * t241 + t244;
t144 = -t184 * t243 + t242;
t143 = t184 * t242 - t243;
t142 = -t184 * t244 - t241;
t140 = rSges(3,3) * t200 + t228 * t203;
t139 = -rSges(3,3) * t203 + t228 * t200;
t132 = t229 * t203;
t131 = t229 * t200;
t130 = rSges(4,3) * t200 + t227 * t203;
t129 = -rSges(4,3) * t203 + t227 * t200;
t120 = rSges(5,3) * t200 + t226 * t203;
t119 = -rSges(5,3) * t203 + t226 * t200;
t112 = V_base(5) * rSges(2,3) - t177 * t186 + t235;
t111 = t178 * t186 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t110 = t203 * t237 + t147;
t109 = t200 * t237 + t146;
t108 = t177 * V_base(4) - t178 * V_base(5) + V_base(3);
t107 = -rSges(6,3) * t184 + (rSges(6,1) * t201 - rSges(6,2) * t198) * t183;
t106 = -Icges(6,5) * t184 + (Icges(6,1) * t201 - Icges(6,4) * t198) * t183;
t105 = -Icges(6,6) * t184 + (Icges(6,4) * t201 - Icges(6,2) * t198) * t183;
t104 = -Icges(6,3) * t184 + (Icges(6,5) * t201 - Icges(6,6) * t198) * t183;
t97 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t245;
t96 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t246;
t95 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t245;
t94 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t246;
t93 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t245;
t92 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t246;
t91 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t245;
t90 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t246;
t89 = t176 * t181 + (-t139 - t179) * t186 + t235;
t88 = t140 * t186 - t176 * t182 + t216;
t87 = t139 * t182 - t140 * t181 + t215;
t86 = t157 * t158 + (-t129 + t240) * t186 + t231;
t85 = t130 * t186 - t157 * t159 + t210;
t84 = t129 * t159 - t130 * t158 + t211;
t83 = t146 * t151 + (-t119 + t232) * t186 + t230;
t82 = t120 * t186 - t147 * t151 + t208;
t81 = t119 * t147 - t120 * t146 + t209;
t80 = t107 * t109 + t146 * t152 - t160 * t96 + (-t131 + t232) * t186 + t230;
t79 = -t107 * t110 + t132 * t186 - t147 * t152 + t160 * t97 + t208;
t78 = -t109 * t97 + t110 * t96 + t131 * t147 - t132 * t146 + t209;
t1 = m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t182 * (t212 * t200 + t205 * t203) / 0.2e1 + t181 * (t205 * t200 - t212 * t203) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t159 * (t213 * t200 + t206 * t203) / 0.2e1 + t158 * (t206 * t200 - t213 * t203) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t147 * (t214 * t200 + t207 * t203) / 0.2e1 + t146 * (t207 * t200 - t214 * t203) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t110 * ((t144 * t93 + t145 * t95 + t91 * t245) * t110 + (t144 * t92 + t145 * t94 + t90 * t245) * t109 + (t104 * t245 + t105 * t144 + t106 * t145) * t160) / 0.2e1 + t109 * ((t142 * t93 + t143 * t95 + t91 * t246) * t110 + (t142 * t92 + t143 * t94 + t90 * t246) * t109 + (t104 * t246 + t105 * t142 + t106 * t143) * t160) / 0.2e1 + t160 * ((-t104 * t160 - t90 * t109 - t91 * t110) * t184 + ((-t198 * t93 + t201 * t95) * t110 + (-t198 * t92 + t201 * t94) * t109 + (-t105 * t198 + t106 * t201) * t160) * t183) / 0.2e1 + ((-t171 * t200 + t174 * t203 + Icges(1,4)) * V_base(5) + (-t172 * t200 + t175 * t203 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t171 * t203 + t174 * t200 + Icges(1,2)) * V_base(5) + (t172 * t203 + t175 * t200 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t202 + t138 * t199) * t182 + (t135 * t202 + t137 * t199) * t181 + (t126 * t191 + t128 * t190) * t159 + (t125 * t191 + t127 * t190) * t158 + (t116 * t184 + t118 * t183) * t147 + (t115 * t184 + t117 * t183) * t146 + (t149 * t184 + t150 * t183 + t155 * t191 + t156 * t190 + t170 * t202 + t173 * t199 + Icges(2,3)) * t186) * t186 / 0.2e1 + V_base(4) * t186 * (Icges(2,5) * t203 - Icges(2,6) * t200) + V_base(5) * t186 * (Icges(2,5) * t200 + Icges(2,6) * t203) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
