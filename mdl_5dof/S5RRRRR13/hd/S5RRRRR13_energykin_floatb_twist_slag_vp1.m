% Calculate kinetic energy for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:18
% EndTime: 2024-09-27 17:30:19
% DurationCPUTime: 0.42s
% Computational Cost: add. (2037->283), mult. (1330->395), div. (0->0), fcn. (1132->16), ass. (0->139)
t246 = -pkin(6) - pkin(7);
t213 = sin(qJ(1));
t245 = pkin(1) * t213;
t215 = cos(qJ(1));
t244 = pkin(1) * t215;
t209 = qJ(1) + qJ(2);
t200 = sin(t209);
t243 = pkin(2) * t200;
t202 = cos(t209);
t242 = pkin(2) * t202;
t211 = cos(pkin(5));
t241 = t211 * pkin(9);
t214 = cos(qJ(4));
t240 = pkin(4) * t214;
t239 = Icges(2,4) * t213;
t238 = Icges(3,4) * t200;
t204 = qJ(3) + t209;
t192 = sin(t204);
t237 = Icges(4,4) * t192;
t210 = sin(pkin(5));
t236 = t192 * t210;
t193 = cos(t204);
t235 = t193 * t210;
t212 = sin(qJ(4));
t234 = t211 * t212;
t233 = t211 * t214;
t208 = qJ(4) + qJ(5);
t232 = qJD(4) * t210;
t231 = -pkin(8) + t246;
t198 = V_base(6) + qJD(1);
t230 = t198 * t244 + V_base(2);
t229 = V_base(4) * t245 + V_base(3);
t228 = V_base(5) * pkin(6) + V_base(1);
t156 = t192 * t232 + V_base(4);
t225 = pkin(4) * t234 - pkin(10) * t210;
t191 = qJD(2) + t198;
t224 = t191 * t242 + t230;
t223 = V_base(4) * t243 + t229;
t222 = -t242 - t244;
t183 = qJD(3) + t191;
t166 = qJD(4) * t211 + t183;
t221 = V_base(5) * pkin(7) - t198 * t245 + t228;
t220 = V_base(5) * pkin(8) - t191 * t243 + t221;
t141 = pkin(3) * t193 + pkin(9) * t236;
t219 = t183 * t141 + (t231 - t241) * V_base(4) + t224;
t140 = pkin(3) * t192 - pkin(9) * t235;
t218 = V_base(4) * t140 + (-t141 + t222) * V_base(5) + t223;
t217 = -t140 * t183 + V_base(5) * t241 + t220;
t203 = Icges(2,4) * t215;
t201 = cos(t208);
t199 = sin(t208);
t197 = pkin(5) - t208;
t196 = pkin(5) + t208;
t190 = cos(t196);
t189 = sin(t197);
t188 = Icges(3,4) * t202;
t185 = cos(t197) / 0.2e1;
t184 = sin(t196) / 0.2e1;
t182 = Icges(4,4) * t193;
t179 = rSges(2,1) * t215 - rSges(2,2) * t213;
t178 = rSges(2,1) * t213 + rSges(2,2) * t215;
t177 = Icges(2,1) * t215 - t239;
t176 = Icges(2,1) * t213 + t203;
t175 = -Icges(2,2) * t213 + t203;
t174 = Icges(2,2) * t215 + t239;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = rSges(3,1) * t202 - rSges(3,2) * t200;
t164 = rSges(3,1) * t200 + rSges(3,2) * t202;
t163 = Icges(3,1) * t202 - t238;
t162 = Icges(3,1) * t200 + t188;
t161 = -Icges(3,2) * t200 + t188;
t160 = Icges(3,2) * t202 + t238;
t155 = -t193 * t232 + V_base(5);
t154 = t185 - t190 / 0.2e1;
t153 = t185 + t190 / 0.2e1;
t152 = t184 - t189 / 0.2e1;
t151 = t184 + t189 / 0.2e1;
t150 = rSges(4,1) * t193 - rSges(4,2) * t192;
t149 = rSges(4,1) * t192 + rSges(4,2) * t193;
t148 = Icges(4,1) * t193 - t237;
t147 = Icges(4,1) * t192 + t182;
t146 = -Icges(4,2) * t192 + t182;
t145 = Icges(4,2) * t193 + t237;
t142 = qJD(5) * t211 + t166;
t139 = t210 * t212 * pkin(4) + pkin(10) * t211;
t138 = rSges(5,3) * t211 + (rSges(5,1) * t212 + rSges(5,2) * t214) * t210;
t137 = Icges(5,5) * t211 + (Icges(5,1) * t212 + Icges(5,4) * t214) * t210;
t136 = Icges(5,6) * t211 + (Icges(5,4) * t212 + Icges(5,2) * t214) * t210;
t135 = Icges(5,3) * t211 + (Icges(5,5) * t212 + Icges(5,6) * t214) * t210;
t133 = -t192 * t234 + t193 * t214;
t132 = -t192 * t233 - t193 * t212;
t131 = t192 * t214 + t193 * t234;
t130 = -t192 * t212 + t193 * t233;
t129 = qJD(5) * t236 + t156;
t128 = V_base(5) + (-qJD(4) - qJD(5)) * t235;
t127 = V_base(5) * rSges(2,3) - t178 * t198 + t228;
t126 = t179 * t198 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t124 = t178 * V_base(4) - t179 * V_base(5) + V_base(3);
t123 = -t152 * t192 + t193 * t201;
t122 = -t153 * t192 - t193 * t199;
t121 = t152 * t193 + t192 * t201;
t120 = t153 * t193 - t192 * t199;
t119 = rSges(6,1) * t154 + rSges(6,2) * t151 + rSges(6,3) * t211;
t118 = Icges(6,1) * t154 + Icges(6,4) * t151 + Icges(6,5) * t211;
t117 = Icges(6,4) * t154 + Icges(6,2) * t151 + Icges(6,6) * t211;
t116 = Icges(6,5) * t154 + Icges(6,6) * t151 + Icges(6,3) * t211;
t115 = V_base(5) * rSges(3,3) - t164 * t191 + t221;
t114 = t165 * t191 + (-rSges(3,3) + t246) * V_base(4) + t230;
t113 = t164 * V_base(4) + (-t165 - t244) * V_base(5) + t229;
t112 = -t192 * t225 + t193 * t240;
t111 = t192 * t240 + t193 * t225;
t110 = rSges(5,1) * t133 + rSges(5,2) * t132 + rSges(5,3) * t236;
t109 = rSges(5,1) * t131 + rSges(5,2) * t130 - rSges(5,3) * t235;
t108 = Icges(5,1) * t133 + Icges(5,4) * t132 + Icges(5,5) * t236;
t107 = Icges(5,1) * t131 + Icges(5,4) * t130 - Icges(5,5) * t235;
t106 = Icges(5,4) * t133 + Icges(5,2) * t132 + Icges(5,6) * t236;
t105 = Icges(5,4) * t131 + Icges(5,2) * t130 - Icges(5,6) * t235;
t104 = Icges(5,5) * t133 + Icges(5,6) * t132 + Icges(5,3) * t236;
t103 = Icges(5,5) * t131 + Icges(5,6) * t130 - Icges(5,3) * t235;
t102 = V_base(5) * rSges(4,3) - t149 * t183 + t220;
t101 = t150 * t183 + (-rSges(4,3) + t231) * V_base(4) + t224;
t100 = t149 * V_base(4) + (-t150 + t222) * V_base(5) + t223;
t99 = rSges(6,1) * t123 + rSges(6,2) * t122 + rSges(6,3) * t236;
t98 = rSges(6,1) * t121 + rSges(6,2) * t120 - rSges(6,3) * t235;
t97 = Icges(6,1) * t123 + Icges(6,4) * t122 + Icges(6,5) * t236;
t96 = Icges(6,1) * t121 + Icges(6,4) * t120 - Icges(6,5) * t235;
t95 = Icges(6,4) * t123 + Icges(6,2) * t122 + Icges(6,6) * t236;
t94 = Icges(6,4) * t121 + Icges(6,2) * t120 - Icges(6,6) * t235;
t93 = Icges(6,5) * t123 + Icges(6,6) * t122 + Icges(6,3) * t236;
t92 = Icges(6,5) * t121 + Icges(6,6) * t120 - Icges(6,3) * t235;
t91 = -t109 * t166 + t138 * t155 + t217;
t90 = t110 * t166 - t138 * t156 + t219;
t89 = t109 * t156 - t110 * t155 + t218;
t88 = -t111 * t166 + t119 * t128 + t139 * t155 - t142 * t98 + t217;
t87 = t112 * t166 - t119 * t129 - t139 * t156 + t142 * t99 + t219;
t86 = t111 * t156 - t112 * t155 - t128 * t99 + t129 * t98 + t218;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t124 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(3) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t156 * ((t104 * t236 + t106 * t132 + t108 * t133) * t156 + (t103 * t236 + t105 * t132 + t107 * t133) * t155 + (t132 * t136 + t133 * t137 + t135 * t236) * t166) / 0.2e1 + t155 * ((-t104 * t235 + t106 * t130 + t108 * t131) * t156 + (-t103 * t235 + t105 * t130 + t107 * t131) * t155 + (t130 * t136 + t131 * t137 - t135 * t235) * t166) / 0.2e1 + t166 * ((t103 * t155 + t104 * t156 + t135 * t166) * t211 + ((t106 * t214 + t108 * t212) * t156 + (t105 * t214 + t107 * t212) * t155 + (t136 * t214 + t137 * t212) * t166) * t210) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t129 * ((t122 * t95 + t123 * t97 + t236 * t93) * t129 + (t122 * t94 + t123 * t96 + t236 * t92) * t128 + (t116 * t236 + t117 * t122 + t118 * t123) * t142) / 0.2e1 + t128 * ((t120 * t95 + t121 * t97 - t235 * t93) * t129 + (t120 * t94 + t121 * t96 - t235 * t92) * t128 + (-t116 * t235 + t117 * t120 + t118 * t121) * t142) / 0.2e1 + t142 * ((t151 * t95 + t154 * t97 + t211 * t93) * t129 + (t151 * t94 + t154 * t96 + t211 * t92) * t128 + (t116 * t211 + t117 * t151 + t118 * t154) * t142) / 0.2e1 + ((-t145 * t192 + t147 * t193 - t160 * t200 + t162 * t202 - t174 * t213 + t176 * t215 + Icges(1,4)) * V_base(5) + (-t146 * t192 + t148 * t193 - t161 * t200 + t163 * t202 - t175 * t213 + t177 * t215 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t145 * t193 + t147 * t192 + t160 * t202 + t162 * t200 + t174 * t215 + t176 * t213 + Icges(1,2)) * V_base(5) + (t146 * t193 + t148 * t192 + t161 * t202 + t163 * t200 + t175 * t215 + t177 * t213 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t213 + Icges(2,6) * t215) * V_base(5) + (Icges(2,5) * t215 - Icges(2,6) * t213) * V_base(4) + Icges(2,3) * t198 / 0.2e1) * t198 + ((Icges(3,5) * t200 + Icges(3,6) * t202) * V_base(5) + (Icges(3,5) * t202 - Icges(3,6) * t200) * V_base(4) + Icges(3,3) * t191 / 0.2e1) * t191 + ((Icges(4,5) * t192 + Icges(4,6) * t193) * V_base(5) + (Icges(4,5) * t193 - Icges(4,6) * t192) * V_base(4) + Icges(4,3) * t183 / 0.2e1) * t183;
T = t1;
