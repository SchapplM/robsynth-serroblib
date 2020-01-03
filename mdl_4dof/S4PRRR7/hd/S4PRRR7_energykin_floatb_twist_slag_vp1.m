% Calculate kinetic energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:02
% DurationCPUTime: 1.89s
% Computational Cost: add. (1239->278), mult. (2890->417), div. (0->0), fcn. (3384->10), ass. (0->120)
t228 = cos(qJ(3));
t201 = cos(pkin(4));
t227 = pkin(5) * t201;
t198 = sin(pkin(8));
t226 = Icges(2,4) * t198;
t199 = sin(pkin(4));
t225 = t198 * t199;
t200 = cos(pkin(8));
t224 = t199 * t200;
t203 = sin(qJ(3));
t223 = t199 * t203;
t206 = cos(qJ(2));
t222 = t199 * t206;
t204 = sin(qJ(2));
t221 = t201 * t204;
t220 = t201 * t206;
t219 = qJD(2) * t199;
t218 = V_base(5) * qJ(1) + V_base(1);
t214 = qJD(1) + V_base(3);
t213 = t199 * t228;
t180 = t198 * t219 + V_base(4);
t191 = qJD(2) * t201 + V_base(6);
t166 = t198 * t220 + t200 * t204;
t148 = qJD(3) * t166 + t180;
t179 = -t200 * t219 + V_base(5);
t164 = t198 * t204 - t200 * t220;
t147 = qJD(3) * t164 + t179;
t168 = -qJD(3) * t222 + t191;
t174 = t198 * pkin(1) - pkin(5) * t224;
t212 = -V_base(6) * t174 + V_base(5) * t227 + t218;
t175 = t200 * pkin(1) + pkin(5) * t225;
t211 = V_base(4) * t174 - V_base(5) * t175 + t214;
t210 = V_base(6) * t175 + V_base(2) + (-qJ(1) - t227) * V_base(4);
t165 = t198 * t206 + t200 * t221;
t142 = t165 * pkin(2) + t164 * pkin(6);
t173 = (pkin(2) * t204 - pkin(6) * t206) * t199;
t209 = -t191 * t142 + t179 * t173 + t212;
t167 = -t198 * t221 + t200 * t206;
t143 = t167 * pkin(2) + t166 * pkin(6);
t208 = t180 * t142 - t179 * t143 + t211;
t207 = t191 * t143 - t180 * t173 + t210;
t205 = cos(qJ(4));
t202 = sin(qJ(4));
t196 = Icges(2,4) * t200;
t188 = t200 * rSges(2,1) - t198 * rSges(2,2);
t187 = t198 * rSges(2,1) + t200 * rSges(2,2);
t186 = Icges(2,1) * t200 - t226;
t185 = Icges(2,1) * t198 + t196;
t184 = -Icges(2,2) * t198 + t196;
t183 = Icges(2,2) * t200 + t226;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = t201 * t203 + t204 * t213;
t171 = -t201 * t228 + t204 * t223;
t161 = t201 * rSges(3,3) + (rSges(3,1) * t204 + rSges(3,2) * t206) * t199;
t160 = Icges(3,5) * t201 + (Icges(3,1) * t204 + Icges(3,4) * t206) * t199;
t159 = Icges(3,6) * t201 + (Icges(3,4) * t204 + Icges(3,2) * t206) * t199;
t158 = Icges(3,3) * t201 + (Icges(3,5) * t204 + Icges(3,6) * t206) * t199;
t157 = V_base(5) * rSges(2,3) - V_base(6) * t187 + t218;
t156 = V_base(6) * t188 + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t154 = t172 * t205 - t202 * t222;
t153 = -t172 * t202 - t205 * t222;
t152 = t167 * t228 + t198 * t223;
t151 = t167 * t203 - t198 * t213;
t150 = t165 * t228 - t200 * t223;
t149 = t165 * t203 + t200 * t213;
t146 = V_base(4) * t187 - V_base(5) * t188 + t214;
t145 = qJD(4) * t171 + t168;
t144 = t172 * pkin(3) + t171 * pkin(7);
t141 = t172 * rSges(4,1) - t171 * rSges(4,2) - rSges(4,3) * t222;
t140 = Icges(4,1) * t172 - Icges(4,4) * t171 - Icges(4,5) * t222;
t139 = Icges(4,4) * t172 - Icges(4,2) * t171 - Icges(4,6) * t222;
t138 = Icges(4,5) * t172 - Icges(4,6) * t171 - Icges(4,3) * t222;
t137 = t167 * rSges(3,1) - t166 * rSges(3,2) + rSges(3,3) * t225;
t136 = t165 * rSges(3,1) - t164 * rSges(3,2) - rSges(3,3) * t224;
t135 = Icges(3,1) * t167 - Icges(3,4) * t166 + Icges(3,5) * t225;
t134 = Icges(3,1) * t165 - Icges(3,4) * t164 - Icges(3,5) * t224;
t133 = Icges(3,4) * t167 - Icges(3,2) * t166 + Icges(3,6) * t225;
t132 = Icges(3,4) * t165 - Icges(3,2) * t164 - Icges(3,6) * t224;
t131 = Icges(3,5) * t167 - Icges(3,6) * t166 + Icges(3,3) * t225;
t130 = Icges(3,5) * t165 - Icges(3,6) * t164 - Icges(3,3) * t224;
t127 = t152 * t205 + t166 * t202;
t126 = -t152 * t202 + t166 * t205;
t125 = t150 * t205 + t164 * t202;
t124 = -t150 * t202 + t164 * t205;
t123 = qJD(4) * t151 + t148;
t122 = qJD(4) * t149 + t147;
t121 = t152 * pkin(3) + t151 * pkin(7);
t120 = t150 * pkin(3) + t149 * pkin(7);
t119 = t154 * rSges(5,1) + t153 * rSges(5,2) + t171 * rSges(5,3);
t118 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t171;
t117 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t171;
t116 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t171;
t115 = t152 * rSges(4,1) - t151 * rSges(4,2) + t166 * rSges(4,3);
t114 = t150 * rSges(4,1) - t149 * rSges(4,2) + t164 * rSges(4,3);
t113 = Icges(4,1) * t152 - Icges(4,4) * t151 + Icges(4,5) * t166;
t112 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t164;
t111 = Icges(4,4) * t152 - Icges(4,2) * t151 + Icges(4,6) * t166;
t110 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t164;
t109 = Icges(4,5) * t152 - Icges(4,6) * t151 + Icges(4,3) * t166;
t108 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t164;
t107 = -t191 * t136 + t179 * t161 + t212;
t106 = t191 * t137 - t180 * t161 + t210;
t105 = t127 * rSges(5,1) + t126 * rSges(5,2) + t151 * rSges(5,3);
t104 = t125 * rSges(5,1) + t124 * rSges(5,2) + t149 * rSges(5,3);
t103 = Icges(5,1) * t127 + Icges(5,4) * t126 + Icges(5,5) * t151;
t102 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t149;
t101 = Icges(5,4) * t127 + Icges(5,2) * t126 + Icges(5,6) * t151;
t100 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t149;
t99 = Icges(5,5) * t127 + Icges(5,6) * t126 + Icges(5,3) * t151;
t98 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t149;
t97 = t180 * t136 - t179 * t137 + t211;
t96 = -t168 * t114 + t147 * t141 + t209;
t95 = t168 * t115 - t148 * t141 + t207;
t94 = t148 * t114 - t147 * t115 + t208;
t93 = -t145 * t104 + t122 * t119 - t168 * t120 + t147 * t144 + t209;
t92 = t145 * t105 - t123 * t119 + t168 * t121 - t148 * t144 + t207;
t91 = t123 * t104 - t122 * t105 + t148 * t120 - t147 * t121 + t208;
t1 = m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(2) * (t146 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + t180 * ((t131 * t225 - t166 * t133 + t167 * t135) * t180 + (t130 * t225 - t166 * t132 + t167 * t134) * t179 + (t158 * t225 - t166 * t159 + t167 * t160) * t191) / 0.2e1 + t179 * ((-t131 * t224 - t164 * t133 + t165 * t135) * t180 + (-t130 * t224 - t164 * t132 + t165 * t134) * t179 + (-t158 * t224 - t164 * t159 + t165 * t160) * t191) / 0.2e1 + t191 * ((t130 * t179 + t131 * t180 + t158 * t191) * t201 + ((t133 * t206 + t135 * t204) * t180 + (t132 * t206 + t134 * t204) * t179 + (t159 * t206 + t160 * t204) * t191) * t199) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t148 * ((t166 * t109 - t151 * t111 + t152 * t113) * t148 + (t166 * t108 - t151 * t110 + t152 * t112) * t147 + (t166 * t138 - t151 * t139 + t152 * t140) * t168) / 0.2e1 + t147 * ((t164 * t109 - t149 * t111 + t150 * t113) * t148 + (t164 * t108 - t149 * t110 + t150 * t112) * t147 + (t164 * t138 - t149 * t139 + t150 * t140) * t168) / 0.2e1 + t168 * ((-t109 * t222 - t171 * t111 + t172 * t113) * t148 + (-t108 * t222 - t171 * t110 + t172 * t112) * t147 + (-t138 * t222 - t171 * t139 + t172 * t140) * t168) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t123 * ((t126 * t101 + t127 * t103 + t151 * t99) * t123 + (t126 * t100 + t127 * t102 + t151 * t98) * t122 + (t151 * t116 + t126 * t117 + t127 * t118) * t145) / 0.2e1 + t122 * ((t124 * t101 + t125 * t103 + t149 * t99) * t123 + (t124 * t100 + t125 * t102 + t149 * t98) * t122 + (t149 * t116 + t124 * t117 + t125 * t118) * t145) / 0.2e1 + t145 * ((t153 * t101 + t154 * t103 + t171 * t99) * t123 + (t153 * t100 + t154 * t102 + t171 * t98) * t122 + (t171 * t116 + t153 * t117 + t154 * t118) * t145) / 0.2e1 + ((-t198 * t183 + t200 * t185 + Icges(1,4)) * V_base(5) + (-t198 * t184 + t200 * t186 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t200 * t183 + t198 * t185 + Icges(1,2)) * V_base(5) + (t200 * t184 + t198 * t186 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t198 + Icges(2,6) * t200 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t200 - Icges(2,6) * t198 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
