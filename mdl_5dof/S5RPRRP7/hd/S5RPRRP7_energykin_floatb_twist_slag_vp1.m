% Calculate kinetic energy for
% S5RPRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:30
% EndTime: 2019-12-31 18:44:32
% DurationCPUTime: 2.05s
% Computational Cost: add. (1277->238), mult. (1404->330), div. (0->0), fcn. (1306->8), ass. (0->118)
t244 = Icges(5,1) + Icges(6,1);
t243 = -Icges(5,4) + Icges(6,5);
t242 = Icges(6,4) + Icges(5,5);
t241 = Icges(5,2) + Icges(6,3);
t240 = -Icges(6,6) + Icges(5,6);
t239 = -Icges(5,3) - Icges(6,2);
t238 = rSges(6,1) + pkin(4);
t237 = rSges(6,3) + qJ(5);
t181 = qJ(1) + pkin(8);
t175 = sin(t181);
t176 = cos(t181);
t185 = cos(qJ(4));
t182 = sin(qJ(4));
t186 = cos(qJ(3));
t211 = t182 * t186;
t124 = t175 * t211 + t176 * t185;
t210 = t185 * t186;
t125 = t175 * t210 - t176 * t182;
t183 = sin(qJ(3));
t213 = t175 * t183;
t236 = t241 * t124 + t243 * t125 - t240 * t213;
t126 = -t175 * t185 + t176 * t211;
t127 = t175 * t182 + t176 * t210;
t212 = t176 * t183;
t235 = t241 * t126 + t243 * t127 - t240 * t212;
t234 = -t240 * t124 + t242 * t125 - t239 * t213;
t233 = -t240 * t126 + t242 * t127 - t239 * t212;
t232 = t243 * t124 + t244 * t125 + t242 * t213;
t231 = t243 * t126 + t244 * t127 + t242 * t212;
t230 = t240 * t186 + (t241 * t182 + t243 * t185) * t183;
t229 = t239 * t186 + (-t240 * t182 + t242 * t185) * t183;
t228 = -t242 * t186 + (t243 * t182 + t244 * t185) * t183;
t184 = sin(qJ(1));
t221 = pkin(1) * t184;
t187 = cos(qJ(1));
t220 = pkin(1) * t187;
t219 = -pkin(5) - qJ(2);
t218 = rSges(6,2) * t213 + t237 * t124 + t238 * t125;
t217 = Icges(2,4) * t184;
t216 = Icges(3,4) * t175;
t215 = Icges(4,4) * t183;
t214 = Icges(4,4) * t186;
t209 = rSges(6,2) * t212 + t237 * t126 + t238 * t127;
t208 = -rSges(6,2) * t186 + (t237 * t182 + t238 * t185) * t183;
t207 = qJD(4) * t183;
t177 = V_base(6) + qJD(1);
t206 = t177 * t220 + V_base(2);
t205 = V_base(5) * pkin(5) + V_base(1);
t155 = qJD(3) * t175 + V_base(4);
t149 = pkin(2) * t175 - pkin(6) * t176;
t202 = -t149 - t221;
t201 = V_base(5) * qJ(2) + t205;
t200 = V_base(4) * t221 + qJD(2) + V_base(3);
t199 = pkin(3) * t186 + pkin(7) * t183;
t154 = -qJD(3) * t176 + V_base(5);
t198 = rSges(4,1) * t186 - rSges(4,2) * t183;
t197 = Icges(4,1) * t186 - t215;
t196 = -Icges(4,2) * t183 + t214;
t195 = Icges(4,5) * t186 - Icges(4,6) * t183;
t194 = (-Icges(4,3) * t176 + t175 * t195) * t154 + (Icges(4,3) * t175 + t176 * t195) * t155 + (Icges(4,5) * t183 + Icges(4,6) * t186) * t177;
t150 = pkin(2) * t176 + pkin(6) * t175;
t193 = t177 * t150 + t219 * V_base(4) + t206;
t137 = t199 * t175;
t169 = pkin(3) * t183 - pkin(7) * t186;
t192 = t154 * t169 + (-t137 + t202) * t177 + t201;
t191 = V_base(4) * t149 + (-t150 - t220) * V_base(5) + t200;
t138 = t199 * t176;
t190 = t177 * t138 - t155 * t169 + t193;
t189 = t155 * t137 - t154 * t138 + t191;
t112 = -Icges(4,6) * t176 + t175 * t196;
t113 = Icges(4,6) * t175 + t176 * t196;
t114 = -Icges(4,5) * t176 + t175 * t197;
t115 = Icges(4,5) * t175 + t176 * t197;
t159 = Icges(4,2) * t186 + t215;
t162 = Icges(4,1) * t183 + t214;
t188 = (-t113 * t183 + t115 * t186) * t155 + (-t112 * t183 + t114 * t186) * t154 + (-t159 * t183 + t162 * t186) * t177;
t179 = Icges(2,4) * t187;
t174 = Icges(3,4) * t176;
t168 = rSges(2,1) * t187 - t184 * rSges(2,2);
t167 = t184 * rSges(2,1) + rSges(2,2) * t187;
t166 = rSges(4,1) * t183 + rSges(4,2) * t186;
t165 = -qJD(4) * t186 + t177;
t164 = Icges(2,1) * t187 - t217;
t163 = Icges(2,1) * t184 + t179;
t161 = -Icges(2,2) * t184 + t179;
t160 = Icges(2,2) * t187 + t217;
t153 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t152 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t151 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t148 = rSges(3,1) * t176 - rSges(3,2) * t175;
t147 = rSges(3,1) * t175 + rSges(3,2) * t176;
t146 = Icges(3,1) * t176 - t216;
t145 = Icges(3,1) * t175 + t174;
t144 = -Icges(3,2) * t175 + t174;
t143 = Icges(3,2) * t176 + t216;
t135 = -rSges(5,3) * t186 + (rSges(5,1) * t185 - rSges(5,2) * t182) * t183;
t123 = t176 * t207 + t155;
t122 = t175 * t207 + t154;
t119 = rSges(4,3) * t175 + t176 * t198;
t118 = -rSges(4,3) * t176 + t175 * t198;
t117 = V_base(5) * rSges(2,3) - t167 * t177 + t205;
t116 = t168 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t109 = t167 * V_base(4) - t168 * V_base(5) + V_base(3);
t105 = V_base(5) * rSges(3,3) + (-t147 - t221) * t177 + t201;
t104 = t148 * t177 + (-rSges(3,3) + t219) * V_base(4) + t206;
t103 = V_base(4) * t147 + (-t148 - t220) * V_base(5) + t200;
t102 = rSges(5,1) * t127 - rSges(5,2) * t126 + rSges(5,3) * t212;
t100 = rSges(5,1) * t125 - rSges(5,2) * t124 + rSges(5,3) * t213;
t86 = t154 * t166 + (-t118 + t202) * t177 + t201;
t85 = t119 * t177 - t155 * t166 + t193;
t84 = t155 * t118 - t154 * t119 + t191;
t83 = -t100 * t165 + t122 * t135 + t192;
t82 = t102 * t165 - t123 * t135 + t190;
t81 = t123 * t100 - t122 * t102 + t189;
t80 = qJD(5) * t126 + t122 * t208 - t165 * t218 + t192;
t79 = qJD(5) * t124 - t123 * t208 + t165 * t209 + t190;
t78 = qJD(5) * t183 * t182 - t122 * t209 + t123 * t218 + t189;
t1 = m(1) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(2) * (t109 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t155 * (t194 * t175 + t188 * t176) / 0.2e1 + t154 * (t188 * t175 - t194 * t176) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + ((t230 * t124 + t228 * t125 + t229 * t213) * t165 + (t235 * t124 + t231 * t125 + t233 * t213) * t123 + (t236 * t124 + t232 * t125 + t234 * t213) * t122) * t122 / 0.2e1 + ((t230 * t126 + t228 * t127 + t229 * t212) * t165 + (t235 * t126 + t231 * t127 + t233 * t212) * t123 + (t236 * t126 + t232 * t127 + t234 * t212) * t122) * t123 / 0.2e1 + ((-t234 * t122 - t233 * t123 - t229 * t165) * t186 + ((t230 * t182 + t228 * t185) * t165 + (t235 * t182 + t231 * t185) * t123 + (t236 * t182 + t232 * t185) * t122) * t183) * t165 / 0.2e1 + ((t113 * t186 + t115 * t183) * t155 + (t112 * t186 + t114 * t183) * t154 + (t186 * t159 + t183 * t162 + Icges(2,3) + Icges(3,3)) * t177) * t177 / 0.2e1 + ((-t143 * t175 + t145 * t176 - t184 * t160 + t163 * t187 + Icges(1,4)) * V_base(5) + (-t175 * t144 + t176 * t146 - t184 * t161 + t187 * t164 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t176 * t143 + t175 * t145 + t187 * t160 + t184 * t163 + Icges(1,2)) * V_base(5) + (t144 * t176 + t146 * t175 + t161 * t187 + t184 * t164 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t177 * (Icges(2,5) * t184 + Icges(3,5) * t175 + Icges(2,6) * t187 + Icges(3,6) * t176) + V_base(4) * t177 * (Icges(2,5) * t187 + Icges(3,5) * t176 - Icges(2,6) * t184 - Icges(3,6) * t175) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
