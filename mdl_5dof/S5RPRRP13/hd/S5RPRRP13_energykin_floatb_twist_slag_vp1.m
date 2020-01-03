% Calculate kinetic energy for
% S5RPRRP13
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:28
% EndTime: 2019-12-31 18:58:30
% DurationCPUTime: 2.28s
% Computational Cost: add. (749->227), mult. (1414->312), div. (0->0), fcn. (1318->6), ass. (0->116)
t258 = Icges(2,4) + Icges(3,6);
t257 = Icges(2,1) + Icges(3,2);
t256 = Icges(5,1) + Icges(6,1);
t255 = -Icges(3,4) + Icges(2,5);
t254 = Icges(5,4) - Icges(6,5);
t253 = Icges(6,4) + Icges(5,5);
t252 = Icges(3,5) - Icges(2,6);
t251 = Icges(2,2) + Icges(3,3);
t250 = Icges(5,2) + Icges(6,3);
t249 = Icges(6,6) - Icges(5,6);
t248 = Icges(5,3) + Icges(6,2);
t247 = rSges(6,1) + pkin(4);
t246 = rSges(6,3) + qJ(5);
t184 = cos(qJ(1));
t245 = t258 * t184;
t181 = sin(qJ(1));
t244 = t258 * t181;
t180 = sin(qJ(3));
t182 = cos(qJ(4));
t211 = t184 * t182;
t179 = sin(qJ(4));
t215 = t181 * t179;
t131 = t180 * t215 - t211;
t214 = t181 * t182;
t216 = t179 * t184;
t132 = t180 * t214 + t216;
t183 = cos(qJ(3));
t213 = t181 * t183;
t243 = t250 * t131 - t254 * t132 - t249 * t213;
t133 = t180 * t216 + t214;
t134 = -t180 * t211 + t215;
t212 = t183 * t184;
t242 = -t250 * t133 - t254 * t134 + t249 * t212;
t241 = t249 * t131 + t253 * t132 - t248 * t213;
t240 = -t249 * t133 + t253 * t134 + t248 * t212;
t239 = -t254 * t131 + t256 * t132 - t253 * t213;
t238 = t254 * t133 + t256 * t134 + t253 * t212;
t237 = (t250 * t179 - t254 * t182) * t183 + t249 * t180;
t236 = (t249 * t179 + t253 * t182) * t183 + t248 * t180;
t235 = (-t254 * t179 + t256 * t182) * t183 + t253 * t180;
t234 = -t251 * t184 - t244;
t233 = t251 * t181 - t245;
t232 = t257 * t181 + t245;
t231 = t257 * t184 - t244;
t220 = Icges(4,4) * t180;
t196 = Icges(4,2) * t183 + t220;
t118 = Icges(4,6) * t184 + t181 * t196;
t119 = Icges(4,6) * t181 - t184 * t196;
t219 = Icges(4,4) * t183;
t197 = Icges(4,1) * t180 + t219;
t122 = Icges(4,5) * t184 + t181 * t197;
t123 = Icges(4,5) * t181 - t184 * t197;
t150 = -Icges(4,2) * t180 + t219;
t155 = Icges(4,1) * t183 - t220;
t168 = qJD(3) * t181 + V_base(5);
t169 = qJD(3) * t184 + V_base(4);
t172 = V_base(6) + qJD(1);
t228 = (t118 * t183 + t122 * t180) * t169 + (t119 * t183 + t123 * t180) * t168 + (t150 * t183 + t155 * t180) * t172;
t223 = pkin(6) * t181;
t222 = pkin(6) * t184;
t210 = -rSges(6,2) * t213 + t246 * t131 + t247 * t132;
t209 = rSges(6,2) * t212 - t246 * t133 + t247 * t134;
t208 = rSges(6,2) * t180 + (t246 * t179 + t247 * t182) * t183;
t207 = qJD(4) * t183;
t159 = t181 * pkin(1) - qJ(2) * t184;
t206 = V_base(4) * t159 + V_base(3);
t205 = V_base(5) * pkin(5) + V_base(1);
t202 = -t159 - t223;
t201 = qJD(2) * t181 + t205;
t200 = V_base(5) * pkin(2) + t201;
t199 = pkin(3) * t180 - pkin(7) * t183;
t198 = rSges(4,1) * t180 + rSges(4,2) * t183;
t195 = Icges(4,5) * t180 + Icges(4,6) * t183;
t163 = pkin(1) * t184 + t181 * qJ(2);
t191 = -qJD(2) * t184 + t172 * t163 + V_base(2);
t190 = (Icges(4,3) * t184 + t181 * t195) * t169 + (Icges(4,3) * t181 - t184 * t195) * t168 + (Icges(4,5) * t183 - Icges(4,6) * t180) * t172;
t189 = V_base(4) * t223 + (-t163 - t222) * V_base(5) + t206;
t188 = t172 * t222 + (-pkin(2) - pkin(5)) * V_base(4) + t191;
t138 = t199 * t184;
t166 = pkin(3) * t183 + pkin(7) * t180;
t187 = t168 * t166 + (t138 + t202) * t172 + t200;
t137 = t199 * t181;
t186 = -t168 * t137 - t169 * t138 + t189;
t185 = t172 * t137 - t169 * t166 + t188;
t165 = rSges(2,1) * t184 - t181 * rSges(2,2);
t164 = -rSges(3,2) * t184 + t181 * rSges(3,3);
t162 = rSges(4,1) * t183 - rSges(4,2) * t180;
t161 = t181 * rSges(2,1) + rSges(2,2) * t184;
t160 = -t181 * rSges(3,2) - rSges(3,3) * t184;
t158 = qJD(4) * t180 + t172;
t142 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t141 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t140 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = -t181 * t207 + t169;
t129 = t184 * t207 + t168;
t127 = t181 * rSges(4,3) - t184 * t198;
t126 = rSges(5,3) * t180 + (rSges(5,1) * t182 - rSges(5,2) * t179) * t183;
t124 = rSges(4,3) * t184 + t181 * t198;
t109 = V_base(5) * rSges(2,3) - t161 * t172 + t205;
t108 = t165 * t172 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t107 = t161 * V_base(4) - t165 * V_base(5) + V_base(3);
t104 = t134 * rSges(5,1) + t133 * rSges(5,2) + rSges(5,3) * t212;
t102 = rSges(5,1) * t132 - rSges(5,2) * t131 - rSges(5,3) * t213;
t88 = V_base(5) * rSges(3,1) + (-t159 - t160) * t172 + t201;
t87 = t172 * t164 + (-rSges(3,1) - pkin(5)) * V_base(4) + t191;
t86 = t160 * V_base(4) + (-t163 - t164) * V_base(5) + t206;
t85 = t162 * t168 + (-t127 + t202) * t172 + t200;
t84 = t172 * t124 - t169 * t162 + t188;
t83 = -t168 * t124 + t169 * t127 + t189;
t82 = -t104 * t158 + t126 * t129 + t187;
t81 = t158 * t102 - t130 * t126 + t185;
t80 = -t129 * t102 + t130 * t104 + t186;
t79 = qJD(5) * t131 + t129 * t208 - t158 * t209 + t187;
t78 = -qJD(5) * t133 - t130 * t208 + t158 * t210 + t185;
t77 = qJD(5) * t183 * t179 - t129 * t210 + t130 * t209 + t186;
t1 = m(1) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t169 * (t228 * t181 + t190 * t184) / 0.2e1 + t168 * (t190 * t181 - t228 * t184) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + ((-t133 * t237 + t134 * t235 + t212 * t236) * t158 + (-t133 * t243 + t239 * t134 + t241 * t212) * t130 + (-t242 * t133 + t238 * t134 + t240 * t212) * t129) * t129 / 0.2e1 + ((t131 * t237 + t132 * t235 - t213 * t236) * t158 + (t243 * t131 + t239 * t132 - t241 * t213) * t130 + (t131 * t242 + t132 * t238 - t213 * t240) * t129) * t130 / 0.2e1 + (((t179 * t237 + t182 * t235) * t158 + (t179 * t243 + t239 * t182) * t130 + (t179 * t242 + t182 * t238) * t129) * t183 + (t129 * t240 + t130 * t241 + t158 * t236) * t180) * t158 / 0.2e1 + ((-t118 * t180 + t122 * t183) * t169 + (-t119 * t180 + t123 * t183) * t168 + (-t150 * t180 + t155 * t183 + Icges(3,1) + Icges(2,3)) * t172) * t172 / 0.2e1 + ((t181 * t234 + t184 * t232 + Icges(1,4)) * V_base(5) + (t233 * t181 + t231 * t184 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t232 * t181 - t234 * t184 + Icges(1,2)) * V_base(5) + (t181 * t231 - t184 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t172 * (t255 * t181 - t252 * t184) + V_base(4) * t172 * (t252 * t181 + t255 * t184) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
