% Calculate kinetic energy for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:39
% EndTime: 2022-01-20 11:48:41
% DurationCPUTime: 2.05s
% Computational Cost: add. (1244->214), mult. (960->291), div. (0->0), fcn. (740->8), ass. (0->116)
t242 = Icges(5,4) + Icges(6,4);
t241 = Icges(5,1) + Icges(6,1);
t240 = Icges(5,2) + Icges(6,2);
t162 = qJ(3) + qJ(4);
t155 = cos(t162);
t239 = t242 * t155;
t153 = sin(t162);
t238 = t242 * t153;
t237 = Icges(5,5) + Icges(6,5);
t236 = Icges(5,6) + Icges(6,6);
t235 = -t240 * t153 + t239;
t234 = t241 * t155 - t238;
t233 = rSges(6,1) + pkin(4);
t163 = qJ(1) + qJ(2);
t154 = sin(t163);
t156 = cos(t163);
t232 = t235 * t154 - t236 * t156;
t231 = t236 * t154 + t235 * t156;
t230 = t234 * t154 - t237 * t156;
t229 = t237 * t154 + t234 * t156;
t228 = Icges(5,3) + Icges(6,3);
t227 = t240 * t155 + t238;
t226 = t241 * t153 + t239;
t225 = -t236 * t153 + t237 * t155;
t224 = rSges(6,3) + qJ(5);
t223 = -rSges(6,2) * t153 + t233 * t155;
t106 = V_base(5) + (-qJD(3) - qJD(4)) * t156;
t133 = qJD(3) * t154 + V_base(4);
t107 = qJD(4) * t154 + t133;
t152 = V_base(6) + qJD(1);
t150 = qJD(2) + t152;
t222 = (-t227 * t153 + t226 * t155) * t150 + (-t231 * t153 + t229 * t155) * t107 + (-t232 * t153 + t230 * t155) * t106;
t221 = (t237 * t153 + t236 * t155) * t150 + (t228 * t154 + t225 * t156) * t107 + (t225 * t154 - t228 * t156) * t106;
t220 = -pkin(5) - pkin(6);
t165 = sin(qJ(1));
t216 = pkin(1) * t165;
t167 = cos(qJ(1));
t215 = pkin(1) * t167;
t164 = sin(qJ(3));
t214 = pkin(3) * t164;
t166 = cos(qJ(3));
t213 = t166 * pkin(3);
t211 = t223 * t154 - t156 * t224;
t210 = t154 * t224 + t223 * t156;
t126 = t154 * pkin(2) - t156 * pkin(7);
t76 = -pkin(8) * t156 + t154 * t213;
t209 = -t126 - t76;
t208 = Icges(2,4) * t165;
t207 = Icges(3,4) * t154;
t206 = Icges(4,4) * t164;
t205 = Icges(4,4) * t166;
t198 = t152 * t215 + V_base(2);
t197 = V_base(4) * t216 + V_base(3);
t196 = V_base(5) * pkin(5) + V_base(1);
t193 = rSges(6,2) * t155 + t233 * t153;
t192 = rSges(4,1) * t166 - rSges(4,2) * t164;
t191 = rSges(5,1) * t155 - rSges(5,2) * t153;
t189 = Icges(4,1) * t166 - t206;
t186 = -Icges(4,2) * t164 + t205;
t183 = Icges(4,5) * t166 - Icges(4,6) * t164;
t180 = V_base(5) * pkin(6) - t152 * t216 + t196;
t132 = -qJD(3) * t156 + V_base(5);
t177 = (-Icges(4,3) * t156 + t154 * t183) * t132 + (Icges(4,3) * t154 + t156 * t183) * t133 + (Icges(4,5) * t164 + Icges(4,6) * t166) * t150;
t176 = t132 * t214 + t180;
t127 = t156 * pkin(2) + t154 * pkin(7);
t175 = t150 * t127 + t220 * V_base(4) + t198;
t174 = V_base(4) * t126 + (-t127 - t215) * V_base(5) + t197;
t77 = pkin(8) * t154 + t156 * t213;
t173 = -t133 * t214 + t150 * t77 + t175;
t172 = -t132 * t77 + t133 * t76 + t174;
t100 = Icges(4,6) * t154 + t156 * t186;
t101 = -Icges(4,5) * t156 + t154 * t189;
t102 = Icges(4,5) * t154 + t156 * t189;
t137 = Icges(4,2) * t166 + t206;
t140 = Icges(4,1) * t164 + t205;
t99 = -Icges(4,6) * t156 + t154 * t186;
t169 = (-t100 * t164 + t102 * t166) * t133 + (t101 * t166 - t164 * t99) * t132 + (-t137 * t164 + t140 * t166) * t150;
t157 = Icges(2,4) * t167;
t149 = Icges(3,4) * t156;
t145 = rSges(2,1) * t167 - rSges(2,2) * t165;
t144 = rSges(2,1) * t165 + rSges(2,2) * t167;
t143 = rSges(4,1) * t164 + rSges(4,2) * t166;
t142 = Icges(2,1) * t167 - t208;
t141 = Icges(2,1) * t165 + t157;
t139 = -Icges(2,2) * t165 + t157;
t138 = Icges(2,2) * t167 + t208;
t131 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t130 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t129 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t125 = rSges(3,1) * t156 - rSges(3,2) * t154;
t124 = rSges(3,1) * t154 + rSges(3,2) * t156;
t123 = rSges(5,1) * t153 + rSges(5,2) * t155;
t121 = Icges(3,1) * t156 - t207;
t120 = Icges(3,1) * t154 + t149;
t117 = -Icges(3,2) * t154 + t149;
t116 = Icges(3,2) * t156 + t207;
t104 = rSges(4,3) * t154 + t156 * t192;
t103 = -rSges(4,3) * t156 + t154 * t192;
t96 = V_base(5) * rSges(2,3) - t144 * t152 + t196;
t95 = t145 * t152 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t94 = t144 * V_base(4) - t145 * V_base(5) + V_base(3);
t93 = rSges(5,3) * t154 + t156 * t191;
t91 = -rSges(5,3) * t156 + t154 * t191;
t73 = V_base(5) * rSges(3,3) - t124 * t150 + t180;
t72 = t125 * t150 + (-rSges(3,3) + t220) * V_base(4) + t198;
t69 = t124 * V_base(4) + (-t125 - t215) * V_base(5) + t197;
t68 = t132 * t143 + (-t103 - t126) * t150 + t180;
t67 = t104 * t150 - t133 * t143 + t175;
t66 = t103 * t133 - t104 * t132 + t174;
t65 = t106 * t123 + (-t91 + t209) * t150 + t176;
t64 = -t107 * t123 + t150 * t93 + t173;
t63 = -t106 * t93 + t107 * t91 + t172;
t62 = qJD(5) * t154 + t193 * t106 + (t209 - t211) * t150 + t176;
t61 = -qJD(5) * t156 - t107 * t193 + t150 * t210 + t173;
t60 = -t106 * t210 + t107 * t211 + t172;
t1 = m(1) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(2) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t133 * (t177 * t154 + t169 * t156) / 0.2e1 + t132 * (t169 * t154 - t177 * t156) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + (t222 * t154 - t221 * t156) * t106 / 0.2e1 + (t221 * t154 + t222 * t156) * t107 / 0.2e1 + ((-t116 * t154 + t120 * t156 - t138 * t165 + t141 * t167 + Icges(1,4)) * V_base(5) + (-t154 * t117 + t156 * t121 - t165 * t139 + t167 * t142 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t156 * t116 + t154 * t120 + t167 * t138 + t165 * t141 + Icges(1,2)) * V_base(5) + (t117 * t156 + t121 * t154 + t139 * t167 + t142 * t165 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t166 + t102 * t164) * t133 + (t101 * t164 + t166 * t99) * t132 + (t229 * t153 + t231 * t155) * t107 + (t230 * t153 + t232 * t155) * t106 + (t166 * t137 + t164 * t140 + t226 * t153 + t227 * t155 + Icges(3,3)) * t150) * t150 / 0.2e1 + V_base(4) * t150 * (Icges(3,5) * t156 - Icges(3,6) * t154) + V_base(5) * t150 * (Icges(3,5) * t154 + Icges(3,6) * t156) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t165 + Icges(2,6) * t167) * V_base(5) + (Icges(2,5) * t167 - Icges(2,6) * t165) * V_base(4) + Icges(2,3) * t152 / 0.2e1) * t152;
T = t1;
