% Calculate kinetic energy for
% S5RPRRP3
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:30
% EndTime: 2022-01-23 09:29:31
% DurationCPUTime: 1.96s
% Computational Cost: add. (1212->215), mult. (960->285), div. (0->0), fcn. (740->8), ass. (0->116)
t244 = Icges(5,4) + Icges(6,4);
t243 = Icges(5,1) + Icges(6,1);
t242 = Icges(5,2) + Icges(6,2);
t162 = qJ(3) + qJ(4);
t155 = cos(t162);
t241 = t244 * t155;
t154 = sin(t162);
t240 = t244 * t154;
t239 = Icges(5,5) + Icges(6,5);
t238 = Icges(5,6) + Icges(6,6);
t237 = -t242 * t154 + t241;
t236 = t243 * t155 - t240;
t235 = rSges(6,1) + pkin(4);
t161 = qJ(1) + pkin(8);
t151 = sin(t161);
t152 = cos(t161);
t234 = t237 * t151 - t238 * t152;
t233 = t238 * t151 + t237 * t152;
t232 = t236 * t151 - t239 * t152;
t231 = t239 * t151 + t236 * t152;
t230 = Icges(5,3) + Icges(6,3);
t229 = t242 * t155 + t240;
t228 = t243 * t154 + t241;
t227 = -t238 * t154 + t239 * t155;
t226 = rSges(6,3) + qJ(5);
t225 = -rSges(6,2) * t154 + t235 * t155;
t106 = V_base(5) + (-qJD(3) - qJD(4)) * t152;
t133 = qJD(3) * t151 + V_base(4);
t107 = qJD(4) * t151 + t133;
t153 = V_base(6) + qJD(1);
t222 = (-t229 * t154 + t228 * t155) * t153 + (-t233 * t154 + t231 * t155) * t107 + (-t234 * t154 + t232 * t155) * t106;
t221 = (t239 * t154 + t238 * t155) * t153 + (t230 * t151 + t227 * t152) * t107 + (t227 * t151 - t230 * t152) * t106;
t164 = sin(qJ(1));
t217 = pkin(1) * t164;
t166 = cos(qJ(1));
t216 = pkin(1) * t166;
t163 = sin(qJ(3));
t215 = pkin(3) * t163;
t165 = cos(qJ(3));
t214 = t165 * pkin(3);
t213 = -pkin(5) - qJ(2);
t211 = t225 * t151 - t226 * t152;
t210 = t226 * t151 + t225 * t152;
t209 = Icges(2,4) * t164;
t208 = Icges(3,4) * t151;
t207 = Icges(4,4) * t163;
t206 = Icges(4,4) * t165;
t199 = t153 * t216 + V_base(2);
t198 = V_base(5) * pkin(5) + V_base(1);
t118 = t151 * pkin(2) - t152 * pkin(6);
t195 = -t118 - t217;
t194 = rSges(6,2) * t155 + t235 * t154;
t193 = V_base(5) * qJ(2) + t198;
t192 = V_base(4) * t217 + qJD(2) + V_base(3);
t76 = -pkin(7) * t152 + t214 * t151;
t191 = t195 - t76;
t132 = -qJD(3) * t152 + V_base(5);
t190 = t132 * t215 + t193;
t189 = rSges(4,1) * t165 - rSges(4,2) * t163;
t188 = rSges(5,1) * t155 - rSges(5,2) * t154;
t186 = Icges(4,1) * t165 - t207;
t183 = -Icges(4,2) * t163 + t206;
t180 = Icges(4,5) * t165 - Icges(4,6) * t163;
t175 = (-Icges(4,3) * t152 + t180 * t151) * t132 + (Icges(4,3) * t151 + t180 * t152) * t133 + (Icges(4,5) * t163 + Icges(4,6) * t165) * t153;
t119 = t152 * pkin(2) + t151 * pkin(6);
t174 = t153 * t119 + t213 * V_base(4) + t199;
t173 = V_base(4) * t118 + (-t119 - t216) * V_base(5) + t192;
t77 = pkin(7) * t151 + t214 * t152;
t172 = -t133 * t215 + t153 * t77 + t174;
t171 = -t132 * t77 + t133 * t76 + t173;
t100 = Icges(4,5) * t151 + t186 * t152;
t137 = Icges(4,2) * t165 + t207;
t140 = Icges(4,1) * t163 + t206;
t97 = -Icges(4,6) * t152 + t183 * t151;
t98 = Icges(4,6) * t151 + t183 * t152;
t99 = -Icges(4,5) * t152 + t186 * t151;
t168 = (t100 * t165 - t163 * t98) * t133 + (-t163 * t97 + t165 * t99) * t132 + (-t137 * t163 + t140 * t165) * t153;
t157 = Icges(2,4) * t166;
t149 = Icges(3,4) * t152;
t145 = rSges(2,1) * t166 - rSges(2,2) * t164;
t144 = rSges(2,1) * t164 + rSges(2,2) * t166;
t143 = rSges(4,1) * t163 + rSges(4,2) * t165;
t142 = Icges(2,1) * t166 - t209;
t141 = Icges(2,1) * t164 + t157;
t139 = -Icges(2,2) * t164 + t157;
t138 = Icges(2,2) * t166 + t209;
t131 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t130 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t129 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t127 = rSges(5,1) * t154 + rSges(5,2) * t155;
t117 = rSges(3,1) * t152 - rSges(3,2) * t151;
t116 = rSges(3,1) * t151 + rSges(3,2) * t152;
t115 = Icges(3,1) * t152 - t208;
t114 = Icges(3,1) * t151 + t149;
t113 = -Icges(3,2) * t151 + t149;
t112 = Icges(3,2) * t152 + t208;
t104 = rSges(4,3) * t151 + t189 * t152;
t103 = -rSges(4,3) * t152 + t189 * t151;
t102 = V_base(5) * rSges(2,3) - t144 * t153 + t198;
t101 = t145 * t153 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t94 = t144 * V_base(4) - t145 * V_base(5) + V_base(3);
t93 = rSges(5,3) * t151 + t188 * t152;
t91 = -rSges(5,3) * t152 + t188 * t151;
t73 = V_base(5) * rSges(3,3) + (-t116 - t217) * t153 + t193;
t72 = t117 * t153 + (-rSges(3,3) + t213) * V_base(4) + t199;
t69 = t116 * V_base(4) + (-t117 - t216) * V_base(5) + t192;
t68 = t132 * t143 + (-t103 + t195) * t153 + t193;
t67 = t104 * t153 - t133 * t143 + t174;
t66 = t103 * t133 - t104 * t132 + t173;
t65 = t106 * t127 + (t191 - t91) * t153 + t190;
t64 = -t107 * t127 + t153 * t93 + t172;
t63 = -t106 * t93 + t107 * t91 + t171;
t62 = qJD(5) * t151 + t194 * t106 + (t191 - t211) * t153 + t190;
t61 = -qJD(5) * t152 - t194 * t107 + t210 * t153 + t172;
t60 = -t210 * t106 + t211 * t107 + t171;
t1 = m(1) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t94 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t133 * (t175 * t151 + t168 * t152) / 0.2e1 + t132 * (t168 * t151 - t175 * t152) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + (t222 * t151 - t221 * t152) * t106 / 0.2e1 + (t221 * t151 + t222 * t152) * t107 / 0.2e1 + ((-t112 * t151 + t114 * t152 - t138 * t164 + t141 * t166 + Icges(1,4)) * V_base(5) + (-t151 * t113 + t152 * t115 - t164 * t139 + t166 * t142 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t152 * t112 + t151 * t114 + t166 * t138 + t164 * t141 + Icges(1,2)) * V_base(5) + (t113 * t152 + t115 * t151 + t139 * t166 + t142 * t164 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t163 + t165 * t98) * t133 + (t163 * t99 + t165 * t97) * t132 + (t231 * t154 + t233 * t155) * t107 + (t232 * t154 + t234 * t155) * t106 + (t165 * t137 + t163 * t140 + t228 * t154 + t229 * t155 + Icges(2,3) + Icges(3,3)) * t153) * t153 / 0.2e1 + t153 * V_base(5) * (Icges(2,5) * t164 + Icges(3,5) * t151 + Icges(2,6) * t166 + Icges(3,6) * t152) + t153 * V_base(4) * (Icges(2,5) * t166 + Icges(3,5) * t152 - Icges(2,6) * t164 - Icges(3,6) * t151) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
