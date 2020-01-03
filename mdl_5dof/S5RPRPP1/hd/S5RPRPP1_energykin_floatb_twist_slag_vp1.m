% Calculate kinetic energy for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:35
% EndTime: 2019-12-31 18:08:37
% DurationCPUTime: 1.74s
% Computational Cost: add. (1136->212), mult. (929->265), div. (0->0), fcn. (709->8), ass. (0->112)
t244 = Icges(5,4) - Icges(6,5);
t243 = Icges(5,1) + Icges(6,1);
t242 = Icges(5,2) + Icges(6,3);
t161 = qJ(3) + pkin(8);
t155 = cos(t161);
t241 = t244 * t155;
t153 = sin(t161);
t240 = t244 * t153;
t239 = Icges(6,4) + Icges(5,5);
t238 = Icges(5,6) - Icges(6,6);
t237 = t242 * t153 - t241;
t236 = t243 * t155 - t240;
t235 = rSges(6,1) + pkin(4);
t234 = rSges(6,3) + qJ(5);
t162 = qJ(1) + pkin(7);
t154 = sin(t162);
t156 = cos(t162);
t233 = t237 * t154 + t238 * t156;
t232 = -t238 * t154 + t237 * t156;
t231 = t236 * t154 - t239 * t156;
t230 = t239 * t154 + t236 * t156;
t229 = -t242 * t155 - t240;
t228 = t243 * t153 + t241;
t227 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t226 = Icges(4,5) * t166 - Icges(4,6) * t164 - t238 * t153 + t239 * t155;
t225 = t234 * t153 + t235 * t155;
t206 = Icges(4,4) * t166;
t184 = -Icges(4,2) * t164 + t206;
t100 = Icges(4,6) * t154 + t156 * t184;
t207 = Icges(4,4) * t164;
t187 = Icges(4,1) * t166 - t207;
t101 = -Icges(4,5) * t156 + t154 * t187;
t102 = Icges(4,5) * t154 + t156 * t187;
t134 = -qJD(3) * t156 + V_base(5);
t135 = qJD(3) * t154 + V_base(4);
t139 = Icges(4,2) * t166 + t207;
t142 = Icges(4,1) * t164 + t206;
t157 = V_base(6) + qJD(1);
t99 = -Icges(4,6) * t156 + t154 * t184;
t222 = (-t139 * t164 + t142 * t166 + t229 * t153 + t228 * t155) * t157 + (-t100 * t164 + t102 * t166 + t232 * t153 + t230 * t155) * t135 + (t101 * t166 + t233 * t153 + t231 * t155 - t164 * t99) * t134;
t221 = (Icges(4,5) * t164 + Icges(4,6) * t166 + t239 * t153 + t238 * t155) * t157 + (t227 * t154 + t226 * t156) * t135 + (t226 * t154 - t227 * t156) * t134;
t165 = sin(qJ(1));
t217 = pkin(1) * t165;
t167 = cos(qJ(1));
t216 = pkin(1) * t167;
t215 = pkin(3) * t164;
t214 = pkin(3) * t166;
t213 = -pkin(5) - qJ(2);
t211 = -rSges(6,2) * t156 + t225 * t154;
t210 = rSges(6,2) * t154 + t225 * t156;
t209 = Icges(2,4) * t165;
t208 = Icges(3,4) * t154;
t201 = t235 * t153 - t234 * t155;
t200 = qJD(5) * t153;
t199 = t157 * t216 + V_base(2);
t198 = V_base(5) * pkin(5) + V_base(1);
t129 = pkin(2) * t154 - pkin(6) * t156;
t195 = -t129 - t217;
t194 = V_base(5) * qJ(2) + t198;
t193 = V_base(4) * t217 + qJD(2) + V_base(3);
t78 = -qJ(4) * t156 + t154 * t214;
t192 = t195 - t78;
t191 = rSges(4,1) * t166 - rSges(4,2) * t164;
t190 = rSges(5,1) * t155 - rSges(5,2) * t153;
t178 = qJD(4) * t154 + t134 * t215 + t194;
t130 = pkin(2) * t156 + pkin(6) * t154;
t174 = t157 * t130 + t213 * V_base(4) + t199;
t173 = V_base(4) * t129 + (-t130 - t216) * V_base(5) + t193;
t172 = t135 * t78 + t173;
t79 = qJ(4) * t154 + t156 * t214;
t171 = -qJD(4) * t156 + t157 * t79 + t174;
t159 = Icges(2,4) * t167;
t151 = Icges(3,4) * t156;
t147 = rSges(2,1) * t167 - t165 * rSges(2,2);
t146 = t165 * rSges(2,1) + rSges(2,2) * t167;
t145 = rSges(4,1) * t164 + rSges(4,2) * t166;
t144 = Icges(2,1) * t167 - t209;
t143 = Icges(2,1) * t165 + t159;
t141 = -Icges(2,2) * t165 + t159;
t140 = Icges(2,2) * t167 + t209;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = rSges(3,1) * t156 - rSges(3,2) * t154;
t127 = rSges(3,1) * t154 + rSges(3,2) * t156;
t126 = rSges(5,1) * t153 + rSges(5,2) * t155;
t123 = Icges(3,1) * t156 - t208;
t122 = Icges(3,1) * t154 + t151;
t119 = -Icges(3,2) * t154 + t151;
t118 = Icges(3,2) * t156 + t208;
t106 = rSges(4,3) * t154 + t156 * t191;
t105 = -rSges(4,3) * t156 + t154 * t191;
t104 = V_base(5) * rSges(2,3) - t146 * t157 + t198;
t103 = t147 * t157 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t96 = t146 * V_base(4) - t147 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t154 + t156 * t190;
t93 = -rSges(5,3) * t156 + t154 * t190;
t76 = V_base(5) * rSges(3,3) + (-t127 - t217) * t157 + t194;
t75 = t128 * t157 + (-rSges(3,3) + t213) * V_base(4) + t199;
t73 = V_base(4) * t127 + (-t128 - t216) * V_base(5) + t193;
t72 = t134 * t145 + (-t105 + t195) * t157 + t194;
t71 = t106 * t157 - t135 * t145 + t174;
t70 = t135 * t105 - t134 * t106 + t173;
t69 = t126 * t134 + (t192 - t93) * t157 + t178;
t68 = t157 * t95 + (-t126 - t215) * t135 + t171;
t67 = t135 * t93 + (-t79 - t95) * t134 + t172;
t66 = t156 * t200 + t201 * t134 + (t192 - t211) * t157 + t178;
t65 = t154 * t200 + t210 * t157 + (-t201 - t215) * t135 + t171;
t64 = -qJD(5) * t155 + t211 * t135 + (-t79 - t210) * t134 + t172;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t104 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t222 * t154 - t221 * t156) * t134 / 0.2e1 + (t221 * t154 + t222 * t156) * t135 / 0.2e1 + ((-t118 * t154 + t122 * t156 - t165 * t140 + t143 * t167 + Icges(1,4)) * V_base(5) + (-t154 * t119 + t156 * t123 - t165 * t141 + t167 * t144 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t156 * t118 + t154 * t122 + t167 * t140 + t165 * t143 + Icges(1,2)) * V_base(5) + (t119 * t156 + t123 * t154 + t141 * t167 + t165 * t144 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t166 + t102 * t164 + t230 * t153 - t232 * t155) * t135 + (t101 * t164 + t231 * t153 - t233 * t155 + t166 * t99) * t134 + (t166 * t139 + t164 * t142 + t228 * t153 - t229 * t155 + Icges(2,3) + Icges(3,3)) * t157) * t157 / 0.2e1 + t157 * V_base(5) * (Icges(2,5) * t165 + Icges(3,5) * t154 + Icges(2,6) * t167 + Icges(3,6) * t156) + t157 * V_base(4) * (Icges(2,5) * t167 + Icges(3,5) * t156 - Icges(2,6) * t165 - Icges(3,6) * t154) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
