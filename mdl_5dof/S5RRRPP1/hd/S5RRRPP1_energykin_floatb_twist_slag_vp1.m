% Calculate kinetic energy for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:20
% EndTime: 2019-12-31 20:49:22
% DurationCPUTime: 1.81s
% Computational Cost: add. (1168->211), mult. (929->271), div. (0->0), fcn. (709->8), ass. (0->112)
t242 = Icges(5,4) - Icges(6,5);
t241 = Icges(5,1) + Icges(6,1);
t240 = Icges(5,2) + Icges(6,3);
t162 = qJ(3) + pkin(8);
t155 = cos(t162);
t239 = t242 * t155;
t154 = sin(t162);
t238 = t242 * t154;
t237 = Icges(6,4) + Icges(5,5);
t236 = Icges(5,6) - Icges(6,6);
t235 = t240 * t154 - t239;
t234 = t241 * t155 - t238;
t233 = rSges(6,1) + pkin(4);
t232 = rSges(6,3) + qJ(5);
t163 = qJ(1) + qJ(2);
t157 = sin(t163);
t158 = cos(t163);
t231 = t235 * t157 + t236 * t158;
t230 = -t236 * t157 + t235 * t158;
t229 = t234 * t157 - t237 * t158;
t228 = t237 * t157 + t234 * t158;
t227 = -t240 * t155 - t238;
t226 = t241 * t154 + t239;
t225 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t165 = sin(qJ(3));
t167 = cos(qJ(3));
t224 = Icges(4,5) * t167 - Icges(4,6) * t165 - t236 * t154 + t237 * t155;
t223 = t232 * t154 + t233 * t155;
t205 = Icges(4,4) * t167;
t186 = -Icges(4,2) * t165 + t205;
t101 = -Icges(4,6) * t158 + t186 * t157;
t102 = Icges(4,6) * t157 + t186 * t158;
t206 = Icges(4,4) * t165;
t189 = Icges(4,1) * t167 - t206;
t103 = -Icges(4,5) * t158 + t189 * t157;
t104 = Icges(4,5) * t157 + t189 * t158;
t134 = -qJD(3) * t158 + V_base(5);
t135 = qJD(3) * t157 + V_base(4);
t139 = Icges(4,2) * t167 + t206;
t142 = Icges(4,1) * t165 + t205;
t156 = V_base(6) + qJD(1);
t152 = qJD(2) + t156;
t222 = (-t139 * t165 + t142 * t167 + t227 * t154 + t226 * t155) * t152 + (-t102 * t165 + t104 * t167 + t230 * t154 + t228 * t155) * t135 + (-t101 * t165 + t103 * t167 + t231 * t154 + t229 * t155) * t134;
t221 = (Icges(4,5) * t165 + Icges(4,6) * t167 + t237 * t154 + t236 * t155) * t152 + (t225 * t157 + t224 * t158) * t135 + (t224 * t157 - t225 * t158) * t134;
t220 = -pkin(5) - pkin(6);
t166 = sin(qJ(1));
t216 = pkin(1) * t166;
t168 = cos(qJ(1));
t215 = pkin(1) * t168;
t214 = pkin(3) * t165;
t213 = pkin(3) * t167;
t211 = -rSges(6,2) * t158 + t223 * t157;
t210 = rSges(6,2) * t157 + t223 * t158;
t129 = pkin(2) * t157 - pkin(7) * t158;
t78 = -qJ(4) * t158 + t213 * t157;
t209 = -t129 - t78;
t208 = Icges(2,4) * t166;
t207 = Icges(3,4) * t157;
t200 = t233 * t154 - t232 * t155;
t199 = qJD(5) * t154;
t198 = t156 * t215 + V_base(2);
t197 = V_base(4) * t216 + V_base(3);
t196 = V_base(5) * pkin(5) + V_base(1);
t193 = rSges(4,1) * t167 - rSges(4,2) * t165;
t192 = rSges(5,1) * t155 - rSges(5,2) * t154;
t180 = V_base(5) * pkin(6) - t156 * t216 + t196;
t130 = pkin(2) * t158 + pkin(7) * t157;
t176 = t152 * t130 + t220 * V_base(4) + t198;
t175 = qJD(4) * t157 + t134 * t214 + t180;
t174 = V_base(4) * t129 + (-t130 - t215) * V_base(5) + t197;
t173 = t135 * t78 + t174;
t79 = qJ(4) * t157 + t213 * t158;
t172 = -qJD(4) * t158 + t152 * t79 + t176;
t159 = Icges(2,4) * t168;
t151 = Icges(3,4) * t158;
t147 = rSges(2,1) * t168 - t166 * rSges(2,2);
t146 = t166 * rSges(2,1) + rSges(2,2) * t168;
t145 = rSges(4,1) * t165 + rSges(4,2) * t167;
t144 = Icges(2,1) * t168 - t208;
t143 = Icges(2,1) * t166 + t159;
t141 = -Icges(2,2) * t166 + t159;
t140 = Icges(2,2) * t168 + t208;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = rSges(3,1) * t158 - rSges(3,2) * t157;
t127 = rSges(3,1) * t157 + rSges(3,2) * t158;
t126 = Icges(3,1) * t158 - t207;
t125 = Icges(3,1) * t157 + t151;
t124 = -Icges(3,2) * t157 + t151;
t123 = Icges(3,2) * t158 + t207;
t120 = rSges(5,1) * t154 + rSges(5,2) * t155;
t106 = rSges(4,3) * t157 + t193 * t158;
t105 = -rSges(4,3) * t158 + t193 * t157;
t98 = V_base(5) * rSges(2,3) - t146 * t156 + t196;
t97 = t147 * t156 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t96 = t146 * V_base(4) - t147 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t157 + t192 * t158;
t93 = -rSges(5,3) * t158 + t192 * t157;
t75 = V_base(5) * rSges(3,3) - t127 * t152 + t180;
t74 = t128 * t152 + (-rSges(3,3) + t220) * V_base(4) + t198;
t73 = V_base(4) * t127 + (-t128 - t215) * V_base(5) + t197;
t72 = t134 * t145 + (-t105 - t129) * t152 + t180;
t71 = t106 * t152 - t135 * t145 + t176;
t70 = t135 * t105 - t134 * t106 + t174;
t69 = t120 * t134 + (-t93 + t209) * t152 + t175;
t68 = t152 * t95 + (-t120 - t214) * t135 + t172;
t67 = t135 * t93 + (-t79 - t95) * t134 + t173;
t66 = t158 * t199 + t200 * t134 + (t209 - t211) * t152 + t175;
t65 = t157 * t199 + t210 * t152 + (-t200 - t214) * t135 + t172;
t64 = -qJD(5) * t155 + t211 * t135 + (-t79 - t210) * t134 + t173;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t222 * t157 - t221 * t158) * t134 / 0.2e1 + (t221 * t157 + t222 * t158) * t135 / 0.2e1 + ((-t123 * t157 + t125 * t158 - t166 * t140 + t143 * t168 + Icges(1,4)) * V_base(5) + (-t157 * t124 + t158 * t126 - t166 * t141 + t168 * t144 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t158 * t123 + t157 * t125 + t168 * t140 + t166 * t143 + Icges(1,2)) * V_base(5) + (t124 * t158 + t126 * t157 + t141 * t168 + t166 * t144 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t102 * t167 + t104 * t165 + t228 * t154 - t230 * t155) * t135 + (t101 * t167 + t103 * t165 + t229 * t154 - t231 * t155) * t134 + (t167 * t139 + t165 * t142 + t226 * t154 - t227 * t155 + Icges(3,3)) * t152) * t152 / 0.2e1 + V_base(4) * t152 * (Icges(3,5) * t158 - Icges(3,6) * t157) + V_base(5) * t152 * (Icges(3,5) * t157 + Icges(3,6) * t158) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t166 + Icges(2,6) * t168) * V_base(5) + (Icges(2,5) * t168 - Icges(2,6) * t166) * V_base(4) + Icges(2,3) * t156 / 0.2e1) * t156;
T = t1;
