% Calculate kinetic energy for
% S5RRRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:05
% EndTime: 2019-12-31 20:53:06
% DurationCPUTime: 1.78s
% Computational Cost: add. (953->187), mult. (925->243), div. (0->0), fcn. (705->6), ass. (0->100)
t241 = Icges(4,4) + Icges(5,6) - Icges(6,6);
t240 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t239 = -Icges(4,2) - Icges(6,2) - Icges(5,3);
t167 = cos(qJ(3));
t238 = t241 * t167;
t165 = sin(qJ(3));
t237 = t241 * t165;
t236 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t235 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t234 = t240 * t167 - t237;
t233 = t239 * t165 + t238;
t232 = rSges(6,3) + qJ(5);
t164 = qJ(1) + qJ(2);
t159 = sin(t164);
t160 = cos(t164);
t231 = t233 * t159 + t235 * t160;
t230 = -t235 * t159 + t233 * t160;
t229 = t234 * t159 + t236 * t160;
t228 = -t236 * t159 + t234 * t160;
t227 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t226 = t240 * t165 + t238;
t225 = t239 * t167 - t237;
t224 = t235 * t165 - t236 * t167;
t223 = rSges(6,1) + pkin(4);
t222 = rSges(6,2) * t165 + t232 * t167;
t129 = -qJD(3) * t160 + V_base(5);
t130 = qJD(3) * t159 + V_base(4);
t158 = V_base(6) + qJD(1);
t157 = qJD(2) + t158;
t221 = (t225 * t165 + t226 * t167) * t157 + (-t230 * t165 + t228 * t167) * t130 + (-t231 * t165 + t229 * t167) * t129;
t220 = (-t236 * t165 - t235 * t167) * t157 + (t227 * t159 + t224 * t160) * t130 + (t224 * t159 - t227 * t160) * t129;
t219 = -pkin(5) - pkin(6);
t166 = sin(qJ(1));
t215 = pkin(1) * t166;
t168 = cos(qJ(1));
t214 = pkin(1) * t168;
t213 = Icges(2,4) * t166;
t212 = Icges(3,4) * t159;
t204 = t223 * t159 + t222 * t160;
t203 = t222 * t159 - t223 * t160;
t190 = pkin(3) * t167 + qJ(4) * t165;
t111 = t190 * t159;
t124 = pkin(2) * t159 - pkin(7) * t160;
t202 = -t111 - t124;
t201 = qJD(4) * t165;
t200 = qJD(5) * t167;
t199 = t158 * t214 + V_base(2);
t198 = V_base(4) * t215 + V_base(3);
t197 = V_base(5) * pkin(5) + V_base(1);
t194 = -rSges(6,2) * t167 + t232 * t165;
t193 = rSges(4,1) * t167 - rSges(4,2) * t165;
t192 = -rSges(5,2) * t167 + rSges(5,3) * t165;
t180 = V_base(5) * pkin(6) - t158 * t215 + t197;
t125 = pkin(2) * t160 + pkin(7) * t159;
t176 = t157 * t125 + t219 * V_base(4) + t199;
t146 = pkin(3) * t165 - qJ(4) * t167;
t175 = t129 * t146 + t160 * t201 + t180;
t174 = V_base(4) * t124 + (-t125 - t214) * V_base(5) + t198;
t112 = t190 * t160;
t173 = t157 * t112 + t159 * t201 + t176;
t172 = -qJD(4) * t167 + t130 * t111 + t174;
t161 = Icges(2,4) * t168;
t156 = Icges(3,4) * t160;
t151 = rSges(2,1) * t168 - t166 * rSges(2,2);
t149 = t166 * rSges(2,1) + rSges(2,2) * t168;
t148 = rSges(4,1) * t165 + rSges(4,2) * t167;
t147 = -rSges(5,2) * t165 - rSges(5,3) * t167;
t145 = Icges(2,1) * t168 - t213;
t144 = Icges(2,1) * t166 + t161;
t142 = -Icges(2,2) * t166 + t161;
t141 = Icges(2,2) * t168 + t213;
t128 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t127 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t126 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t123 = rSges(3,1) * t160 - rSges(3,2) * t159;
t122 = rSges(3,1) * t159 + rSges(3,2) * t160;
t121 = Icges(3,1) * t160 - t212;
t120 = Icges(3,1) * t159 + t156;
t119 = -Icges(3,2) * t159 + t156;
t118 = Icges(3,2) * t160 + t212;
t108 = -rSges(5,1) * t160 + t192 * t159;
t106 = rSges(5,1) * t159 + t192 * t160;
t104 = rSges(4,3) * t159 + t193 * t160;
t103 = -rSges(4,3) * t160 + t193 * t159;
t83 = V_base(5) * rSges(2,3) - t149 * t158 + t197;
t82 = t151 * t158 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t81 = t149 * V_base(4) - t151 * V_base(5) + V_base(3);
t79 = V_base(5) * rSges(3,3) - t122 * t157 + t180;
t78 = t123 * t157 + (-rSges(3,3) + t219) * V_base(4) + t199;
t77 = V_base(4) * t122 + (-t123 - t214) * V_base(5) + t198;
t76 = t129 * t148 + (-t103 - t124) * t157 + t180;
t75 = t104 * t157 - t130 * t148 + t176;
t74 = t130 * t103 - t129 * t104 + t174;
t73 = t129 * t147 + (-t108 + t202) * t157 + t175;
t72 = t106 * t157 + (-t146 - t147) * t130 + t173;
t71 = t130 * t108 + (-t106 - t112) * t129 + t172;
t70 = t160 * t200 + t194 * t129 + (t202 - t203) * t157 + t175;
t69 = t159 * t200 + t204 * t157 + (-t146 - t194) * t130 + t173;
t68 = qJD(5) * t165 + t203 * t130 + (-t112 - t204) * t129 + t172;
t1 = m(1) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(3) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (t221 * t159 - t220 * t160) * t129 / 0.2e1 + (t220 * t159 + t221 * t160) * t130 / 0.2e1 + ((-t118 * t159 + t120 * t160 - t166 * t141 + t144 * t168 + Icges(1,4)) * V_base(5) + (-t119 * t159 + t121 * t160 - t166 * t142 + t145 * t168 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t118 * t160 + t159 * t120 + t141 * t168 + t166 * t144 + Icges(1,2)) * V_base(5) + (t119 * t160 + t121 * t159 + t142 * t168 + t166 * t145 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t228 * t165 + t230 * t167) * t130 + (t229 * t165 + t231 * t167) * t129 + (t226 * t165 - t225 * t167 + Icges(3,3)) * t157) * t157 / 0.2e1 + V_base(4) * t157 * (Icges(3,5) * t160 - Icges(3,6) * t159) + V_base(5) * t157 * (Icges(3,5) * t159 + Icges(3,6) * t160) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t166 + Icges(2,6) * t168) * V_base(5) + (Icges(2,5) * t168 - Icges(2,6) * t166) * V_base(4) + Icges(2,3) * t158 / 0.2e1) * t158;
T = t1;
