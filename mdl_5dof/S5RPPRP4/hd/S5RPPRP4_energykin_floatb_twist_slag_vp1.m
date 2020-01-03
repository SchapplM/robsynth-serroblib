% Calculate kinetic energy for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:00
% EndTime: 2019-12-31 17:52:02
% DurationCPUTime: 1.81s
% Computational Cost: add. (748->190), mult. (1282->231), div. (0->0), fcn. (1326->6), ass. (0->100)
t236 = Icges(5,4) + Icges(6,4);
t235 = Icges(5,1) + Icges(6,1);
t234 = Icges(5,2) + Icges(6,2);
t154 = cos(qJ(4));
t233 = t236 * t154;
t153 = sin(qJ(4));
t232 = t236 * t153;
t231 = Icges(5,5) + Icges(6,5);
t230 = Icges(5,6) + Icges(6,6);
t229 = t234 * t153 - t233;
t228 = t235 * t154 - t232;
t227 = rSges(6,1) + pkin(4);
t226 = Icges(2,4) - Icges(3,5);
t191 = sin(pkin(7));
t192 = cos(pkin(7));
t197 = sin(qJ(1));
t198 = cos(qJ(1));
t108 = -t197 * t191 - t198 * t192;
t109 = t198 * t191 - t197 * t192;
t225 = -t230 * t108 + t229 * t109;
t224 = t229 * t108 + t230 * t109;
t223 = t231 * t108 + t228 * t109;
t222 = t228 * t108 - t231 * t109;
t221 = Icges(2,1) + Icges(3,1);
t220 = Icges(3,4) + Icges(2,5);
t219 = Icges(2,2) + Icges(3,3);
t218 = Icges(2,6) - Icges(3,6);
t217 = Icges(5,3) + Icges(6,3);
t216 = -t234 * t154 - t232;
t215 = t235 * t153 + t233;
t214 = t230 * t153 - t231 * t154;
t213 = rSges(6,3) + qJ(5);
t212 = t226 * t197;
t211 = t226 * t198;
t210 = rSges(6,2) * t153 - t227 * t154;
t209 = -t219 * t198 - t212;
t208 = t219 * t197 - t211;
t207 = t221 * t197 + t211;
t206 = t221 * t198 - t212;
t104 = -qJD(4) * t108 + V_base(5);
t105 = qJD(4) * t109 + V_base(4);
t146 = V_base(6) + qJD(1);
t205 = (t216 * t153 + t215 * t154) * t146 + (t224 * t153 + t222 * t154) * t105 + (t225 * t153 + t223 * t154) * t104;
t204 = (-t231 * t153 - t230 * t154) * t146 + (t214 * t108 + t217 * t109) * t105 + (-t217 * t108 + t214 * t109) * t104;
t194 = -t213 * t108 + t210 * t109;
t193 = t210 * t108 + t213 * t109;
t190 = Icges(4,4) * t108;
t134 = t197 * pkin(1) - t198 * qJ(2);
t185 = V_base(4) * t134 + V_base(3);
t184 = V_base(5) * pkin(5) + V_base(1);
t181 = t198 * pkin(2);
t180 = t197 * pkin(2);
t177 = rSges(6,2) * t154 + t227 * t153;
t176 = qJD(2) * t197 + t184;
t175 = -t134 - t180;
t137 = t198 * pkin(1) + t197 * qJ(2);
t174 = -t137 - t181;
t173 = V_base(4) * t180 - qJD(3) + t185;
t172 = -rSges(5,1) * t154 + rSges(5,2) * t153;
t99 = -pkin(3) * t109 - pkin(6) * t108;
t164 = -t99 + t175;
t163 = -qJD(2) * t198 + t146 * t137 + V_base(2);
t162 = -V_base(5) * qJ(3) + t176;
t159 = V_base(4) * qJ(3) + t146 * t181 + t163;
t100 = -pkin(3) * t108 + pkin(6) * t109;
t158 = -V_base(4) * pkin(5) + t146 * t100 + t159;
t157 = V_base(4) * t99 + (-t100 + t174) * V_base(5) + t173;
t139 = t198 * rSges(2,1) - t197 * rSges(2,2);
t138 = t198 * rSges(3,1) + t197 * rSges(3,3);
t136 = t197 * rSges(2,1) + t198 * rSges(2,2);
t135 = t197 * rSges(3,1) - t198 * rSges(3,3);
t133 = -t153 * rSges(5,1) - rSges(5,2) * t154;
t113 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t112 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t111 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t106 = Icges(4,4) * t109;
t103 = V_base(5) * rSges(2,3) - t136 * t146 + t184;
t102 = t139 * t146 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t101 = t136 * V_base(4) - t139 * V_base(5) + V_base(3);
t98 = -rSges(4,1) * t108 - rSges(4,2) * t109;
t97 = -rSges(4,1) * t109 + rSges(4,2) * t108;
t96 = -Icges(4,1) * t108 - t106;
t95 = -Icges(4,1) * t109 + t190;
t94 = -Icges(4,2) * t109 - t190;
t93 = Icges(4,2) * t108 - t106;
t88 = V_base(5) * rSges(3,2) + (-t134 - t135) * t146 + t176;
t87 = t146 * t138 + (-pkin(5) - rSges(3,2)) * V_base(4) + t163;
t86 = t109 * rSges(5,3) + t172 * t108;
t84 = -t108 * rSges(5,3) + t172 * t109;
t70 = t135 * V_base(4) + (-t137 - t138) * V_base(5) + t185;
t67 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t97 + t175) * t146 + t176;
t66 = t146 * t98 + (rSges(4,3) - pkin(5)) * V_base(4) + t159;
t65 = V_base(4) * t97 + (-t98 + t174) * V_base(5) + t173;
t64 = t104 * t133 + (-t84 + t164) * t146 + t162;
t63 = -t105 * t133 + t146 * t86 + t158;
t62 = -t104 * t86 + t105 * t84 + t157;
t61 = qJD(5) * t109 - t177 * t104 + (t164 - t194) * t146 + t162;
t60 = -qJD(5) * t108 + t177 * t105 + t193 * t146 + t158;
t59 = -t193 * t104 + t194 * t105 + t157;
t1 = m(1) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(3) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (-t204 * t108 + t205 * t109) * t104 / 0.2e1 + (t205 * t108 + t204 * t109) * t105 / 0.2e1 + ((-t108 * t95 - t109 * t93 + t197 * t209 + t207 * t198 + Icges(1,4)) * V_base(5) + (-t108 * t96 - t109 * t94 + t208 * t197 + t206 * t198 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t108 * t93 - t109 * t95 + t207 * t197 - t209 * t198 + Icges(1,2)) * V_base(5) + (t108 * t94 - t109 * t96 + t197 * t206 - t198 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t222 * t153 - t224 * t154) * t105 + (t223 * t153 - t225 * t154) * t104 + (t215 * t153 - t216 * t154 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t146) * t146 / 0.2e1 + t146 * V_base(5) * (Icges(4,5) * t109 - Icges(4,6) * t108 + t220 * t197 + t218 * t198) + t146 * V_base(4) * (Icges(4,5) * t108 + Icges(4,6) * t109 - t218 * t197 + t220 * t198) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
