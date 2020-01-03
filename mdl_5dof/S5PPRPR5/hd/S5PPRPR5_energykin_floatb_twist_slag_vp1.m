% Calculate kinetic energy for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:24
% DurationCPUTime: 1.34s
% Computational Cost: add. (581->194), mult. (1006->235), div. (0->0), fcn. (986->6), ass. (0->98)
t213 = Icges(2,4) - Icges(3,5);
t212 = Icges(4,4) + Icges(5,6);
t211 = Icges(2,1) + Icges(3,1);
t210 = Icges(4,1) + Icges(5,2);
t209 = Icges(3,4) + Icges(2,5);
t208 = -Icges(5,4) + Icges(4,5);
t207 = Icges(5,5) - Icges(4,6);
t206 = Icges(2,2) + Icges(3,3);
t205 = Icges(4,2) + Icges(5,3);
t204 = Icges(2,6) - Icges(3,6);
t180 = sin(pkin(7));
t181 = cos(pkin(7));
t184 = sin(qJ(3));
t185 = cos(qJ(3));
t109 = -t180 * t185 + t181 * t184;
t203 = t212 * t109;
t202 = t213 * t180;
t201 = t213 * t181;
t108 = -t180 * t184 - t181 * t185;
t200 = t212 * t108;
t199 = t205 * t109 + t200;
t198 = -t205 * t108 + t203;
t197 = t210 * t108 + t203;
t196 = t210 * t109 - t200;
t193 = -t206 * t181 - t202;
t192 = t206 * t180 - t201;
t191 = t211 * t180 + t201;
t190 = t211 * t181 - t202;
t100 = -qJD(5) * t109 + V_base(5);
t101 = -qJD(5) * t108 + V_base(4);
t146 = sin(qJ(5));
t147 = cos(qJ(5));
t177 = Icges(6,4) * t147;
t132 = Icges(6,2) * t146 - t177;
t178 = Icges(6,4) * t146;
t133 = -Icges(6,1) * t147 + t178;
t143 = V_base(6) - qJD(3);
t158 = Icges(6,2) * t147 + t178;
t69 = -Icges(6,6) * t108 + t109 * t158;
t70 = -Icges(6,6) * t109 - t108 * t158;
t159 = Icges(6,1) * t146 + t177;
t71 = -Icges(6,5) * t108 + t109 * t159;
t72 = -Icges(6,5) * t109 - t108 * t159;
t187 = (t146 * t71 + t147 * t69) * t101 + (t146 * t72 + t147 * t70) * t100 + (t132 * t147 + t133 * t146) * t143;
t183 = pkin(6) * t108;
t182 = pkin(6) * t109;
t175 = V_base(5) * qJ(1) + V_base(1);
t171 = qJD(1) + V_base(3);
t170 = t181 * pkin(2);
t169 = t180 * pkin(2);
t166 = qJD(2) * t180 + t175;
t125 = pkin(1) * t180 - qJ(2) * t181;
t165 = V_base(4) * t125 + t171;
t128 = pkin(1) * t181 + qJ(2) * t180;
t164 = -t128 - t170;
t163 = V_base(4) * t169 + t165;
t162 = rSges(6,1) * t146 + rSges(6,2) * t147;
t157 = Icges(6,5) * t146 + Icges(6,6) * t147;
t96 = -pkin(3) * t108 + qJ(4) * t109;
t155 = -t96 + t164;
t94 = -pkin(3) * t109 - qJ(4) * t108;
t154 = V_base(4) * t94 + t163;
t153 = -qJD(2) * t181 + V_base(6) * t128 + V_base(2);
t152 = -(-Icges(6,3) * t109 - t108 * t157) * t100 - (-Icges(6,3) * t108 + t109 * t157) * t101 - (-Icges(6,5) * t147 + Icges(6,6) * t146) * t143;
t151 = V_base(4) * pkin(5) + V_base(6) * t170 + t153;
t150 = -qJD(4) * t108 + t143 * t96 + t151;
t149 = (-t169 - t125) * V_base(6) + t166;
t148 = qJD(4) * t109 + t149;
t134 = -rSges(6,1) * t147 + t146 * rSges(6,2);
t130 = rSges(2,1) * t181 - rSges(2,2) * t180;
t129 = rSges(3,1) * t181 + rSges(3,3) * t180;
t127 = rSges(2,1) * t180 + rSges(2,2) * t181;
t126 = rSges(3,1) * t180 - rSges(3,3) * t181;
t112 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t111 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t110 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t99 = V_base(5) * rSges(2,3) - t127 * V_base(6) + t175;
t98 = t130 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t97 = -rSges(4,1) * t108 - rSges(4,2) * t109;
t95 = -rSges(4,1) * t109 + rSges(4,2) * t108;
t93 = rSges(5,2) * t109 - rSges(5,3) * t108;
t92 = rSges(5,2) * t108 + rSges(5,3) * t109;
t91 = t127 * V_base(4) - t130 * V_base(5) + t171;
t76 = V_base(5) * rSges(3,2) + (-t125 - t126) * V_base(6) + t166;
t75 = V_base(6) * t129 + (-qJ(1) - rSges(3,2)) * V_base(4) + t153;
t74 = -t109 * rSges(6,3) - t108 * t162;
t73 = -t108 * rSges(6,3) + t109 * t162;
t66 = t126 * V_base(4) + (-t128 - t129) * V_base(5) + t165;
t65 = -t143 * t95 + (-pkin(5) - rSges(4,3)) * V_base(5) + t149;
t64 = t143 * t97 + (rSges(4,3) - qJ(1)) * V_base(4) + t151;
t63 = V_base(4) * t95 + (-t97 + t164) * V_base(5) + t163;
t62 = (-pkin(5) - rSges(5,1)) * V_base(5) + (-t93 - t94) * t143 + t148;
t61 = t143 * t92 + (rSges(5,1) - qJ(1)) * V_base(4) + t150;
t60 = V_base(4) * t93 + (-t92 + t155) * V_base(5) + t154;
t59 = t100 * t134 + (-pkin(4) - pkin(5)) * V_base(5) + (-t74 - t94 + t182) * t143 + t148;
t58 = -t101 * t134 + (pkin(4) - qJ(1)) * V_base(4) + (t73 - t183) * t143 + t150;
t57 = -V_base(4) * t182 - t100 * t73 + t101 * t74 + (t155 + t183) * V_base(5) + t154;
t1 = m(1) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(2) * (t91 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t101 * (t152 * t108 + t187 * t109) / 0.2e1 + t100 * (-t187 * t108 + t152 * t109) / 0.2e1 + ((t146 * t69 - t147 * t71) * t101 + (t146 * t70 - t147 * t72) * t100 + (t146 * t132 - t147 * t133 + Icges(5,1) + Icges(4,3)) * t143) * t143 / 0.2e1 + ((t108 * t196 + t109 * t198 + t180 * t193 + t181 * t191 + Icges(1,4)) * V_base(5) + (t197 * t108 + t199 * t109 + t192 * t180 + t190 * t181 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t198 * t108 + t196 * t109 + t191 * t180 - t193 * t181 + Icges(1,2)) * V_base(5) + (-t199 * t108 + t197 * t109 + t190 * t180 - t192 * t181 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t143 * (t207 * t108 + t208 * t109) + V_base(4) * t143 * (t208 * t108 - t207 * t109) + (Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + (t204 * V_base(5) + t209 * V_base(4)) * t181 + (-t204 * V_base(4) + t209 * V_base(5)) * t180) * V_base(6);
T = t1;
