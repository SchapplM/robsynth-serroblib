% Calculate kinetic energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:57
% EndTime: 2019-12-31 16:45:58
% DurationCPUTime: 1.22s
% Computational Cost: add. (393->149), mult. (646->187), div. (0->0), fcn. (478->4), ass. (0->79)
t190 = Icges(4,4) + Icges(5,4);
t189 = Icges(4,1) + Icges(5,1);
t188 = Icges(4,2) + Icges(5,2);
t108 = cos(qJ(3));
t187 = t190 * t108;
t106 = sin(qJ(3));
t186 = t190 * t106;
t185 = Icges(4,5) + Icges(5,5);
t184 = Icges(4,6) + Icges(5,6);
t183 = t188 * t108 + t186;
t182 = t189 * t106 + t187;
t181 = rSges(5,1) + pkin(3);
t180 = Icges(2,4) + Icges(3,6);
t107 = sin(qJ(1));
t109 = cos(qJ(1));
t179 = t183 * t107 + t184 * t109;
t178 = t184 * t107 - t183 * t109;
t177 = t182 * t107 + t185 * t109;
t176 = t185 * t107 - t182 * t109;
t175 = -t188 * t106 + t187;
t174 = t189 * t108 - t186;
t173 = Icges(2,1) + Icges(3,2);
t172 = -Icges(3,4) + Icges(2,5);
t171 = Icges(3,5) - Icges(2,6);
t170 = Icges(2,2) + Icges(3,3);
t169 = Icges(4,3) + Icges(5,3);
t168 = t185 * t106 + t184 * t108;
t167 = rSges(5,3) + qJ(4);
t166 = t180 * t109;
t165 = rSges(5,2) * t108 + t181 * t106;
t164 = t180 * t107;
t100 = V_base(6) + qJD(1);
t97 = qJD(3) * t107 + V_base(5);
t98 = qJD(3) * t109 + V_base(4);
t163 = t100 * (t174 * t106 + t175 * t108) + (t177 * t106 + t179 * t108) * t98 + (t176 * t106 + t178 * t108) * t97;
t162 = -t170 * t109 - t164;
t161 = t170 * t107 - t166;
t160 = t173 * t107 + t166;
t159 = t173 * t109 - t164;
t156 = (t168 * t107 + t169 * t109) * t98 + (t169 * t107 - t168 * t109) * t97 + (-t184 * t106 + t185 * t108) * t100;
t148 = pkin(5) * t109;
t147 = t107 * pkin(5);
t145 = t165 * t107 + t167 * t109;
t144 = t167 * t107 - t165 * t109;
t88 = t107 * pkin(1) - qJ(2) * t109;
t136 = V_base(4) * t88 + V_base(3);
t135 = V_base(5) * pkin(4) + V_base(1);
t132 = -rSges(5,2) * t106 + t181 * t108;
t131 = -t88 - t147;
t130 = qJD(2) * t107 + t135;
t129 = V_base(5) * pkin(2) + t130;
t128 = rSges(4,1) * t106 + rSges(4,2) * t108;
t93 = pkin(1) * t109 + t107 * qJ(2);
t114 = -qJD(2) * t109 + t100 * t93 + V_base(2);
t111 = V_base(4) * t147 + (-t93 - t148) * V_base(5) + t136;
t110 = t100 * t148 + (-pkin(2) - pkin(4)) * V_base(4) + t114;
t95 = rSges(2,1) * t109 - t107 * rSges(2,2);
t94 = -rSges(3,2) * t109 + t107 * rSges(3,3);
t92 = rSges(4,1) * t108 - rSges(4,2) * t106;
t90 = t107 * rSges(2,1) + rSges(2,2) * t109;
t89 = -t107 * rSges(3,2) - rSges(3,3) * t109;
t69 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t68 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t67 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t62 = t107 * rSges(4,3) - t109 * t128;
t60 = rSges(4,3) * t109 + t107 * t128;
t46 = V_base(5) * rSges(2,3) - t100 * t90 + t135;
t45 = t100 * t95 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t44 = t90 * V_base(4) - t95 * V_base(5) + V_base(3);
t43 = V_base(5) * rSges(3,1) + (-t88 - t89) * t100 + t130;
t42 = t100 * t94 + (-rSges(3,1) - pkin(4)) * V_base(4) + t114;
t41 = t89 * V_base(4) + (-t93 - t94) * V_base(5) + t136;
t40 = t92 * t97 + (t131 - t62) * t100 + t129;
t39 = t100 * t60 - t98 * t92 + t110;
t38 = -t97 * t60 + t98 * t62 + t111;
t37 = qJD(4) * t109 + t132 * t97 + (t131 - t144) * t100 + t129;
t36 = qJD(4) * t107 + t100 * t145 - t132 * t98 + t110;
t35 = t144 * t98 - t145 * t97 + t111;
t1 = m(1) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(2) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(3) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(4) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(5) * (t35 ^ 2 + t36 ^ 2 + t37 ^ 2) / 0.2e1 + (t156 * t107 - t163 * t109) * t97 / 0.2e1 + (t163 * t107 + t156 * t109) * t98 / 0.2e1 + ((t107 * t162 + t160 * t109 + Icges(1,4)) * V_base(5) + (t161 * t107 + t159 * t109 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t160 * t107 - t162 * t109 + Icges(1,2)) * V_base(5) + (t107 * t159 - t109 * t161 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t179 * t106 + t177 * t108) * t98 + (-t178 * t106 + t176 * t108) * t97 + (-t175 * t106 + t174 * t108 + Icges(3,1) + Icges(2,3)) * t100) * t100 / 0.2e1 + t100 * V_base(5) * (t172 * t107 - t171 * t109) + t100 * V_base(4) * (t171 * t107 + t172 * t109) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
