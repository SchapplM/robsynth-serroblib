% Calculate kinetic energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:29
% EndTime: 2019-12-31 16:42:30
% DurationCPUTime: 1.18s
% Computational Cost: add. (654->160), mult. (642->203), div. (0->0), fcn. (472->6), ass. (0->83)
t180 = Icges(4,4) + Icges(5,4);
t179 = Icges(4,1) + Icges(5,1);
t178 = Icges(4,2) + Icges(5,2);
t117 = cos(qJ(3));
t177 = t180 * t117;
t115 = sin(qJ(3));
t176 = t180 * t115;
t175 = Icges(4,5) + Icges(5,5);
t174 = Icges(4,6) + Icges(5,6);
t173 = -t115 * t178 + t177;
t172 = t117 * t179 - t176;
t171 = rSges(5,1) + pkin(3);
t113 = qJ(1) + pkin(6);
t107 = sin(t113);
t108 = cos(t113);
t170 = t173 * t107 - t174 * t108;
t169 = t174 * t107 + t173 * t108;
t168 = t172 * t107 - t175 * t108;
t167 = t175 * t107 + t172 * t108;
t166 = t117 * t178 + t176;
t165 = t115 * t179 + t177;
t164 = Icges(4,3) + Icges(5,3);
t163 = -t174 * t115 + t175 * t117;
t162 = rSges(5,3) + qJ(4);
t161 = -rSges(5,2) * t115 + t171 * t117;
t109 = V_base(6) + qJD(1);
t85 = -qJD(3) * t108 + V_base(5);
t86 = qJD(3) * t107 + V_base(4);
t158 = (-t169 * t115 + t167 * t117) * t86 + (-t170 * t115 + t168 * t117) * t85 + (-t166 * t115 + t165 * t117) * t109;
t157 = (t164 * t107 + t163 * t108) * t86 + (t163 * t107 - t164 * t108) * t85 + (t175 * t115 + t174 * t117) * t109;
t116 = sin(qJ(1));
t153 = pkin(1) * t116;
t118 = cos(qJ(1));
t152 = pkin(1) * t118;
t150 = -pkin(4) - qJ(2);
t148 = t161 * t107 - t162 * t108;
t147 = t162 * t107 + t161 * t108;
t146 = Icges(2,4) * t116;
t145 = Icges(3,4) * t107;
t140 = t109 * t152 + V_base(2);
t139 = V_base(5) * pkin(4) + V_base(1);
t80 = pkin(2) * t107 - pkin(5) * t108;
t136 = -t80 - t153;
t135 = rSges(5,2) * t117 + t171 * t115;
t134 = V_base(5) * qJ(2) + t139;
t133 = V_base(4) * t153 + qJD(2) + V_base(3);
t132 = rSges(4,1) * t117 - rSges(4,2) * t115;
t81 = pkin(2) * t108 + pkin(5) * t107;
t122 = t109 * t81 + t150 * V_base(4) + t140;
t121 = V_base(4) * t80 + (-t81 - t152) * V_base(5) + t133;
t111 = Icges(2,4) * t118;
t105 = Icges(3,4) * t108;
t102 = rSges(2,1) * t118 - t116 * rSges(2,2);
t101 = t116 * rSges(2,1) + rSges(2,2) * t118;
t100 = rSges(4,1) * t115 + rSges(4,2) * t117;
t98 = Icges(2,1) * t118 - t146;
t97 = Icges(2,1) * t116 + t111;
t94 = -Icges(2,2) * t116 + t111;
t93 = Icges(2,2) * t118 + t146;
t84 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t83 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t82 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t79 = rSges(3,1) * t108 - rSges(3,2) * t107;
t78 = rSges(3,1) * t107 + rSges(3,2) * t108;
t77 = Icges(3,1) * t108 - t145;
t76 = Icges(3,1) * t107 + t105;
t75 = -Icges(3,2) * t107 + t105;
t74 = Icges(3,2) * t108 + t145;
t69 = rSges(4,3) * t107 + t108 * t132;
t67 = -rSges(4,3) * t108 + t107 * t132;
t65 = V_base(5) * rSges(2,3) - t101 * t109 + t139;
t64 = t102 * t109 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t51 = t101 * V_base(4) - t102 * V_base(5) + V_base(3);
t48 = V_base(5) * rSges(3,3) + (-t78 - t153) * t109 + t134;
t47 = t109 * t79 + (-rSges(3,3) + t150) * V_base(4) + t140;
t46 = V_base(4) * t78 + (-t79 - t152) * V_base(5) + t133;
t45 = t100 * t85 + (t136 - t67) * t109 + t134;
t44 = -t100 * t86 + t109 * t69 + t122;
t43 = t86 * t67 - t85 * t69 + t121;
t42 = qJD(4) * t107 + t135 * t85 + (t136 - t148) * t109 + t134;
t41 = -qJD(4) * t108 + t109 * t147 - t135 * t86 + t122;
t40 = -t147 * t85 + t148 * t86 + t121;
t1 = m(1) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + (t158 * t107 - t157 * t108) * t85 / 0.2e1 + (t157 * t107 + t158 * t108) * t86 / 0.2e1 + ((-t107 * t74 + t108 * t76 - t116 * t93 + t118 * t97 + Icges(1,4)) * V_base(5) + (-t107 * t75 + t108 * t77 - t116 * t94 + t118 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t107 * t76 + t108 * t74 + t116 * t97 + t118 * t93 + Icges(1,2)) * V_base(5) + (t107 * t77 + t108 * t75 + t116 * t98 + t118 * t94 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t167 * t115 + t169 * t117) * t86 + (t168 * t115 + t170 * t117) * t85 + (t165 * t115 + t166 * t117 + Icges(2,3) + Icges(3,3)) * t109) * t109 / 0.2e1 + t109 * V_base(5) * (Icges(2,5) * t116 + Icges(3,5) * t107 + Icges(2,6) * t118 + Icges(3,6) * t108) + t109 * V_base(4) * (Icges(2,5) * t118 + Icges(3,5) * t108 - Icges(2,6) * t116 - Icges(3,6) * t107) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
