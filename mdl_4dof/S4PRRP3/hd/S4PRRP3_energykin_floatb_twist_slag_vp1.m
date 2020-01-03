% Calculate kinetic energy for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:41
% EndTime: 2019-12-31 16:26:43
% DurationCPUTime: 1.25s
% Computational Cost: add. (643->159), mult. (642->206), div. (0->0), fcn. (472->6), ass. (0->83)
t181 = Icges(4,4) + Icges(5,4);
t180 = Icges(4,1) + Icges(5,1);
t179 = Icges(4,2) + Icges(5,2);
t118 = cos(qJ(3));
t178 = t181 * t118;
t117 = sin(qJ(3));
t177 = t181 * t117;
t176 = Icges(4,5) + Icges(5,5);
t175 = Icges(4,6) + Icges(5,6);
t174 = -t179 * t117 + t178;
t173 = t180 * t118 - t177;
t172 = rSges(5,1) + pkin(3);
t113 = pkin(6) + qJ(2);
t107 = sin(t113);
t108 = cos(t113);
t171 = t174 * t107 - t175 * t108;
t170 = t175 * t107 + t174 * t108;
t169 = t173 * t107 - t176 * t108;
t168 = t176 * t107 + t173 * t108;
t167 = t179 * t118 + t177;
t166 = t180 * t117 + t178;
t165 = Icges(4,3) + Icges(5,3);
t164 = -t175 * t117 + t176 * t118;
t163 = rSges(5,3) + qJ(4);
t162 = -rSges(5,2) * t117 + t172 * t118;
t110 = V_base(6) + qJD(2);
t91 = -qJD(3) * t108 + V_base(5);
t92 = qJD(3) * t107 + V_base(4);
t159 = (-t170 * t117 + t168 * t118) * t92 + (-t171 * t117 + t169 * t118) * t91 + (-t167 * t117 + t166 * t118) * t110;
t158 = (t165 * t107 + t164 * t108) * t92 + (t164 * t107 - t165 * t108) * t91 + (t176 * t117 + t175 * t118) * t110;
t115 = cos(pkin(6));
t154 = pkin(1) * t115;
t152 = -pkin(4) - qJ(1);
t150 = t162 * t107 - t163 * t108;
t149 = t163 * t107 + t162 * t108;
t114 = sin(pkin(6));
t148 = Icges(2,4) * t114;
t147 = Icges(3,4) * t107;
t136 = pkin(1) * V_base(6);
t142 = t115 * t136 + V_base(2);
t141 = V_base(5) * qJ(1) + V_base(1);
t137 = qJD(1) + V_base(3);
t135 = rSges(5,2) * t118 + t172 * t117;
t134 = V_base(4) * t114 * pkin(1) + t137;
t133 = rSges(4,1) * t118 - rSges(4,2) * t117;
t123 = V_base(5) * pkin(4) - t114 * t136 + t141;
t81 = pkin(2) * t108 + pkin(5) * t107;
t122 = t110 * t81 + t152 * V_base(4) + t142;
t80 = pkin(2) * t107 - pkin(5) * t108;
t121 = V_base(4) * t80 + (-t81 - t154) * V_base(5) + t134;
t109 = Icges(2,4) * t115;
t105 = Icges(3,4) * t108;
t102 = t117 * rSges(4,1) + rSges(4,2) * t118;
t94 = rSges(2,1) * t115 - rSges(2,2) * t114;
t93 = rSges(2,1) * t114 + rSges(2,2) * t115;
t90 = Icges(2,1) * t115 - t148;
t89 = Icges(2,1) * t114 + t109;
t88 = -Icges(2,2) * t114 + t109;
t87 = Icges(2,2) * t115 + t148;
t84 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t83 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t82 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t79 = rSges(3,1) * t108 - rSges(3,2) * t107;
t78 = rSges(3,1) * t107 + rSges(3,2) * t108;
t77 = Icges(3,1) * t108 - t147;
t76 = Icges(3,1) * t107 + t105;
t75 = -Icges(3,2) * t107 + t105;
t74 = Icges(3,2) * t108 + t147;
t69 = V_base(5) * rSges(2,3) - t93 * V_base(6) + t141;
t68 = t94 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t67 = t107 * rSges(4,3) + t133 * t108;
t65 = -t108 * rSges(4,3) + t133 * t107;
t51 = t93 * V_base(4) - t94 * V_base(5) + t137;
t48 = V_base(5) * rSges(3,3) - t110 * t78 + t123;
t47 = t110 * t79 + (-rSges(3,3) + t152) * V_base(4) + t142;
t46 = t78 * V_base(4) + (-t79 - t154) * V_base(5) + t134;
t45 = t102 * t91 + (-t65 - t80) * t110 + t123;
t44 = -t102 * t92 + t110 * t67 + t122;
t43 = t65 * t92 - t67 * t91 + t121;
t42 = qJD(4) * t107 + t135 * t91 + (-t80 - t150) * t110 + t123;
t41 = -qJD(4) * t108 + t149 * t110 - t135 * t92 + t122;
t40 = -t149 * t91 + t150 * t92 + t121;
t1 = m(1) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(2) * (t51 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + (t107 * t159 - t158 * t108) * t91 / 0.2e1 + (t107 * t158 + t159 * t108) * t92 / 0.2e1 + ((t168 * t117 + t170 * t118) * t92 + (t169 * t117 + t171 * t118) * t91 + (t166 * t117 + t167 * t118 + Icges(3,3)) * t110) * t110 / 0.2e1 + ((-t107 * t74 + t108 * t76 - t114 * t87 + t115 * t89 + Icges(1,4)) * V_base(5) + (-t107 * t75 + t108 * t77 - t114 * t88 + t115 * t90 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t107 * t76 + t108 * t74 + t114 * t89 + t115 * t87 + Icges(1,2)) * V_base(5) + (t107 * t77 + t108 * t75 + t114 * t90 + t115 * t88 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t110 * (Icges(3,5) * t108 - Icges(3,6) * t107) + V_base(5) * t110 * (Icges(3,5) * t107 + Icges(3,6) * t108) + ((Icges(2,5) * t114 + Icges(2,6) * t115 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t115 - Icges(2,6) * t114 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
