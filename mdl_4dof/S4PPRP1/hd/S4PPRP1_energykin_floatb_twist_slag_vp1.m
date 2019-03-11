% Calculate kinetic energy for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:03
% EndTime: 2019-03-08 18:12:04
% DurationCPUTime: 0.73s
% Computational Cost: add. (337->135), mult. (554->159), div. (0->0), fcn. (478->4), ass. (0->68)
t157 = Icges(2,4) - Icges(3,5);
t156 = Icges(4,4) - Icges(5,5);
t155 = Icges(2,1) + Icges(3,1);
t154 = Icges(4,1) + Icges(5,1);
t153 = Icges(2,2) + Icges(3,3);
t152 = Icges(4,2) + Icges(5,3);
t111 = cos(pkin(5));
t126 = sin(pkin(5));
t132 = sin(qJ(3));
t133 = cos(qJ(3));
t77 = -t111 * t133 - t126 * t132;
t151 = t156 * t77;
t78 = t111 * t132 - t126 * t133;
t150 = t156 * t78;
t149 = t157 * t126;
t148 = t157 * t111;
t147 = -rSges(5,1) - pkin(3);
t146 = -t152 * t77 + t150;
t145 = t152 * t78 + t151;
t144 = t154 * t78 - t151;
t143 = t154 * t77 + t150;
t142 = -t111 * t153 - t149;
t141 = t126 * t153 - t148;
t140 = t126 * t155 + t148;
t139 = t111 * t155 - t149;
t138 = rSges(5,3) + qJ(4);
t137 = Icges(3,4) + Icges(2,5);
t136 = Icges(5,4) + Icges(4,5);
t135 = Icges(2,6) - Icges(3,6);
t134 = Icges(4,6) - Icges(5,6);
t131 = pkin(2) * t111;
t130 = -t138 * t77 + t147 * t78;
t129 = t138 * t78 + t147 * t77;
t124 = V_base(5) * qJ(1) + V_base(1);
t121 = qJD(1) + V_base(3);
t120 = t126 * pkin(2);
t97 = t111 * pkin(1) + t126 * qJ(2);
t119 = -t97 - t131;
t94 = t126 * pkin(1) - t111 * qJ(2);
t117 = V_base(4) * t94 + t121;
t116 = qJD(2) * t126 + t124;
t115 = V_base(4) * t120 + t117;
t114 = -qJD(2) * t111 + V_base(6) * t97 + V_base(2);
t113 = V_base(4) * pkin(4) + V_base(6) * t131 + t114;
t112 = (-t120 - t94) * V_base(6) + t116;
t108 = V_base(6) - qJD(3);
t99 = t111 * rSges(2,1) - t126 * rSges(2,2);
t98 = t111 * rSges(3,1) + t126 * rSges(3,3);
t96 = t126 * rSges(2,1) + t111 * rSges(2,2);
t95 = t126 * rSges(3,1) - t111 * rSges(3,3);
t81 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t80 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t79 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t72 = V_base(5) * rSges(2,3) - t96 * V_base(6) + t124;
t71 = t99 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t70 = -rSges(4,1) * t77 - rSges(4,2) * t78;
t67 = -rSges(4,1) * t78 + rSges(4,2) * t77;
t64 = t96 * V_base(4) - t99 * V_base(5) + t121;
t51 = V_base(5) * rSges(3,2) + (-t94 - t95) * V_base(6) + t116;
t50 = V_base(6) * t98 + (-rSges(3,2) - qJ(1)) * V_base(4) + t114;
t49 = t95 * V_base(4) + (-t97 - t98) * V_base(5) + t117;
t48 = -t108 * t67 + (-pkin(4) - rSges(4,3)) * V_base(5) + t112;
t47 = t108 * t70 + (rSges(4,3) - qJ(1)) * V_base(4) + t113;
t46 = V_base(4) * t67 + (t119 - t70) * V_base(5) + t115;
t45 = qJD(4) * t78 + (-pkin(4) - rSges(5,2)) * V_base(5) - t130 * t108 + t112;
t44 = -qJD(4) * t77 + (rSges(5,2) - qJ(1)) * V_base(4) + t129 * t108 + t113;
t43 = t130 * V_base(4) + (t119 - t129) * V_base(5) + t115;
t1 = m(1) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t64 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(3) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(4) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(5) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + ((Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t108 + (t134 * V_base(4) + t136 * V_base(5)) * t78 + (-t134 * V_base(5) + t136 * V_base(4)) * t77) * t108 + ((Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (-t135 * V_base(4) + t137 * V_base(5)) * t126 + (t135 * V_base(5) + t137 * V_base(4)) * t111) * V_base(6) + ((t140 * t111 + t142 * t126 + t144 * t77 + t146 * t78) * V_base(5) + (t139 * t111 + t141 * t126 + t143 * t77 + t145 * t78) * V_base(4)) * V_base(4) / 0.2e1 + ((-t142 * t111 + t140 * t126 + t144 * t78 - t146 * t77) * V_base(5) + (-t141 * t111 + t139 * t126 + t143 * t78 - t145 * t77) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
