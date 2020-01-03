% Calculate kinetic energy for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:38
% EndTime: 2019-12-31 16:19:39
% DurationCPUTime: 1.42s
% Computational Cost: add. (446->201), mult. (870->284), div. (0->0), fcn. (758->6), ass. (0->97)
t178 = Icges(2,4) + Icges(3,6);
t177 = Icges(2,1) + Icges(3,2);
t176 = -Icges(3,4) + Icges(2,5);
t175 = Icges(3,5) - Icges(2,6);
t174 = Icges(2,2) + Icges(3,3);
t127 = cos(pkin(6));
t173 = t178 * t127;
t126 = sin(pkin(6));
t172 = t178 * t126;
t171 = t174 * t127 + t172;
t170 = -t174 * t126 + t173;
t169 = t177 * t126 + t173;
t168 = t177 * t127 - t172;
t129 = sin(qJ(3));
t131 = cos(qJ(3));
t160 = Icges(4,4) * t131;
t111 = -Icges(4,2) * t129 + t160;
t161 = Icges(4,4) * t129;
t112 = Icges(4,1) * t131 - t161;
t115 = qJD(3) * t126 + V_base(5);
t116 = qJD(3) * t127 + V_base(4);
t138 = Icges(4,2) * t131 + t161;
t69 = Icges(4,6) * t127 + t138 * t126;
t70 = Icges(4,6) * t126 - t138 * t127;
t139 = Icges(4,1) * t129 + t160;
t71 = Icges(4,5) * t127 + t139 * t126;
t72 = Icges(4,5) * t126 - t139 * t127;
t165 = (t129 * t71 + t131 * t69) * t116 + (t129 * t72 + t131 * t70) * t115 + (t111 * t131 + t112 * t129) * V_base(6);
t164 = pkin(4) * t126;
t163 = pkin(4) * t127;
t157 = t126 * t131;
t156 = t127 * t131;
t128 = sin(qJ(4));
t155 = t128 * t129;
t130 = cos(qJ(4));
t154 = t129 * t130;
t153 = qJD(4) * t131;
t152 = V_base(5) * qJ(1) + V_base(1);
t148 = qJD(1) + V_base(3);
t104 = pkin(1) * t126 - qJ(2) * t127;
t147 = -t104 - t164;
t146 = V_base(4) * t104 + t148;
t145 = qJD(2) * t126 + t152;
t144 = V_base(5) * pkin(2) + t145;
t143 = pkin(3) * t129 - pkin(5) * t131;
t142 = rSges(4,1) * t129 + rSges(4,2) * t131;
t137 = Icges(4,5) * t129 + Icges(4,6) * t131;
t107 = pkin(1) * t127 + qJ(2) * t126;
t135 = -qJD(2) * t127 + V_base(6) * t107 + V_base(2);
t134 = (Icges(4,5) * t131 - Icges(4,6) * t129) * V_base(6) + t115 * (Icges(4,3) * t126 - t137 * t127) + t116 * (Icges(4,3) * t127 + t137 * t126);
t133 = V_base(4) * t164 + (-t107 - t163) * V_base(5) + t146;
t132 = V_base(6) * t163 + (-pkin(2) - qJ(1)) * V_base(4) + t135;
t117 = qJD(4) * t129 + V_base(6);
t114 = pkin(3) * t131 + t129 * pkin(5);
t113 = rSges(4,1) * t131 - t129 * rSges(4,2);
t109 = rSges(2,1) * t127 - rSges(2,2) * t126;
t108 = -rSges(3,2) * t127 + rSges(3,3) * t126;
t106 = rSges(2,1) * t126 + rSges(2,2) * t127;
t105 = -rSges(3,2) * t126 - rSges(3,3) * t127;
t91 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t90 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t89 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t86 = t143 * t127;
t85 = t143 * t126;
t84 = t126 * t128 - t127 * t154;
t83 = t126 * t130 + t127 * t155;
t82 = t126 * t154 + t127 * t128;
t81 = -t126 * t155 + t127 * t130;
t80 = -t126 * t153 + t116;
t79 = t127 * t153 + t115;
t78 = t129 * rSges(5,3) + (rSges(5,1) * t130 - rSges(5,2) * t128) * t131;
t77 = Icges(5,5) * t129 + (Icges(5,1) * t130 - Icges(5,4) * t128) * t131;
t76 = Icges(5,6) * t129 + (Icges(5,4) * t130 - Icges(5,2) * t128) * t131;
t75 = Icges(5,3) * t129 + (Icges(5,5) * t130 - Icges(5,6) * t128) * t131;
t74 = t126 * rSges(4,3) - t142 * t127;
t73 = t127 * rSges(4,3) + t142 * t126;
t66 = V_base(5) * rSges(2,3) - t106 * V_base(6) + t152;
t65 = t109 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t64 = t106 * V_base(4) - t109 * V_base(5) + t148;
t63 = t84 * rSges(5,1) + t83 * rSges(5,2) + rSges(5,3) * t156;
t62 = t82 * rSges(5,1) + t81 * rSges(5,2) - rSges(5,3) * t157;
t61 = V_base(5) * rSges(3,1) + (-t104 - t105) * V_base(6) + t145;
t60 = t108 * V_base(6) + (-rSges(3,1) - qJ(1)) * V_base(4) + t135;
t59 = Icges(5,1) * t84 + Icges(5,4) * t83 + Icges(5,5) * t156;
t58 = Icges(5,1) * t82 + Icges(5,4) * t81 - Icges(5,5) * t157;
t57 = Icges(5,4) * t84 + Icges(5,2) * t83 + Icges(5,6) * t156;
t56 = Icges(5,4) * t82 + Icges(5,2) * t81 - Icges(5,6) * t157;
t55 = Icges(5,5) * t84 + Icges(5,6) * t83 + Icges(5,3) * t156;
t54 = Icges(5,5) * t82 + Icges(5,6) * t81 - Icges(5,3) * t157;
t53 = t105 * V_base(4) + (-t107 - t108) * V_base(5) + t146;
t52 = t113 * t115 + (t147 - t74) * V_base(6) + t144;
t51 = -t113 * t116 + t73 * V_base(6) + t132;
t50 = -t115 * t73 + t116 * t74 + t133;
t49 = t114 * t115 - t117 * t63 + t78 * t79 + (t147 + t86) * V_base(6) + t144;
t48 = -t114 * t116 + t117 * t62 - t78 * t80 + t85 * V_base(6) + t132;
t47 = -t115 * t85 - t116 * t86 - t62 * t79 + t63 * t80 + t133;
t1 = m(1) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(2) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + t116 * (t165 * t126 + t134 * t127) / 0.2e1 + t115 * (t134 * t126 - t165 * t127) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t80 * ((-t54 * t157 + t81 * t56 + t82 * t58) * t80 + (-t55 * t157 + t81 * t57 + t82 * t59) * t79 + (-t75 * t157 + t81 * t76 + t82 * t77) * t117) / 0.2e1 + t79 * ((t54 * t156 + t83 * t56 + t84 * t58) * t80 + (t55 * t156 + t83 * t57 + t84 * t59) * t79 + (t75 * t156 + t83 * t76 + t84 * t77) * t117) / 0.2e1 + t117 * ((t117 * t75 + t54 * t80 + t55 * t79) * t129 + ((-t128 * t56 + t130 * t58) * t80 + (-t128 * t57 + t130 * t59) * t79 + (-t128 * t76 + t130 * t77) * t117) * t131) / 0.2e1 + ((-t171 * t126 + t169 * t127 + Icges(1,4)) * V_base(5) + (-t170 * t126 + t168 * t127 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t169 * t126 + t171 * t127 + Icges(1,2)) * V_base(5) + (t168 * t126 + t170 * t127 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t129 * t69 + t131 * t71) * t116 + (-t129 * t70 + t131 * t72) * t115 + (-t129 * t111 + t131 * t112 + Icges(3,1) + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (t175 * t126 + t176 * t127 + Icges(1,5)) + V_base(6) * V_base(5) * (t176 * t126 - t175 * t127 + Icges(1,6));
T = t1;
