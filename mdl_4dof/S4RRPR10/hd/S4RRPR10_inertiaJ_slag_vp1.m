% Calculate joint inertia matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:13
% DurationCPUTime: 1.19s
% Computational Cost: add. (1029->202), mult. (2485->314), div. (0->0), fcn. (2585->6), ass. (0->110)
t158 = -Icges(4,4) + Icges(3,5);
t157 = Icges(4,5) - Icges(3,6);
t156 = Icges(4,1) + Icges(3,3);
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t155 = t157 * t92 + t158 * t95;
t154 = (Icges(3,6) / 0.2e1 - Icges(4,5) / 0.2e1) * t95 + (Icges(3,5) / 0.2e1 - Icges(4,4) / 0.2e1) * t92;
t93 = sin(qJ(1));
t153 = -t93 / 0.2e1;
t141 = t93 / 0.2e1;
t96 = cos(qJ(1));
t152 = -t96 / 0.2e1;
t151 = t96 / 0.2e1;
t146 = -t155 * t93 + t156 * t96;
t145 = t155 * t96 + t156 * t93;
t135 = t95 * t96;
t94 = cos(qJ(4));
t131 = t96 * t94;
t91 = sin(qJ(4));
t138 = t93 * t91;
t60 = t92 * t131 - t138;
t132 = t96 * t91;
t137 = t93 * t94;
t61 = t92 * t132 + t137;
t30 = t61 * rSges(5,1) + t60 * rSges(5,2) + rSges(5,3) * t135;
t144 = t93 * pkin(3) + pkin(6) * t135 + t30;
t89 = t93 ^ 2;
t90 = t96 ^ 2;
t143 = m(4) / 0.2e1;
t142 = m(5) / 0.2e1;
t140 = t96 * pkin(3);
t139 = t92 * t96;
t136 = t93 * t95;
t134 = t96 * rSges(4,1);
t133 = t96 * rSges(3,3);
t117 = qJ(3) * t92;
t128 = pkin(2) * t135 + t96 * t117;
t130 = t89 * (pkin(2) * t95 + t117) + t96 * t128;
t71 = t92 * pkin(2) - t95 * qJ(3);
t129 = t92 * rSges(4,2) + t95 * rSges(4,3) - t71;
t127 = t96 * pkin(1) + t93 * pkin(5);
t125 = t89 + t90;
t124 = Icges(3,4) * t92;
t123 = Icges(3,4) * t95;
t122 = Icges(5,5) * t95;
t121 = Icges(4,6) * t92;
t120 = Icges(4,6) * t95;
t119 = Icges(5,6) * t95;
t118 = Icges(5,3) * t95;
t24 = Icges(5,5) * t61 + Icges(5,6) * t60 + t96 * t118;
t26 = Icges(5,4) * t61 + Icges(5,2) * t60 + t96 * t119;
t28 = Icges(5,1) * t61 + Icges(5,4) * t60 + t96 * t122;
t10 = t92 * t24 + (-t26 * t94 - t28 * t91) * t95;
t39 = Icges(5,3) * t92 + (-Icges(5,5) * t91 - Icges(5,6) * t94) * t95;
t42 = Icges(5,6) * t92 + (-Icges(5,4) * t91 - Icges(5,2) * t94) * t95;
t45 = Icges(5,5) * t92 + (-Icges(5,1) * t91 - Icges(5,4) * t94) * t95;
t12 = t39 * t135 + t60 * t42 + t61 * t45;
t116 = t10 / 0.2e1 + t12 / 0.2e1;
t62 = t92 * t137 + t132;
t63 = t92 * t138 - t131;
t25 = Icges(5,5) * t63 + Icges(5,6) * t62 + t93 * t118;
t27 = Icges(5,4) * t63 + Icges(5,2) * t62 + t93 * t119;
t29 = Icges(5,1) * t63 + Icges(5,4) * t62 + t93 * t122;
t11 = t92 * t25 + (-t27 * t94 - t29 * t91) * t95;
t13 = t39 * t136 + t62 * t42 + t63 * t45;
t115 = t11 / 0.2e1 + t13 / 0.2e1;
t114 = t127 + t128;
t54 = t92 * rSges(5,3) + (-rSges(5,1) * t91 - rSges(5,2) * t94) * t95;
t113 = -pkin(6) * t92 - t54 - t71;
t111 = rSges(3,1) * t95 - rSges(3,2) * t92;
t110 = -t63 * rSges(5,1) - t62 * rSges(5,2);
t109 = -t94 * t42 - t91 * t45;
t104 = Icges(3,1) * t95 - t124;
t103 = -Icges(3,2) * t92 + t123;
t100 = -Icges(4,2) * t95 + t121;
t99 = Icges(4,3) * t92 - t120;
t98 = rSges(3,1) * t135 - rSges(3,2) * t139 + t93 * rSges(3,3);
t97 = t93 * rSges(4,1) - rSges(4,2) * t135 + rSges(4,3) * t139;
t86 = t96 * pkin(5);
t75 = t96 * rSges(2,1) - t93 * rSges(2,2);
t74 = -t93 * rSges(2,1) - t96 * rSges(2,2);
t73 = t92 * rSges(3,1) + t95 * rSges(3,2);
t38 = t92 * t39;
t37 = t129 * t96;
t36 = t129 * t93;
t35 = t98 + t127;
t34 = t133 + t86 + (-pkin(1) - t111) * t93;
t33 = t97 + t114;
t32 = t134 + t86 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t95 + (-rSges(4,3) - qJ(3)) * t92) * t93;
t31 = rSges(5,3) * t136 - t110;
t23 = t113 * t96;
t22 = t113 * t93;
t21 = t96 * t98 + (t111 * t93 - t133) * t93;
t20 = -t54 * t135 + t92 * t30;
t19 = t54 * t136 - t92 * t31;
t18 = t114 + t144;
t17 = t140 + t86 + (-t117 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(6)) * t95) * t93 + t110;
t16 = (t109 * t95 + t38) * t92;
t15 = t96 * t97 + (-t134 + (-rSges(4,2) * t95 + rSges(4,3) * t92) * t93) * t93 + t130;
t14 = (-t30 * t93 + t31 * t96) * t95;
t9 = t144 * t96 + (pkin(6) * t136 - t140 + t31) * t93 + t130;
t8 = t25 * t136 + t62 * t27 + t63 * t29;
t7 = t24 * t136 + t62 * t26 + t63 * t28;
t6 = t25 * t135 + t60 * t27 + t61 * t29;
t5 = t24 * t135 + t60 * t26 + t61 * t28;
t4 = t7 * t93 - t8 * t96;
t3 = t5 * t93 - t6 * t96;
t2 = t13 * t92 + (t7 * t96 + t8 * t93) * t95;
t1 = t12 * t92 + (t5 * t96 + t6 * t93) * t95;
t40 = [Icges(2,3) + t38 + m(2) * (t74 ^ 2 + t75 ^ 2) + m(3) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + (t109 + t121 + t124 + (Icges(3,2) + Icges(4,3)) * t95) * t95 + (t120 + t123 + (Icges(3,1) + Icges(4,2)) * t92) * t92; ((t103 * t153 + t99 * t141) * t95 + (t100 * t141 + t104 * t153) * t92 - t115 + t154 * t96) * t96 + ((t103 * t151 + t99 * t152) * t95 + (t100 * t152 + t104 * t151) * t92 + t116 + t154 * t93) * t93 + m(4) * (t37 * t32 + t36 * t33) + m(5) * (t23 * t17 + t22 * t18) + m(3) * (-t34 * t96 - t35 * t93) * t73 + (-t157 * t95 + t158 * t92) * (t89 / 0.2e1 + t90 / 0.2e1); m(5) * (t22 ^ 2 + t23 ^ 2 + t9 ^ 2) + m(4) * (t15 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(3) * (t125 * t73 ^ 2 + t21 ^ 2) + (t145 * t89 + t3) * t93 + (-t4 + t146 * t90 + (t145 * t96 + t146 * t93) * t93) * t96; 0.2e1 * ((t32 * t96 + t33 * t93) * t143 + (t17 * t96 + t18 * t93) * t142) * t92; m(5) * (-t95 * t9 + (t22 * t93 + t23 * t96) * t92) + m(4) * (-t95 * t15 + (t36 * t93 + t37 * t96) * t92); 0.2e1 * (t143 + t142) * (t125 * t92 ^ 2 + t95 ^ 2); t16 + m(5) * (t19 * t17 + t20 * t18) + (t115 * t93 + t116 * t96) * t95; m(5) * (t14 * t9 + t19 * t23 + t20 * t22) + t1 * t141 + t92 * (t10 * t93 - t11 * t96) / 0.2e1 + t2 * t152 + (t4 * t141 + t3 * t151) * t95; m(5) * (-t14 * t95 + (t19 * t96 + t20 * t93) * t92); m(5) * (t14 ^ 2 + t19 ^ 2 + t20 ^ 2) + t92 * t16 + (t96 * t1 + t93 * t2 + t92 * (t10 * t96 + t11 * t93)) * t95;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t40(1), t40(2), t40(4), t40(7); t40(2), t40(3), t40(5), t40(8); t40(4), t40(5), t40(6), t40(9); t40(7), t40(8), t40(9), t40(10);];
Mq = res;
