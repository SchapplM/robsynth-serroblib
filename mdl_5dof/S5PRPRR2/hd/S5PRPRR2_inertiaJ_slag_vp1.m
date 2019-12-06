% Calculate joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:46
% EndTime: 2019-12-05 15:44:50
% DurationCPUTime: 0.96s
% Computational Cost: add. (3550->165), mult. (3458->273), div. (0->0), fcn. (3650->10), ass. (0->91)
t99 = sin(pkin(8));
t96 = t99 ^ 2;
t100 = cos(pkin(8));
t97 = t100 ^ 2;
t154 = t96 + t97;
t153 = Icges(3,3) + Icges(4,3);
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t98 = qJ(2) + pkin(9);
t91 = sin(t98);
t92 = cos(t98);
t152 = Icges(3,5) * t105 + Icges(4,5) * t92 - Icges(3,6) * t103 - Icges(4,6) * t91;
t93 = qJ(4) + t98;
t88 = sin(t93);
t89 = cos(t93);
t116 = Icges(5,4) * t89 - Icges(5,2) * t88;
t118 = Icges(5,1) * t89 - Icges(5,4) * t88;
t123 = -(Icges(5,6) * t99 + t116 * t100) * t88 + (Icges(5,5) * t99 + t118 * t100) * t89;
t124 = (-Icges(5,6) * t100 + t116 * t99) * t88 - (-Icges(5,5) * t100 + t118 * t99) * t89;
t114 = Icges(5,5) * t89 - Icges(5,6) * t88;
t54 = -Icges(5,3) * t100 + t114 * t99;
t55 = Icges(5,3) * t99 + t114 * t100;
t144 = t88 * t99;
t104 = cos(qJ(5));
t132 = t100 * t104;
t102 = sin(qJ(5));
t136 = t99 * t102;
t76 = -t89 * t136 - t132;
t133 = t100 * t102;
t135 = t99 * t104;
t77 = t89 * t135 - t133;
t31 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t144;
t33 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t144;
t35 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t144;
t14 = t31 * t144 + t76 * t33 + t77 * t35;
t137 = t100 * t88;
t78 = -t89 * t133 + t135;
t79 = t89 * t132 + t136;
t32 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t137;
t34 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t137;
t36 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t137;
t15 = t32 * t144 + t76 * t34 + t77 * t36;
t8 = -t14 * t100 + t15 * t99;
t150 = -t97 * t54 - (t123 * t99 + (t124 - t55) * t100) * t99 - t8;
t149 = t153 * t100 - t152 * t99;
t148 = t152 * t100 + t153 * t99;
t16 = t31 * t137 + t78 * t33 + t79 * t35;
t17 = t32 * t137 + t78 * t34 + t79 * t36;
t9 = -t16 * t100 + t17 * t99;
t146 = (t96 * t55 + t9 + (t124 * t100 + (t123 - t54) * t99) * t100) * t99;
t145 = pkin(2) * t103;
t48 = -Icges(6,3) * t89 + (Icges(6,5) * t104 - Icges(6,6) * t102) * t88;
t143 = t89 * t48;
t28 = t154 * (rSges(5,1) * t89 - rSges(5,2) * t88);
t51 = -t89 * rSges(6,3) + (rSges(6,1) * t104 - rSges(6,2) * t102) * t88;
t141 = -t88 * pkin(4) + t89 * pkin(7) - t51;
t140 = t154 * t105 * pkin(2);
t131 = -t91 * rSges(4,1) - t92 * rSges(4,2) - t145;
t130 = pkin(3) * t154 * t92 + t140;
t37 = t77 * rSges(6,1) + t76 * rSges(6,2) + rSges(6,3) * t144;
t38 = t79 * rSges(6,1) + t78 * rSges(6,2) + rSges(6,3) * t137;
t21 = t100 * t38 + t99 * t37 + t154 * (pkin(4) * t89 + pkin(7) * t88);
t19 = -t89 * t31 + (-t102 * t33 + t104 * t35) * t88;
t20 = -t89 * t32 + (-t102 * t34 + t104 * t36) * t88;
t49 = -Icges(6,6) * t89 + (Icges(6,4) * t104 - Icges(6,2) * t102) * t88;
t50 = -Icges(6,5) * t89 + (Icges(6,1) * t104 - Icges(6,4) * t102) * t88;
t3 = -(t76 * t49 + t77 * t50) * t89 + (t15 * t100 + (t14 - t143) * t99) * t88;
t4 = -(t78 * t49 + t79 * t50) * t89 + (t16 * t99 + (t17 - t143) * t100) * t88;
t129 = t99 * t4 / 0.2e1 - t89 * (-t19 * t100 + t20 * t99) / 0.2e1 - t100 * t3 / 0.2e1 + t8 * t144 / 0.2e1 + t9 * t137 / 0.2e1;
t127 = -pkin(3) * t91 - t145;
t108 = t150 * t100 + t146;
t81 = t88 * rSges(5,1) + t89 * rSges(5,2);
t107 = t127 - t81;
t106 = t127 + t141;
t86 = t103 * rSges(3,1) + t105 * rSges(3,2);
t63 = t131 * t100;
t62 = t131 * t99;
t45 = t107 * t100;
t44 = t107 * t99;
t41 = t154 * (rSges(3,1) * t105 - rSges(3,2) * t103);
t40 = t141 * t100;
t39 = t141 * t99;
t27 = t106 * t100;
t26 = t106 * t99;
t25 = -t51 * t137 - t89 * t38;
t24 = t51 * t144 + t89 * t37;
t23 = t140 + t154 * (rSges(4,1) * t92 - rSges(4,2) * t91);
t22 = (t100 * t37 - t38 * t99) * t88;
t18 = t130 + t28;
t11 = t21 + t130;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t41 + m(4) * t23 + m(5) * t18 + m(6) * t11; m(6) * (t11 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t18 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t23 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(3) * (t154 * t86 ^ 2 + t41 ^ 2) + t146 + t148 * t99 * t96 + (t149 * t97 + (t148 * t100 + t149 * t99) * t99 + t150) * t100; 0; m(6) * (-t100 * t26 + t99 * t27) + m(5) * (-t100 * t44 + t99 * t45) + m(4) * (-t100 * t62 + t99 * t63); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t154; m(5) * t28 + m(6) * t21; m(6) * (t21 * t11 + t39 * t26 + t40 * t27) + m(5) * (t28 * t18 + (-t100 * t45 - t44 * t99) * t81) + t108; m(6) * (-t39 * t100 + t40 * t99); m(5) * (t154 * t81 ^ 2 + t28 ^ 2) + m(6) * (t21 ^ 2 + t39 ^ 2 + t40 ^ 2) + t108; m(6) * t22; m(6) * (t22 * t11 + t24 * t27 + t25 * t26) + t129; m(6) * (-t25 * t100 + t24 * t99); m(6) * (t22 * t21 + t24 * t40 + t25 * t39) + t129; m(6) * (t22 ^ 2 + t24 ^ 2 + t25 ^ 2) + t4 * t137 + t3 * t144 - t89 * (t89 ^ 2 * t48 + (t20 * t100 + t19 * t99 - (-t102 * t49 + t104 * t50) * t89) * t88);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
