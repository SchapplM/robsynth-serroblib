% Calculate joint inertia matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:21
% EndTime: 2020-01-03 11:57:22
% DurationCPUTime: 0.44s
% Computational Cost: add. (1758->146), mult. (1474->218), div. (0->0), fcn. (1472->10), ass. (0->73)
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t83 = sin(qJ(5));
t85 = cos(qJ(5));
t40 = -Icges(6,3) * t82 + (Icges(6,5) * t85 - Icges(6,6) * t83) * t81;
t41 = -Icges(6,6) * t82 + (Icges(6,4) * t85 - Icges(6,2) * t83) * t81;
t42 = -Icges(6,5) * t82 + (Icges(6,1) * t85 - Icges(6,4) * t83) * t81;
t16 = -t82 * t40 + (-t41 * t83 + t42 * t85) * t81;
t102 = t16 * t82;
t80 = qJ(1) + qJ(2);
t75 = pkin(8) + t80;
t71 = sin(t75);
t101 = t71 * t81;
t100 = t71 * t82;
t72 = cos(t75);
t99 = t72 * t81;
t98 = t72 * t82;
t97 = t82 * t83;
t96 = t82 * t85;
t76 = sin(t80);
t73 = pkin(2) * t76;
t95 = t71 * pkin(3) + t73;
t77 = cos(t80);
t49 = t76 * rSges(3,1) + t77 * rSges(3,2);
t94 = Icges(6,5) * t81;
t93 = Icges(6,6) * t81;
t92 = Icges(6,3) * t81;
t36 = -t71 * t97 - t72 * t85;
t37 = t71 * t96 - t72 * t83;
t23 = t37 * rSges(6,1) + t36 * rSges(6,2) + rSges(6,3) * t101;
t74 = pkin(2) * t77;
t91 = t72 * pkin(3) + t71 * qJ(4) + t74;
t34 = t71 * rSges(4,1) + t72 * rSges(4,2) + t73;
t50 = t77 * rSges(3,1) - t76 * rSges(3,2);
t35 = t72 * rSges(4,1) - t71 * rSges(4,2) + t74;
t17 = Icges(6,5) * t37 + Icges(6,6) * t36 + t71 * t92;
t19 = Icges(6,4) * t37 + Icges(6,2) * t36 + t71 * t93;
t21 = Icges(6,1) * t37 + Icges(6,4) * t36 + t71 * t94;
t3 = -t82 * t17 + (-t19 * t83 + t21 * t85) * t81;
t38 = -t71 * t85 + t72 * t97;
t39 = -t71 * t83 - t72 * t96;
t18 = Icges(6,5) * t39 + Icges(6,6) * t38 - t72 * t92;
t20 = Icges(6,4) * t39 + Icges(6,2) * t38 - t72 * t93;
t22 = Icges(6,1) * t39 + Icges(6,4) * t38 - t72 * t94;
t4 = -t82 * t18 + (-t20 * t83 + t22 * t85) * t81;
t8 = t40 * t101 + t36 * t41 + t37 * t42;
t9 = t38 * t41 + t39 * t42 - t40 * t99;
t88 = -t102 + (t3 + t8) * t101 / 0.2e1 - (t4 + t9) * t99 / 0.2e1;
t24 = t39 * rSges(6,1) + t38 * rSges(6,2) - rSges(6,3) * t99;
t28 = rSges(5,1) * t98 - rSges(5,2) * t99 + t71 * rSges(5,3) + t91;
t12 = pkin(4) * t100 + pkin(7) * t101 - t72 * qJ(4) + t23 + t95;
t27 = -rSges(5,2) * t101 + rSges(5,1) * t100 + (-rSges(5,3) - qJ(4)) * t72 + t95;
t13 = pkin(4) * t98 + pkin(7) * t99 - t24 + t91;
t87 = Icges(5,2) * t82 ^ 2 + Icges(3,3) + Icges(4,3) + t16 + (Icges(5,1) * t81 + 0.2e1 * Icges(5,4) * t82) * t81;
t86 = cos(qJ(1));
t84 = sin(qJ(1));
t79 = t86 * pkin(1);
t78 = t84 * pkin(1);
t60 = t86 * rSges(2,1) - t84 * rSges(2,2);
t59 = t84 * rSges(2,1) + t86 * rSges(2,2);
t45 = t50 + t79;
t44 = t78 + t49;
t43 = -t82 * rSges(6,3) + (rSges(6,1) * t85 - rSges(6,2) * t83) * t81;
t31 = t35 + t79;
t30 = t78 + t34;
t26 = t28 + t79;
t25 = t78 + t27;
t15 = t82 * t24 - t43 * t99;
t14 = -t43 * t101 - t82 * t23;
t11 = t13 + t79;
t10 = t12 + t78;
t5 = (t23 * t72 + t24 * t71) * t81;
t1 = [Icges(2,3) + m(2) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t87; m(3) * (t49 * t44 + t50 * t45) + m(4) * (t34 * t30 + t35 * t31) + m(5) * (t27 * t25 + t28 * t26) + m(6) * (t12 * t10 + t13 * t11) + t87; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(3) * (t49 ^ 2 + t50 ^ 2) + t87; 0; 0; m(4) + m(5) + m(6); m(5) * (-t71 * t25 - t72 * t26) + m(6) * (-t71 * t10 - t72 * t11); m(6) * (-t71 * t12 - t72 * t13) + m(5) * (-t71 * t27 - t72 * t28); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t71 ^ 2 + t72 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t88; m(6) * (t14 * t12 + t15 * t13) + t88; m(6) * t5; m(6) * (-t14 * t71 - t15 * t72); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) - t82 * (-t102 + (t3 * t71 - t4 * t72) * t81) + (-t8 * t82 + (t17 * t101 + t36 * t19 + t37 * t21) * t101 - (t18 * t101 + t36 * t20 + t37 * t22) * t99) * t101 - (-t9 * t82 + (-t17 * t99 + t38 * t19 + t39 * t21) * t101 - (-t18 * t99 + t38 * t20 + t39 * t22) * t99) * t99;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
