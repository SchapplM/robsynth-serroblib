% Calculate joint inertia matrix for
% S5RRPPR1
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:45
% DurationCPUTime: 0.56s
% Computational Cost: add. (1175->99), mult. (718->138), div. (0->0), fcn. (582->10), ass. (0->58)
t73 = qJ(1) + qJ(2);
t67 = pkin(8) + t73;
t60 = sin(t67);
t55 = t60 ^ 2;
t61 = cos(t67);
t56 = t61 ^ 2;
t72 = pkin(9) + qJ(5);
t65 = sin(t72);
t104 = Icges(6,5) * t65;
t66 = cos(t72);
t103 = Icges(6,6) * t66;
t29 = t103 + t104;
t102 = -rSges(6,1) * t66 + rSges(6,2) * t65;
t101 = rSges(5,3) + qJ(4);
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t100 = rSges(5,1) * t75 - rSges(5,2) * t74 + pkin(3);
t99 = t60 * t61;
t35 = t65 * rSges(6,1) + t66 * rSges(6,2);
t96 = m(6) * t35;
t91 = t56 + t55;
t68 = sin(t73);
t69 = cos(t73);
t36 = t68 * rSges(3,1) + t69 * rSges(3,2);
t63 = pkin(2) * t68;
t22 = t60 * rSges(4,1) + t61 * rSges(4,2) + t63;
t88 = t29 * t56 + (t104 / 0.2e1 + t103 / 0.2e1 + t29 / 0.2e1) * t55;
t37 = t69 * rSges(3,1) - t68 * rSges(3,2);
t64 = pkin(2) * t69;
t23 = t61 * rSges(4,1) - t60 * rSges(4,2) + t64;
t82 = Icges(6,5) * t66 - Icges(6,6) * t65;
t81 = -t60 * rSges(6,3) + t102 * t61;
t80 = Icges(5,2) * t75 ^ 2 + Icges(6,2) * t66 ^ 2 + Icges(3,3) + Icges(4,3) + (Icges(5,1) * t74 + 0.2e1 * Icges(5,4) * t75) * t74 + (Icges(6,1) * t65 + 0.2e1 * Icges(6,4) * t66) * t65;
t79 = -t61 * rSges(6,3) - t102 * t60;
t13 = t100 * t61 + t101 * t60 + t64;
t62 = t75 * pkin(4) + pkin(3);
t76 = -pkin(7) - qJ(4);
t8 = t60 * t62 + t61 * t76 + t63 + t79;
t9 = -t60 * t76 + t61 * t62 + t64 - t81;
t12 = t100 * t60 - t101 * t61 + t63;
t78 = cos(qJ(1));
t77 = sin(qJ(1));
t71 = t78 * pkin(1);
t70 = t77 * pkin(1);
t46 = t78 * rSges(2,1) - t77 * rSges(2,2);
t45 = t77 * rSges(2,1) + t78 * rSges(2,2);
t27 = t37 + t71;
t26 = t70 + t36;
t21 = t23 + t71;
t20 = t70 + t22;
t15 = -Icges(6,3) * t60 - t82 * t61;
t14 = -Icges(6,3) * t61 + t82 * t60;
t11 = t13 + t71;
t10 = t70 + t12;
t7 = t71 + t9;
t6 = t70 + t8;
t3 = t60 * t79 - t61 * t81;
t1 = [Icges(2,3) + m(2) * (t45 ^ 2 + t46 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t80; m(3) * (t36 * t26 + t37 * t27) + m(4) * (t22 * t20 + t23 * t21) + m(5) * (t12 * t10 + t13 * t11) + m(6) * (t8 * t6 + t9 * t7) + t80; m(6) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + t80; 0; 0; m(4) + m(5) + m(6); m(5) * (-t60 * t10 - t61 * t11) + m(6) * (-t60 * t6 - t61 * t7); m(6) * (-t60 * t8 - t61 * t9) + m(5) * (-t60 * t12 - t61 * t13); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t91; (t6 * t61 - t60 * t7) * t96 + t88; (-t60 * t9 + t61 * t8) * t96 + t88; m(6) * t3; 0; m(6) * (t91 * t35 ^ 2 + t3 ^ 2) - t61 * (t56 * t14 + t15 * t99) - t60 * (t14 * t99 + t55 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
