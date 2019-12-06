% Calculate joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:45
% DurationCPUTime: 0.42s
% Computational Cost: add. (1570->111), mult. (862->147), div. (0->0), fcn. (686->10), ass. (0->64)
t58 = qJ(1) + qJ(2);
t57 = qJ(3) + t58;
t52 = pkin(9) + t57;
t49 = sin(t52);
t88 = t49 ^ 2;
t50 = cos(t52);
t87 = t50 ^ 2;
t59 = sin(qJ(5));
t91 = Icges(6,5) * t59;
t61 = cos(qJ(5));
t90 = Icges(6,6) * t61;
t39 = t90 + t91;
t89 = t49 * t50;
t42 = t59 * rSges(6,1) + t61 * rSges(6,2);
t84 = m(6) * t42;
t55 = sin(t58);
t83 = pkin(2) * t55;
t56 = cos(t58);
t82 = pkin(2) * t56;
t53 = sin(t57);
t81 = pkin(3) * t53;
t54 = cos(t57);
t80 = pkin(3) * t54;
t60 = sin(qJ(1));
t79 = t60 * pkin(1);
t62 = cos(qJ(1));
t78 = t62 * pkin(1);
t77 = rSges(6,1) * t61;
t76 = rSges(6,2) * t59;
t75 = t50 * rSges(6,3) + t49 * t76;
t72 = t39 * t87 + (t91 / 0.2e1 + t90 / 0.2e1 + t39 / 0.2e1) * t88;
t71 = -pkin(4) - t77;
t33 = -t56 * rSges(3,1) + t55 * rSges(3,2);
t31 = -t54 * rSges(4,1) + t53 * rSges(4,2);
t70 = Icges(6,2) * t61 ^ 2 + Icges(4,3) + Icges(5,3) + (Icges(6,1) * t59 + 0.2e1 * Icges(6,4) * t61) * t59;
t32 = -t55 * rSges(3,1) - t56 * rSges(3,2);
t30 = -t53 * rSges(4,1) - t54 * rSges(4,2);
t66 = Icges(3,3) + t70;
t63 = Icges(6,5) * t61 - Icges(6,6) * t59;
t27 = t31 - t82;
t23 = -t50 * rSges(5,1) + t49 * rSges(5,2) - t80;
t26 = t30 - t83;
t22 = -t49 * rSges(5,1) - t50 * rSges(5,2) - t81;
t21 = t23 - t82;
t20 = t22 - t83;
t10 = t50 * pkin(8) + t71 * t49 + t75 - t81;
t37 = t50 * t76;
t11 = -t80 + t37 + t71 * t50 + (-rSges(6,3) - pkin(8)) * t49;
t8 = t10 - t83;
t9 = t11 - t82;
t44 = -t62 * rSges(2,1) + t60 * rSges(2,2);
t43 = -t60 * rSges(2,1) - t62 * rSges(2,2);
t29 = t33 - t78;
t28 = t32 - t79;
t25 = t27 - t78;
t24 = t26 - t79;
t15 = Icges(6,3) * t49 + t63 * t50;
t14 = Icges(6,3) * t50 - t63 * t49;
t13 = t21 - t78;
t12 = t20 - t79;
t7 = t9 - t78;
t6 = t8 - t79;
t3 = t50 * (t49 * rSges(6,3) + t50 * t77 - t37) - t49 * (-t49 * t77 + t75);
t1 = [Icges(2,3) + m(2) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t66; m(3) * (t32 * t28 + t33 * t29) + m(4) * (t26 * t24 + t27 * t25) + m(5) * (t20 * t12 + t21 * t13) + m(6) * (t8 * t6 + t9 * t7) + t66; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + t66; m(4) * (t30 * t24 + t31 * t25) + m(5) * (t22 * t12 + t23 * t13) + m(6) * (t10 * t6 + t11 * t7) + t70; m(6) * (t10 * t8 + t11 * t9) + m(5) * (t22 * t20 + t23 * t21) + m(4) * (t30 * t26 + t31 * t27) + t70; m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t70; 0; 0; 0; m(5) + m(6); (t49 * t7 - t50 * t6) * t84 + t72; (t49 * t9 - t50 * t8) * t84 + t72; (-t10 * t50 + t11 * t49) * t84 + t72; m(6) * t3; m(6) * (t3 ^ 2 + (t87 + t88) * t42 ^ 2) + t50 * (t87 * t14 + t15 * t89) + t49 * (t14 * t89 + t88 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
