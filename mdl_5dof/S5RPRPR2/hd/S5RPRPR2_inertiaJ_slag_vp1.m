% Calculate joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:47
% EndTime: 2022-01-23 09:18:48
% DurationCPUTime: 0.38s
% Computational Cost: add. (1114->94), mult. (674->133), div. (0->0), fcn. (550->10), ass. (0->57)
t64 = qJ(1) + pkin(8);
t61 = qJ(3) + t64;
t54 = sin(t61);
t51 = t54 ^ 2;
t55 = cos(t61);
t52 = t55 ^ 2;
t63 = pkin(9) + qJ(5);
t57 = sin(t63);
t97 = Icges(6,5) * t57;
t59 = cos(t63);
t96 = Icges(6,6) * t59;
t29 = t96 + t97;
t95 = qJ(4) + rSges(5,3);
t65 = sin(pkin(9));
t66 = cos(pkin(9));
t94 = -rSges(5,1) * t66 + rSges(5,2) * t65 - pkin(3);
t93 = t54 * t55;
t35 = t57 * rSges(6,1) + t59 * rSges(6,2);
t90 = m(6) * t35;
t68 = sin(qJ(1));
t89 = t68 * pkin(1);
t87 = rSges(6,1) * t59;
t85 = rSges(6,2) * t57;
t84 = t55 * rSges(6,3) + t54 * t85;
t83 = t51 + t52;
t60 = cos(t64);
t69 = cos(qJ(1));
t62 = t69 * pkin(1);
t82 = pkin(2) * t60 + t62;
t79 = t29 * t52 + (t97 / 0.2e1 + t96 / 0.2e1 + t29 / 0.2e1) * t51;
t27 = t55 * rSges(4,1) - t54 * rSges(4,2);
t58 = sin(t64);
t78 = -pkin(2) * t58 - t89;
t77 = Icges(5,2) * t66 ^ 2 + Icges(6,2) * t59 ^ 2 + Icges(4,3) + (Icges(5,1) * t65 + 0.2e1 * Icges(5,4) * t66) * t65 + (Icges(6,1) * t57 + 0.2e1 * Icges(6,4) * t59) * t57;
t26 = -t54 * rSges(4,1) - t55 * rSges(4,2);
t71 = Icges(6,5) * t59 - Icges(6,6) * t57;
t70 = t54 * rSges(6,3) + (-t85 + t87) * t55;
t13 = t95 * t54 - t94 * t55;
t12 = t94 * t54 + t95 * t55;
t56 = t66 * pkin(4) + pkin(3);
t67 = -pkin(7) - qJ(4);
t11 = -t54 * t67 + t55 * t56 + t70;
t10 = -t55 * t67 + (-t56 - t87) * t54 + t84;
t42 = t69 * rSges(2,1) - t68 * rSges(2,2);
t41 = -t68 * rSges(2,1) - t69 * rSges(2,2);
t25 = t60 * rSges(3,1) - t58 * rSges(3,2) + t62;
t24 = -t58 * rSges(3,1) - t60 * rSges(3,2) - t89;
t21 = t27 + t82;
t20 = t26 + t78;
t15 = Icges(6,3) * t54 + t71 * t55;
t14 = -Icges(6,3) * t55 + t71 * t54;
t9 = t13 + t82;
t8 = t12 + t78;
t7 = t11 + t82;
t6 = t10 + t78;
t3 = t54 * (t54 * t87 - t84) + t55 * t70;
t1 = [Icges(2,3) + Icges(3,3) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t24 ^ 2 + t25 ^ 2) + m(2) * (t41 ^ 2 + t42 ^ 2) + t77; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t10 * t6 + t11 * t7) + m(5) * (t12 * t8 + t13 * t9) + m(4) * (t26 * t20 + t27 * t21) + t77; 0; m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + t77; m(6) * (t54 * t6 - t55 * t7) + m(5) * (t54 * t8 - t55 * t9); 0; m(6) * (t54 * t10 - t55 * t11) + m(5) * (t54 * t12 - t55 * t13); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t83; (-t54 * t7 - t55 * t6) * t90 + t79; m(6) * t3; (-t10 * t55 - t11 * t54) * t90 + t79; 0; m(6) * (t83 * t35 ^ 2 + t3 ^ 2) + t54 * (-t14 * t93 + t51 * t15) - t55 * (t52 * t14 - t15 * t93);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
