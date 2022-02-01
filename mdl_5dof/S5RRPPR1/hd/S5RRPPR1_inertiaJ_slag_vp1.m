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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:18
% EndTime: 2022-01-20 09:51:18
% DurationCPUTime: 0.43s
% Computational Cost: add. (1175->101), mult. (718->139), div. (0->0), fcn. (582->10), ass. (0->59)
t67 = qJ(1) + qJ(2);
t62 = pkin(8) + t67;
t56 = sin(t62);
t53 = t56 ^ 2;
t57 = cos(t62);
t54 = t57 ^ 2;
t66 = pkin(9) + qJ(5);
t60 = sin(t66);
t99 = Icges(6,5) * t60;
t61 = cos(t66);
t98 = Icges(6,6) * t61;
t29 = t98 + t99;
t97 = qJ(4) + rSges(5,3);
t68 = sin(pkin(9));
t69 = cos(pkin(9));
t96 = -rSges(5,1) * t69 + rSges(5,2) * t68 - pkin(3);
t95 = t56 * t57;
t35 = t60 * rSges(6,1) + t61 * rSges(6,2);
t92 = m(6) * t35;
t63 = sin(t67);
t91 = pkin(2) * t63;
t71 = sin(qJ(1));
t90 = t71 * pkin(1);
t88 = rSges(6,1) * t61;
t86 = rSges(6,2) * t60;
t85 = t57 * rSges(6,3) + t56 * t86;
t84 = t53 + t54;
t81 = t29 * t54 + (t99 / 0.2e1 + t98 / 0.2e1 + t29 / 0.2e1) * t53;
t64 = cos(t67);
t37 = t64 * rSges(3,1) - t63 * rSges(3,2);
t59 = pkin(2) * t64;
t23 = t57 * rSges(4,1) - t56 * rSges(4,2) + t59;
t36 = -t63 * rSges(3,1) - t64 * rSges(3,2);
t75 = Icges(6,5) * t61 - Icges(6,6) * t60;
t74 = t56 * rSges(6,3) + (-t86 + t88) * t57;
t73 = Icges(5,2) * t69 ^ 2 + Icges(6,2) * t61 ^ 2 + Icges(3,3) + Icges(4,3) + (Icges(5,1) * t68 + 0.2e1 * Icges(5,4) * t69) * t68 + (Icges(6,1) * t60 + 0.2e1 * Icges(6,4) * t61) * t60;
t22 = -t56 * rSges(4,1) - t57 * rSges(4,2) - t91;
t13 = t97 * t56 - t96 * t57 + t59;
t58 = t69 * pkin(4) + pkin(3);
t70 = -pkin(7) - qJ(4);
t9 = -t56 * t70 + t57 * t58 + t59 + t74;
t12 = t96 * t56 + t97 * t57 - t91;
t8 = -t91 - t57 * t70 + (-t58 - t88) * t56 + t85;
t72 = cos(qJ(1));
t65 = t72 * pkin(1);
t44 = t72 * rSges(2,1) - t71 * rSges(2,2);
t43 = -t71 * rSges(2,1) - t72 * rSges(2,2);
t27 = t37 + t65;
t26 = t36 - t90;
t21 = t23 + t65;
t20 = t22 - t90;
t15 = Icges(6,3) * t56 + t75 * t57;
t14 = -Icges(6,3) * t57 + t75 * t56;
t11 = t13 + t65;
t10 = t12 - t90;
t7 = t65 + t9;
t6 = t8 - t90;
t3 = t56 * (t56 * t88 - t85) + t57 * t74;
t1 = [Icges(2,3) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(2) * (t43 ^ 2 + t44 ^ 2) + t73; m(6) * (t8 * t6 + t9 * t7) + m(4) * (t22 * t20 + t23 * t21) + m(5) * (t12 * t10 + t13 * t11) + m(3) * (t36 * t26 + t37 * t27) + t73; m(6) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + t73; 0; 0; m(4) + m(5) + m(6); m(6) * (t56 * t6 - t57 * t7) + m(5) * (t56 * t10 - t57 * t11); m(6) * (t56 * t8 - t57 * t9) + m(5) * (t56 * t12 - t57 * t13); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t84; (-t56 * t7 - t57 * t6) * t92 + t81; (-t56 * t9 - t57 * t8) * t92 + t81; m(6) * t3; 0; m(6) * (t84 * t35 ^ 2 + t3 ^ 2) + t56 * (-t14 * t95 + t53 * t15) - t57 * (t54 * t14 - t15 * t95);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
