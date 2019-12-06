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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:11
% DurationCPUTime: 0.42s
% Computational Cost: add. (1175->101), mult. (718->140), div. (0->0), fcn. (582->10), ass. (0->59)
t61 = qJ(1) + qJ(2);
t57 = pkin(8) + t61;
t52 = sin(t57);
t49 = t52 ^ 2;
t53 = cos(t57);
t50 = t53 ^ 2;
t60 = pkin(9) + qJ(5);
t55 = sin(t60);
t96 = Icges(6,5) * t55;
t56 = cos(t60);
t95 = Icges(6,6) * t56;
t29 = t95 + t96;
t62 = sin(pkin(9));
t63 = cos(pkin(9));
t94 = -rSges(5,1) * t63 + rSges(5,2) * t62 - pkin(3);
t93 = rSges(5,3) + qJ(4);
t92 = t52 * t53;
t35 = t55 * rSges(6,1) + t56 * rSges(6,2);
t89 = m(6) * t35;
t58 = sin(t61);
t88 = pkin(2) * t58;
t59 = cos(t61);
t87 = pkin(2) * t59;
t65 = sin(qJ(1));
t86 = t65 * pkin(1);
t66 = cos(qJ(1));
t85 = t66 * pkin(1);
t84 = rSges(6,1) * t56;
t82 = rSges(6,2) * t55;
t81 = t53 * rSges(6,3) + t52 * t82;
t80 = t50 + t49;
t77 = t29 * t50 + (t96 / 0.2e1 + t95 / 0.2e1 + t29 / 0.2e1) * t49;
t37 = -t59 * rSges(3,1) + t58 * rSges(3,2);
t75 = -t63 * pkin(4) - pkin(3) - t84;
t74 = -t52 * rSges(6,3) + t53 * t82;
t36 = -t58 * rSges(3,1) - t59 * rSges(3,2);
t68 = Icges(6,5) * t56 - Icges(6,6) * t55;
t23 = -t53 * rSges(4,1) + t52 * rSges(4,2) - t87;
t67 = Icges(5,2) * t63 ^ 2 + Icges(6,2) * t56 ^ 2 + Icges(3,3) + Icges(4,3) + (Icges(5,1) * t62 + 0.2e1 * Icges(5,4) * t63) * t62 + (Icges(6,1) * t55 + 0.2e1 * Icges(6,4) * t56) * t55;
t22 = -t52 * rSges(4,1) - t53 * rSges(4,2) - t88;
t12 = t94 * t52 + t93 * t53 - t88;
t64 = -pkin(7) - qJ(4);
t8 = t75 * t52 - t53 * t64 + t81 - t88;
t9 = t52 * t64 + t75 * t53 + t74 - t87;
t13 = -t93 * t52 + t94 * t53 - t87;
t44 = -t66 * rSges(2,1) + t65 * rSges(2,2);
t43 = -t65 * rSges(2,1) - t66 * rSges(2,2);
t27 = t37 - t85;
t26 = t36 - t86;
t21 = t23 - t85;
t20 = t22 - t86;
t15 = Icges(6,3) * t52 + t68 * t53;
t14 = Icges(6,3) * t53 - t68 * t52;
t11 = t13 - t85;
t10 = t12 - t86;
t7 = t9 - t85;
t6 = t8 - t86;
t3 = t53 * (t53 * t84 - t74) - t52 * (-t52 * t84 + t81);
t1 = [Icges(2,3) + m(2) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t67; m(3) * (t36 * t26 + t37 * t27) + m(4) * (t22 * t20 + t23 * t21) + m(5) * (t12 * t10 + t13 * t11) + m(6) * (t8 * t6 + t9 * t7) + t67; m(6) * (t8 ^ 2 + t9 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + t67; 0; 0; m(4) + m(5) + m(6); m(5) * (t52 * t10 + t53 * t11) + m(6) * (t52 * t6 + t53 * t7); m(6) * (t52 * t8 + t53 * t9) + m(5) * (t52 * t12 + t53 * t13); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t80; (t52 * t7 - t53 * t6) * t89 + t77; (t52 * t9 - t53 * t8) * t89 + t77; m(6) * t3; 0; m(6) * (t80 * t35 ^ 2 + t3 ^ 2) + t53 * (t50 * t14 + t15 * t92) + t52 * (t14 * t92 + t49 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
