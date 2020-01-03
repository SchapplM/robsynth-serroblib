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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:23
% DurationCPUTime: 0.41s
% Computational Cost: add. (1570->109), mult. (862->145), div. (0->0), fcn. (686->10), ass. (0->62)
t72 = qJ(1) + qJ(2);
t69 = qJ(3) + t72;
t64 = pkin(9) + t69;
t55 = sin(t64);
t96 = t55 ^ 2;
t56 = cos(t64);
t95 = t56 ^ 2;
t73 = sin(qJ(5));
t100 = Icges(6,5) * t73;
t75 = cos(qJ(5));
t99 = Icges(6,6) * t75;
t39 = t99 + t100;
t98 = -rSges(6,1) * t75 + rSges(6,2) * t73;
t97 = t55 * t56;
t42 = t73 * rSges(6,1) + t75 * rSges(6,2);
t92 = m(6) * t42;
t65 = sin(t69);
t66 = cos(t69);
t30 = t65 * rSges(4,1) + t66 * rSges(4,2);
t67 = sin(t72);
t68 = cos(t72);
t32 = t67 * rSges(3,1) + t68 * rSges(3,2);
t57 = pkin(3) * t65;
t22 = t55 * rSges(5,1) + t56 * rSges(5,2) + t57;
t62 = pkin(2) * t67;
t26 = t62 + t30;
t87 = t39 * t95 + (t100 / 0.2e1 + t99 / 0.2e1 + t39 / 0.2e1) * t96;
t33 = t68 * rSges(3,1) - t67 * rSges(3,2);
t31 = t66 * rSges(4,1) - t65 * rSges(4,2);
t20 = t62 + t22;
t86 = Icges(6,2) * t75 ^ 2 + Icges(4,3) + Icges(5,3) + (Icges(6,1) * t73 + 0.2e1 * Icges(6,4) * t75) * t73;
t85 = t98 * t55;
t63 = pkin(2) * t68;
t27 = t31 + t63;
t58 = pkin(3) * t66;
t23 = t56 * rSges(5,1) - t55 * rSges(5,2) + t58;
t81 = Icges(3,3) + t86;
t78 = Icges(6,5) * t75 - Icges(6,6) * t73;
t77 = -t55 * rSges(6,3) + t98 * t56;
t21 = t23 + t63;
t11 = t56 * pkin(4) + t55 * pkin(8) + t58 - t77;
t9 = t11 + t63;
t10 = t55 * pkin(4) + t57 + (-rSges(6,3) - pkin(8)) * t56 - t85;
t8 = t62 + t10;
t76 = cos(qJ(1));
t74 = sin(qJ(1));
t71 = t76 * pkin(1);
t70 = t74 * pkin(1);
t44 = t76 * rSges(2,1) - t74 * rSges(2,2);
t43 = t74 * rSges(2,1) + t76 * rSges(2,2);
t29 = t33 + t71;
t28 = t70 + t32;
t25 = t27 + t71;
t24 = t70 + t26;
t15 = -Icges(6,3) * t55 - t78 * t56;
t14 = -Icges(6,3) * t56 + t78 * t55;
t13 = t21 + t71;
t12 = t70 + t20;
t7 = t71 + t9;
t6 = t70 + t8;
t3 = -t56 * t77 + t55 * (-t56 * rSges(6,3) - t85);
t1 = [Icges(2,3) + m(2) * (t43 ^ 2 + t44 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + t81; m(3) * (t32 * t28 + t33 * t29) + m(4) * (t26 * t24 + t27 * t25) + m(5) * (t20 * t12 + t21 * t13) + m(6) * (t8 * t6 + t9 * t7) + t81; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + t81; m(4) * (t30 * t24 + t31 * t25) + m(5) * (t22 * t12 + t23 * t13) + m(6) * (t10 * t6 + t11 * t7) + t86; m(6) * (t10 * t8 + t11 * t9) + m(5) * (t22 * t20 + t23 * t21) + m(4) * (t30 * t26 + t31 * t27) + t86; m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t86; 0; 0; 0; m(5) + m(6); (-t55 * t7 + t56 * t6) * t92 + t87; (-t55 * t9 + t56 * t8) * t92 + t87; (t10 * t56 - t11 * t55) * t92 + t87; m(6) * t3; m(6) * (t3 ^ 2 + (t95 + t96) * t42 ^ 2) - t56 * (t95 * t14 + t15 * t97) - t55 * (t14 * t97 + t96 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
