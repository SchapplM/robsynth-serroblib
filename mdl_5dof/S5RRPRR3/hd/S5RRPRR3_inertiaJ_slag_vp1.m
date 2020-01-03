% Calculate joint inertia matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:22
% DurationCPUTime: 0.40s
% Computational Cost: add. (1460->101), mult. (804->136), div. (0->0), fcn. (646->10), ass. (0->60)
t70 = qJ(1) + qJ(2);
t65 = pkin(9) + t70;
t64 = qJ(4) + t65;
t55 = sin(t64);
t96 = t55 ^ 2;
t56 = cos(t64);
t95 = t56 ^ 2;
t71 = sin(qJ(5));
t100 = Icges(6,5) * t71;
t73 = cos(qJ(5));
t99 = Icges(6,6) * t73;
t37 = t99 + t100;
t98 = -rSges(6,1) * t73 + rSges(6,2) * t71;
t97 = t55 * t56;
t40 = t71 * rSges(6,1) + t73 * rSges(6,2);
t92 = m(6) * t40;
t26 = t55 * rSges(5,1) + t56 * rSges(5,2);
t60 = sin(t65);
t66 = sin(t70);
t62 = pkin(2) * t66;
t89 = pkin(3) * t60 + t62;
t61 = cos(t65);
t67 = cos(t70);
t63 = pkin(2) * t67;
t88 = pkin(3) * t61 + t63;
t30 = t66 * rSges(3,1) + t67 * rSges(3,2);
t24 = t60 * rSges(4,1) + t61 * rSges(4,2) + t62;
t85 = Icges(6,2) * t73 ^ 2 + Icges(5,3) + (Icges(6,1) * t71 + 0.2e1 * Icges(6,4) * t73) * t71;
t84 = t37 * t95 + (t100 / 0.2e1 + t99 / 0.2e1 + t37 / 0.2e1) * t96;
t31 = t67 * rSges(3,1) - t66 * rSges(3,2);
t27 = t56 * rSges(5,1) - t55 * rSges(5,2);
t20 = t89 + t26;
t83 = t98 * t55;
t25 = t61 * rSges(4,1) - t60 * rSges(4,2) + t63;
t79 = Icges(3,3) + Icges(4,3) + t85;
t76 = Icges(6,5) * t73 - Icges(6,6) * t71;
t75 = -t55 * rSges(6,3) + t98 * t56;
t21 = t27 + t88;
t11 = t56 * pkin(4) + t55 * pkin(8) - t75;
t9 = t11 + t88;
t10 = t55 * pkin(4) + (-rSges(6,3) - pkin(8)) * t56 - t83;
t8 = t10 + t89;
t74 = cos(qJ(1));
t72 = sin(qJ(1));
t69 = t74 * pkin(1);
t68 = t72 * pkin(1);
t42 = t74 * rSges(2,1) - t72 * rSges(2,2);
t41 = t72 * rSges(2,1) + t74 * rSges(2,2);
t29 = t31 + t69;
t28 = t68 + t30;
t23 = t25 + t69;
t22 = t68 + t24;
t15 = -Icges(6,3) * t55 - t76 * t56;
t14 = -Icges(6,3) * t56 + t76 * t55;
t13 = t21 + t69;
t12 = t68 + t20;
t5 = t69 + t9;
t4 = t68 + t8;
t3 = -t56 * t75 + t55 * (-t56 * rSges(6,3) - t83);
t1 = [Icges(2,3) + m(2) * (t41 ^ 2 + t42 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2) + t79; m(3) * (t30 * t28 + t31 * t29) + m(4) * (t24 * t22 + t25 * t23) + m(5) * (t20 * t12 + t21 * t13) + m(6) * (t8 * t4 + t9 * t5) + t79; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t30 ^ 2 + t31 ^ 2) + t79; 0; 0; m(4) + m(5) + m(6); m(5) * (t26 * t12 + t27 * t13) + m(6) * (t10 * t4 + t11 * t5) + t85; m(6) * (t10 * t8 + t11 * t9) + m(5) * (t26 * t20 + t27 * t21) + t85; 0; m(5) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t85; (t4 * t56 - t5 * t55) * t92 + t84; (-t55 * t9 + t56 * t8) * t92 + t84; m(6) * t3; (t10 * t56 - t11 * t55) * t92 + t84; m(6) * (t3 ^ 2 + (t95 + t96) * t40 ^ 2) - t56 * (t95 * t14 + t15 * t97) - t55 * (t14 * t97 + t96 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
