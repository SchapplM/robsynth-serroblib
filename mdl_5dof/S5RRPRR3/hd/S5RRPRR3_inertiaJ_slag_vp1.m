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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:33:57
% EndTime: 2022-01-20 10:33:58
% DurationCPUTime: 0.42s
% Computational Cost: add. (1460->102), mult. (804->137), div. (0->0), fcn. (646->10), ass. (0->61)
t62 = qJ(1) + qJ(2);
t58 = pkin(9) + t62;
t57 = qJ(4) + t58;
t51 = sin(t57);
t90 = t51 ^ 2;
t52 = cos(t57);
t89 = t52 ^ 2;
t63 = sin(qJ(5));
t93 = Icges(6,5) * t63;
t65 = cos(qJ(5));
t92 = Icges(6,6) * t65;
t37 = t92 + t93;
t91 = t51 * t52;
t40 = t63 * rSges(6,1) + t65 * rSges(6,2);
t86 = m(6) * t40;
t59 = sin(t62);
t85 = pkin(2) * t59;
t64 = sin(qJ(1));
t84 = t64 * pkin(1);
t83 = rSges(6,1) * t65;
t82 = rSges(6,2) * t63;
t81 = t52 * rSges(6,3) + t51 * t82;
t55 = cos(t58);
t60 = cos(t62);
t56 = pkin(2) * t60;
t80 = pkin(3) * t55 + t56;
t77 = Icges(6,2) * t65 ^ 2 + Icges(5,3) + (Icges(6,1) * t63 + 0.2e1 * Icges(6,4) * t65) * t63;
t76 = t37 * t89 + (t93 / 0.2e1 + t92 / 0.2e1 + t37 / 0.2e1) * t90;
t31 = t60 * rSges(3,1) - t59 * rSges(3,2);
t27 = t52 * rSges(5,1) - t51 * rSges(5,2);
t54 = sin(t58);
t25 = t55 * rSges(4,1) - t54 * rSges(4,2) + t56;
t75 = -pkin(3) * t54 - t85;
t30 = -t59 * rSges(3,1) - t60 * rSges(3,2);
t26 = -t51 * rSges(5,1) - t52 * rSges(5,2);
t71 = Icges(4,3) + Icges(3,3) + t77;
t68 = Icges(6,5) * t65 - Icges(6,6) * t63;
t67 = t51 * rSges(6,3) + (-t82 + t83) * t52;
t21 = t27 + t80;
t11 = t52 * pkin(4) + t51 * pkin(8) + t67;
t24 = -t54 * rSges(4,1) - t55 * rSges(4,2) - t85;
t10 = t52 * pkin(8) + (-pkin(4) - t83) * t51 + t81;
t9 = t11 + t80;
t20 = t26 + t75;
t8 = t10 + t75;
t66 = cos(qJ(1));
t61 = t66 * pkin(1);
t42 = t66 * rSges(2,1) - t64 * rSges(2,2);
t41 = -t64 * rSges(2,1) - t66 * rSges(2,2);
t29 = t31 + t61;
t28 = t30 - t84;
t23 = t25 + t61;
t22 = t24 - t84;
t15 = Icges(6,3) * t51 + t68 * t52;
t14 = -Icges(6,3) * t52 + t68 * t51;
t13 = t21 + t61;
t12 = t20 - t84;
t5 = t61 + t9;
t4 = t8 - t84;
t3 = t51 * (t51 * t83 - t81) + t52 * t67;
t1 = [Icges(2,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(2) * (t41 ^ 2 + t42 ^ 2) + t71; m(6) * (t8 * t4 + t9 * t5) + m(5) * (t20 * t12 + t21 * t13) + m(3) * (t30 * t28 + t31 * t29) + m(4) * (t24 * t22 + t25 * t23) + t71; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t30 ^ 2 + t31 ^ 2) + t71; 0; 0; m(4) + m(5) + m(6); m(6) * (t10 * t4 + t11 * t5) + m(5) * (t26 * t12 + t27 * t13) + t77; m(6) * (t10 * t8 + t11 * t9) + m(5) * (t26 * t20 + t27 * t21) + t77; 0; m(5) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t77; (-t4 * t52 - t5 * t51) * t86 + t76; (-t51 * t9 - t52 * t8) * t86 + t76; m(6) * t3; (-t10 * t52 - t11 * t51) * t86 + t76; m(6) * (t3 ^ 2 + (t89 + t90) * t40 ^ 2) + t51 * (-t14 * t91 + t90 * t15) - t52 * (t89 * t14 - t15 * t91);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
