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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:05:26
% EndTime: 2022-01-20 10:05:26
% DurationCPUTime: 0.45s
% Computational Cost: add. (1758->147), mult. (1474->218), div. (0->0), fcn. (1472->10), ass. (0->73)
t72 = qJ(1) + qJ(2);
t69 = sin(t72);
t96 = pkin(2) * t69;
t76 = sin(qJ(1));
t95 = t76 * pkin(1);
t73 = sin(pkin(9));
t74 = cos(pkin(9));
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t40 = -Icges(6,3) * t74 + (Icges(6,5) * t77 - Icges(6,6) * t75) * t73;
t41 = -Icges(6,6) * t74 + (Icges(6,4) * t77 - Icges(6,2) * t75) * t73;
t42 = -Icges(6,5) * t74 + (Icges(6,1) * t77 - Icges(6,4) * t75) * t73;
t16 = -t74 * t40 + (-t41 * t75 + t42 * t77) * t73;
t94 = t16 * t74;
t68 = pkin(8) + t72;
t65 = sin(t68);
t93 = t65 * t73;
t66 = cos(t68);
t92 = t66 * t73;
t91 = t66 * t74;
t90 = t74 * t75;
t89 = t74 * t77;
t88 = Icges(6,5) * t73;
t87 = Icges(6,6) * t73;
t86 = Icges(6,3) * t73;
t38 = t65 * t77 - t66 * t90;
t39 = t65 * t75 + t66 * t89;
t24 = t39 * rSges(6,1) + t38 * rSges(6,2) + rSges(6,3) * t92;
t70 = cos(t72);
t67 = pkin(2) * t70;
t85 = t66 * pkin(3) + t65 * qJ(4) + t67;
t82 = t66 * qJ(4) - t96;
t50 = t70 * rSges(3,1) - t69 * rSges(3,2);
t35 = t66 * rSges(4,1) - t65 * rSges(4,2) + t67;
t49 = -t69 * rSges(3,1) - t70 * rSges(3,2);
t36 = -t65 * t90 - t66 * t77;
t37 = t65 * t89 - t66 * t75;
t81 = -t37 * rSges(6,1) - t36 * rSges(6,2);
t17 = Icges(6,5) * t37 + Icges(6,6) * t36 + t65 * t86;
t19 = Icges(6,4) * t37 + Icges(6,2) * t36 + t65 * t87;
t21 = Icges(6,1) * t37 + Icges(6,4) * t36 + t65 * t88;
t3 = -t74 * t17 + (-t19 * t75 + t21 * t77) * t73;
t18 = Icges(6,5) * t39 + Icges(6,6) * t38 + t66 * t86;
t20 = Icges(6,4) * t39 + Icges(6,2) * t38 + t66 * t87;
t22 = Icges(6,1) * t39 + Icges(6,4) * t38 + t66 * t88;
t4 = -t74 * t18 + (-t20 * t75 + t22 * t77) * t73;
t8 = t36 * t41 + t37 * t42 + t40 * t93;
t9 = t38 * t41 + t39 * t42 + t40 * t92;
t80 = -t94 + (t3 + t8) * t93 / 0.2e1 + (t4 + t9) * t92 / 0.2e1;
t34 = -t65 * rSges(4,1) - t66 * rSges(4,2) - t96;
t13 = pkin(4) * t91 + pkin(7) * t92 + t24 + t85;
t28 = rSges(5,1) * t91 - rSges(5,2) * t92 + t65 * rSges(5,3) + t85;
t27 = rSges(5,2) * t93 + t66 * rSges(5,3) + (-rSges(5,1) * t74 - pkin(3)) * t65 + t82;
t79 = Icges(5,2) * t74 ^ 2 + Icges(3,3) + Icges(4,3) + t16 + (Icges(5,1) * t73 + 0.2e1 * Icges(5,4) * t74) * t73;
t12 = (-pkin(4) * t74 - pkin(3) + (-rSges(6,3) - pkin(7)) * t73) * t65 + t81 + t82;
t78 = cos(qJ(1));
t71 = t78 * pkin(1);
t57 = t78 * rSges(2,1) - t76 * rSges(2,2);
t56 = -t76 * rSges(2,1) - t78 * rSges(2,2);
t45 = t50 + t71;
t44 = t49 - t95;
t43 = -t74 * rSges(6,3) + (rSges(6,1) * t77 - rSges(6,2) * t75) * t73;
t31 = t35 + t71;
t30 = t34 - t95;
t26 = t28 + t71;
t25 = t27 - t95;
t23 = rSges(6,3) * t93 - t81;
t15 = -t74 * t24 - t43 * t92;
t14 = t74 * t23 + t43 * t93;
t11 = t71 + t13;
t10 = t12 - t95;
t5 = (t23 * t66 - t24 * t65) * t73;
t1 = [Icges(2,3) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(2) * (t56 ^ 2 + t57 ^ 2) + t79; m(6) * (t12 * t10 + t13 * t11) + m(4) * (t34 * t30 + t35 * t31) + m(5) * (t27 * t25 + t28 * t26) + m(3) * (t49 * t44 + t50 * t45) + t79; m(6) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(3) * (t49 ^ 2 + t50 ^ 2) + t79; 0; 0; m(4) + m(5) + m(6); m(6) * (t65 * t10 - t66 * t11) + m(5) * (t65 * t25 - t66 * t26); m(6) * (t65 * t12 - t66 * t13) + m(5) * (t65 * t27 - t66 * t28); 0; 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t65 ^ 2 + t66 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t80; m(6) * (t14 * t12 + t15 * t13) + t80; m(6) * t5; m(6) * (t14 * t65 - t15 * t66); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) + ((t18 * t92 + t38 * t20 + t39 * t22) * t92 + (t17 * t92 + t38 * t19 + t39 * t21) * t93 - t9 * t74) * t92 + ((t18 * t93 + t36 * t20 + t37 * t22) * t92 + (t17 * t93 + t36 * t19 + t37 * t21) * t93 - t8 * t74) * t93 - t74 * (-t94 + (t3 * t65 + t4 * t66) * t73);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
