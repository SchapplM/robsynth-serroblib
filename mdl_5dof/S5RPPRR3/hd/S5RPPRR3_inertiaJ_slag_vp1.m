% Calculate joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:07
% EndTime: 2022-01-23 09:14:08
% DurationCPUTime: 0.52s
% Computational Cost: add. (1611->143), mult. (1054->210), div. (0->0), fcn. (916->10), ass. (0->73)
t68 = qJ(1) + pkin(8);
t60 = sin(t68);
t62 = cos(t68);
t107 = t60 * t62;
t67 = pkin(9) + qJ(4);
t63 = qJ(5) + t67;
t54 = sin(t63);
t55 = cos(t63);
t87 = rSges(6,1) * t55 - rSges(6,2) * t54;
t57 = t60 ^ 2;
t58 = t62 ^ 2;
t106 = t60 / 0.2e1;
t105 = -t62 / 0.2e1;
t72 = sin(qJ(1));
t104 = t72 * pkin(1);
t70 = cos(pkin(9));
t56 = t70 * pkin(3) + pkin(2);
t61 = cos(t67);
t103 = rSges(5,1) * t61;
t59 = sin(t67);
t101 = rSges(5,2) * t59;
t71 = -pkin(6) - qJ(3);
t75 = t60 * rSges(6,3) + t87 * t62;
t8 = t60 * (-t62 * rSges(6,3) + t87 * t60) + t62 * t75;
t46 = pkin(4) * t61 + t56;
t66 = pkin(7) - t71;
t99 = t62 * t46 + t66 * t60;
t98 = t60 * rSges(5,3) + t62 * t103;
t97 = t57 + t58;
t96 = Icges(5,4) * t59;
t95 = Icges(5,4) * t61;
t94 = Icges(6,4) * t54;
t93 = Icges(6,4) * t55;
t92 = rSges(4,3) + qJ(3);
t36 = Icges(6,5) * t54 + Icges(6,6) * t55;
t78 = -Icges(6,2) * t54 + t93;
t80 = Icges(6,1) * t55 - t94;
t37 = Icges(6,2) * t55 + t94;
t38 = Icges(6,1) * t54 + t93;
t82 = -t37 * t54 + t38 * t55;
t91 = (t55 * (Icges(6,6) * t60 + t78 * t62) + t54 * (Icges(6,5) * t60 + t80 * t62) + t60 * t36 + t82 * t62) * t106 + (t55 * (-Icges(6,6) * t62 + t78 * t60) + t54 * (-Icges(6,5) * t62 + t80 * t60) - t62 * t36 + t82 * t60) * t105;
t76 = Icges(6,5) * t55 - Icges(6,6) * t54;
t20 = -Icges(6,3) * t62 + t76 * t60;
t21 = Icges(6,3) * t60 + t76 * t62;
t90 = -t62 * (-t21 * t107 + t58 * t20) + t60 * (-t20 * t107 + t57 * t21);
t39 = t54 * rSges(6,1) + t55 * rSges(6,2);
t89 = -pkin(4) * t59 - t39;
t88 = -t101 + t103;
t81 = Icges(5,1) * t61 - t96;
t79 = -Icges(5,2) * t59 + t95;
t77 = Icges(5,5) * t61 - Icges(5,6) * t59;
t69 = sin(pkin(9));
t74 = rSges(4,1) * t70 - rSges(4,2) * t69 + pkin(2);
t73 = cos(qJ(1));
t65 = t73 * pkin(1);
t50 = t73 * rSges(2,1) - t72 * rSges(2,2);
t49 = -t72 * rSges(2,1) - t73 * rSges(2,2);
t44 = t59 * rSges(5,1) + t61 * rSges(5,2);
t33 = t62 * rSges(3,1) - t60 * rSges(3,2) + t65;
t32 = -t60 * rSges(3,1) - t62 * rSges(3,2) - t104;
t27 = Icges(5,3) * t60 + t77 * t62;
t26 = -Icges(5,3) * t62 + t77 * t60;
t19 = t89 * t62;
t18 = t89 * t60;
t15 = t92 * t60 + t74 * t62 + t65;
t14 = -t74 * t60 + t92 * t62 - t104;
t13 = -t60 * t71 + t65 + (t56 - t101) * t62 + t98;
t12 = -t104 + (rSges(5,3) - t71) * t62 + (-t56 - t88) * t60;
t11 = t65 + t75 + t99;
t10 = -t104 + (rSges(6,3) + t66) * t62 + (-t46 - t87) * t60;
t9 = t62 * (-t62 * t101 + t98) + (-t62 * rSges(5,3) + t88 * t60) * t60;
t3 = t62 * (-t62 * t56 + t99) + (-t62 * t66 + (t46 - t56) * t60) * t60 + t8;
t1 = [Icges(4,2) * t70 ^ 2 + t55 * t37 + t54 * t38 + t61 * (Icges(5,2) * t61 + t96) + t59 * (Icges(5,1) * t59 + t95) + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t69 + 0.2e1 * Icges(4,4) * t70) * t69 + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(2) * (t49 ^ 2 + t50 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(6) * (t60 * t10 - t62 * t11) + m(5) * (t60 * t12 - t62 * t13) + m(4) * (t60 * t14 - t62 * t15); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t97; (t61 * (Icges(5,6) * t60 + t79 * t62) + t59 * (Icges(5,5) * t60 + t81 * t62)) * t106 + (t61 * (-Icges(5,6) * t62 + t79 * t60) + t59 * (-Icges(5,5) * t62 + t81 * t60)) * t105 + m(6) * (t19 * t10 + t18 * t11) + m(5) * (-t12 * t62 - t13 * t60) * t44 + (t57 / 0.2e1 + t58 / 0.2e1) * (Icges(5,5) * t59 + Icges(5,6) * t61) + t91; m(5) * t9 + m(6) * t3; m(6) * (-t18 * t62 + t19 * t60); m(5) * (t97 * t44 ^ 2 + t9 ^ 2) + t60 * (-t26 * t107 + t57 * t27) - t62 * (-t27 * t107 + t58 * t26) + m(6) * (t18 ^ 2 + t19 ^ 2 + t3 ^ 2) + t90; m(6) * (-t10 * t62 - t11 * t60) * t39 + t91; m(6) * t8; 0; m(6) * (t8 * t3 + (-t18 * t60 - t19 * t62) * t39) + t90; m(6) * (t97 * t39 ^ 2 + t8 ^ 2) + t90;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
