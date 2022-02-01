% Calculate joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:03
% EndTime: 2022-01-23 08:59:04
% DurationCPUTime: 0.60s
% Computational Cost: add. (1496->182), mult. (3579->299), div. (0->0), fcn. (4343->10), ass. (0->92)
t76 = cos(pkin(7));
t107 = 0.2e1 * t76;
t101 = m(6) / 0.2e1;
t102 = m(5) / 0.2e1;
t87 = t102 + t101;
t106 = 0.2e1 * t87;
t105 = t76 ^ 2;
t73 = sin(pkin(7));
t104 = 0.2e1 * t73;
t103 = m(4) / 0.2e1;
t100 = t73 * qJ(3) + pkin(1);
t72 = sin(pkin(8));
t99 = t72 * qJ(4) + pkin(2);
t77 = sin(qJ(5));
t98 = t72 * t77;
t79 = cos(qJ(5));
t97 = t72 * t79;
t75 = cos(pkin(8));
t96 = t73 * t75;
t78 = sin(qJ(1));
t95 = t73 * t78;
t80 = cos(qJ(1));
t94 = t73 * t80;
t74 = cos(pkin(9));
t93 = t76 * t74;
t92 = t78 * t72;
t91 = t78 * t75;
t90 = t80 * t72;
t89 = t80 * t75;
t88 = t78 ^ 2 + t80 ^ 2;
t71 = sin(pkin(9));
t51 = -t76 * t71 + t74 * t96;
t41 = t51 * t79 + t73 * t98;
t42 = -t51 * t77 + t73 * t97;
t50 = t71 * t96 + t93;
t24 = Icges(6,5) * t41 + Icges(6,6) * t42 + Icges(6,3) * t50;
t25 = Icges(6,4) * t41 + Icges(6,2) * t42 + Icges(6,6) * t50;
t26 = Icges(6,1) * t41 + Icges(6,4) * t42 + Icges(6,5) * t50;
t86 = t50 * t24 + t42 * t25 + t41 * t26;
t52 = t73 * t71 + t75 * t93;
t34 = (t52 * t79 + t76 * t98) * t80 + t78 * (t74 * t97 - t77 * t75);
t43 = -t52 * t77 + t76 * t97;
t53 = t74 * t98 + t79 * t75;
t36 = t43 * t80 - t78 * t53;
t57 = t76 * t89 + t92;
t45 = t57 * t71 - t74 * t94;
t18 = t34 * rSges(6,1) + t36 * rSges(6,2) + t45 * rSges(6,3);
t85 = rSges(4,3) * t73 + pkin(2) * t76 + t100;
t84 = qJ(4) * t75 - qJ(2);
t83 = t103 + t87;
t82 = rSges(3,1) * t76 - rSges(3,2) * t73 + pkin(1);
t60 = t74 * pkin(4) + t71 * pkin(6) + pkin(3);
t81 = -t60 * t72 + t84;
t54 = t76 * t92 + t89;
t33 = (t52 * t78 - t74 * t90) * t79 + t54 * t77;
t35 = t43 * t78 + t80 * t53;
t55 = t76 * t91 - t90;
t44 = t55 * t71 - t74 * t95;
t17 = t33 * rSges(6,1) + t35 * rSges(6,2) + t44 * rSges(6,3);
t68 = t73 ^ 2;
t67 = t80 * qJ(2);
t66 = t78 * qJ(2);
t62 = t80 * rSges(2,1) - t78 * rSges(2,2);
t61 = -t78 * rSges(2,1) - t80 * rSges(2,2);
t58 = -t72 * pkin(3) + t84;
t56 = t76 * t90 - t91;
t48 = (t75 * pkin(3) + t99) * t76 + t100;
t47 = t78 * rSges(3,3) + t82 * t80 + t66;
t46 = t80 * rSges(3,3) - t82 * t78 + t67;
t37 = (t60 * t75 + t99) * t76 + pkin(1) + (t71 * pkin(4) - t74 * pkin(6) + qJ(3)) * t73;
t29 = -t55 * rSges(4,1) + t54 * rSges(4,2) - t85 * t78 + t67;
t28 = t57 * rSges(4,1) - t56 * rSges(4,2) + t85 * t80 + t66;
t27 = t41 * rSges(6,1) + t42 * rSges(6,2) + t50 * rSges(6,3);
t20 = t48 * t80 - t58 * t78 + (t57 * t74 + t71 * t94) * rSges(5,1) - t45 * rSges(5,2) + t56 * rSges(5,3);
t19 = -t48 * t78 - t58 * t80 - (t55 * t74 + t71 * t95) * rSges(5,1) + t44 * rSges(5,2) - t54 * rSges(5,3);
t16 = Icges(6,1) * t34 + Icges(6,4) * t36 + Icges(6,5) * t45;
t15 = Icges(6,1) * t33 + Icges(6,4) * t35 + Icges(6,5) * t44;
t14 = Icges(6,4) * t34 + Icges(6,2) * t36 + Icges(6,6) * t45;
t13 = Icges(6,4) * t33 + Icges(6,2) * t35 + Icges(6,6) * t44;
t12 = Icges(6,5) * t34 + Icges(6,6) * t36 + Icges(6,3) * t45;
t11 = Icges(6,5) * t33 + Icges(6,6) * t35 + Icges(6,3) * t44;
t10 = t37 * t80 - t81 * t78 + t18;
t9 = -t37 * t78 - t81 * t80 - t17;
t8 = t50 * t18 - t45 * t27;
t7 = -t50 * t17 + t44 * t27;
t6 = t45 * t17 - t44 * t18;
t5 = t86 * t50;
t4 = t45 * t24 + t36 * t25 + t34 * t26;
t3 = t44 * t24 + t35 * t25 + t33 * t26;
t2 = t50 * t12 + t42 * t14 + t41 * t16;
t1 = t50 * t11 + t42 * t13 + t41 * t15;
t21 = [Icges(5,1) * t51 ^ 2 + Icges(2,3) + (-0.2e1 * Icges(5,4) * t51 + Icges(5,2) * t50) * t50 + (Icges(3,2) + Icges(4,3)) * t105 + m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t46 ^ 2 + t47 ^ 2) + m(2) * (t61 ^ 2 + t62 ^ 2) + ((Icges(4,1) * t75 ^ 2 + Icges(3,1)) * t73 + (-Icges(4,5) * t75 + Icges(3,4)) * t107 + (0.2e1 * Icges(5,5) * t51 + Icges(4,6) * t107 - 0.2e1 * Icges(5,6) * t50 + (-0.2e1 * Icges(4,4) * t75 + (Icges(4,2) + Icges(5,3)) * t72) * t73) * t72) * t73 + t86; m(6) * (-t80 * t10 + t78 * t9) + m(5) * (t78 * t19 - t80 * t20) + m(4) * (-t80 * t28 + t78 * t29) + m(3) * (t78 * t46 - t80 * t47); 0.2e1 * (m(3) / 0.2e1 + t83) * t88; ((t10 * t78 + t80 * t9) * t101 + (t19 * t80 + t20 * t78) * t102 + (t28 * t78 + t29 * t80) * t103) * t104; 0; 0.2e1 * t83 * (t88 * t68 + t105); m(6) * (t54 * t10 + t56 * t9) + m(5) * (t56 * t19 + t54 * t20); (-t54 * t80 + t56 * t78) * t106; t87 * (t54 * t78 + t56 * t80 - t72 * t76) * t104; (t68 * t72 ^ 2 + t54 ^ 2 + t56 ^ 2) * t106; t5 + m(6) * (t8 * t10 + t7 * t9) + (t2 / 0.2e1 + t4 / 0.2e1) * t45 + (t3 / 0.2e1 + t1 / 0.2e1) * t44; m(6) * (t7 * t78 - t8 * t80); m(6) * (-t6 * t76 + (t7 * t80 + t78 * t8) * t73); m(6) * (t6 * t73 * t72 + t8 * t54 + t7 * t56); m(6) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) + t45 * ((t45 * t12 + t36 * t14 + t34 * t16) * t45 + (t45 * t11 + t36 * t13 + t34 * t15) * t44 + t4 * t50) + t44 * ((t44 * t12 + t35 * t14 + t33 * t16) * t45 + (t44 * t11 + t35 * t13 + t33 * t15) * t44 + t3 * t50) + t50 * (t1 * t44 + t2 * t45 + t5);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
