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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:21
% EndTime: 2020-01-03 11:22:23
% DurationCPUTime: 0.65s
% Computational Cost: add. (1496->171), mult. (3769->275), div. (0->0), fcn. (4609->10), ass. (0->86)
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t111 = pkin(2) * t78 + qJ(3) * t75;
t110 = 0.2e1 * t78;
t104 = m(6) / 0.2e1;
t105 = m(5) / 0.2e1;
t91 = t105 + t104;
t109 = 0.2e1 * t91;
t108 = t78 ^ 2;
t107 = 0.2e1 * t75;
t106 = m(4) / 0.2e1;
t74 = sin(pkin(8));
t102 = t74 * t75;
t77 = cos(pkin(8));
t101 = t75 * t77;
t80 = sin(qJ(1));
t100 = t75 * t80;
t82 = cos(qJ(1));
t99 = t75 * t82;
t98 = t80 * t74;
t97 = t80 * t77;
t96 = t82 * t74;
t95 = t82 * t77;
t94 = t82 * pkin(1) + t80 * qJ(2);
t93 = t80 ^ 2 + t82 ^ 2;
t73 = sin(pkin(9));
t76 = cos(pkin(9));
t52 = t76 * t101 - t78 * t73;
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t40 = t81 * t102 - t52 * t79;
t41 = t79 * t102 + t52 * t81;
t51 = t73 * t101 + t78 * t76;
t24 = Icges(6,5) * t41 + Icges(6,6) * t40 + Icges(6,3) * t51;
t25 = Icges(6,4) * t41 + Icges(6,2) * t40 + Icges(6,6) * t51;
t26 = Icges(6,1) * t41 + Icges(6,4) * t40 + Icges(6,5) * t51;
t90 = t51 * t24 + t40 * t25 + t41 * t26;
t54 = t78 * t97 - t96;
t43 = t73 * t100 + t54 * t76;
t53 = t78 * t98 + t95;
t33 = -t43 * t79 + t53 * t81;
t34 = t43 * t81 + t53 * t79;
t42 = -t76 * t100 + t54 * t73;
t17 = t34 * rSges(6,1) + t33 * rSges(6,2) + t42 * rSges(6,3);
t89 = t111 * t82 + t94;
t88 = t106 + t91;
t87 = rSges(3,1) * t78 - rSges(3,2) * t75;
t57 = -t78 * t95 - t98;
t45 = t57 * t76 - t73 * t99;
t55 = t78 * t96 - t97;
t35 = -t45 * t79 - t55 * t81;
t36 = t45 * t81 - t55 * t79;
t86 = -t36 * rSges(6,1) - t35 * rSges(6,2);
t68 = t80 * pkin(1);
t85 = -t82 * qJ(2) + t111 * t80 + t68;
t84 = t54 * pkin(3) + t85;
t83 = -t57 * pkin(3) + t55 * qJ(4) + t89;
t70 = t75 ^ 2;
t59 = t82 * rSges(2,1) - t80 * rSges(2,2);
t58 = t80 * rSges(2,1) + t82 * rSges(2,2);
t47 = t80 * rSges(3,3) + t82 * t87 + t94;
t46 = t68 + (-rSges(3,3) - qJ(2)) * t82 + t87 * t80;
t44 = t57 * t73 + t76 * t99;
t29 = -t57 * rSges(4,1) - t55 * rSges(4,2) + rSges(4,3) * t99 + t89;
t28 = t54 * rSges(4,1) - t53 * rSges(4,2) + rSges(4,3) * t100 + t85;
t27 = t41 * rSges(6,1) + t40 * rSges(6,2) + t51 * rSges(6,3);
t20 = -t45 * rSges(5,1) + t44 * rSges(5,2) + t55 * rSges(5,3) + t83;
t19 = t43 * rSges(5,1) - t42 * rSges(5,2) + (rSges(5,3) + qJ(4)) * t53 + t84;
t18 = t44 * rSges(6,3) - t86;
t16 = Icges(6,1) * t36 + Icges(6,4) * t35 + Icges(6,5) * t44;
t15 = Icges(6,1) * t34 + Icges(6,4) * t33 + Icges(6,5) * t42;
t14 = Icges(6,4) * t36 + Icges(6,2) * t35 + Icges(6,6) * t44;
t13 = Icges(6,4) * t34 + Icges(6,2) * t33 + Icges(6,6) * t42;
t12 = Icges(6,5) * t36 + Icges(6,6) * t35 + Icges(6,3) * t44;
t11 = Icges(6,5) * t34 + Icges(6,6) * t33 + Icges(6,3) * t42;
t10 = -t45 * pkin(4) + (-rSges(6,3) - pkin(6)) * t44 + t83 + t86;
t9 = t43 * pkin(4) + t42 * pkin(6) + t53 * qJ(4) + t17 + t84;
t8 = -t51 * t18 + t44 * t27;
t7 = t51 * t17 - t42 * t27;
t6 = -t44 * t17 + t42 * t18;
t5 = t90 * t51;
t4 = t44 * t24 + t35 * t25 + t36 * t26;
t3 = t42 * t24 + t33 * t25 + t34 * t26;
t2 = t51 * t12 + t40 * t14 + t41 * t16;
t1 = t51 * t11 + t40 * t13 + t41 * t15;
t21 = [Icges(5,1) * t52 ^ 2 + Icges(2,3) + (-0.2e1 * Icges(5,4) * t52 + Icges(5,2) * t51) * t51 + (Icges(3,2) + Icges(4,3)) * t108 + m(2) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + m(6) * (t10 ^ 2 + t9 ^ 2) + ((Icges(4,1) * t77 ^ 2 + Icges(3,1)) * t75 + (-Icges(4,5) * t77 + Icges(3,4)) * t110 + (0.2e1 * Icges(5,5) * t52 + Icges(4,6) * t110 - 0.2e1 * Icges(5,6) * t51 + (-0.2e1 * Icges(4,4) * t77 + (Icges(4,2) + Icges(5,3)) * t74) * t75) * t74) * t75 + t90; m(3) * (-t80 * t46 - t82 * t47) + m(4) * (-t80 * t28 - t82 * t29) + m(5) * (-t80 * t19 - t82 * t20) + m(6) * (-t82 * t10 - t80 * t9); 0.2e1 * (m(3) / 0.2e1 + t88) * t93; ((-t28 * t82 + t29 * t80) * t106 + (-t19 * t82 + t20 * t80) * t105 + (t10 * t80 - t82 * t9) * t104) * t107; 0; 0.2e1 * t88 * (t93 * t70 + t108); m(5) * (-t55 * t19 + t53 * t20) + m(6) * (t53 * t10 - t55 * t9); (-t53 * t82 + t55 * t80) * t109; t91 * (t53 * t80 + t55 * t82 - t78 * t74) * t107; (t70 * t74 ^ 2 + t53 ^ 2 + t55 ^ 2) * t109; t5 + m(6) * (t8 * t10 + t7 * t9) + (t2 / 0.2e1 + t4 / 0.2e1) * t44 + (t1 / 0.2e1 + t3 / 0.2e1) * t42; m(6) * (-t7 * t80 - t8 * t82); m(6) * (-t6 * t78 + (-t7 * t82 + t8 * t80) * t75); m(6) * (t6 * t102 + t8 * t53 - t7 * t55); m(6) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) + t51 * (t1 * t42 + t2 * t44 + t5) + t42 * (t3 * t51 + (t42 * t11 + t33 * t13 + t34 * t15) * t42 + (t42 * t12 + t33 * t14 + t34 * t16) * t44) + t44 * (t4 * t51 + (t44 * t11 + t35 * t13 + t36 * t15) * t42 + (t44 * t12 + t35 * t14 + t36 * t16) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
