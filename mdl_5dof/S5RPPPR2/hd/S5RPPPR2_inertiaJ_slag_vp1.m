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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:02
% DurationCPUTime: 0.59s
% Computational Cost: add. (1496->172), mult. (3769->274), div. (0->0), fcn. (4609->10), ass. (0->86)
t72 = cos(pkin(7));
t103 = 0.2e1 * t72;
t97 = m(6) / 0.2e1;
t98 = m(5) / 0.2e1;
t86 = t98 + t97;
t102 = 0.2e1 * t86;
t101 = t72 ^ 2;
t69 = sin(pkin(7));
t100 = 0.2e1 * t69;
t99 = m(4) / 0.2e1;
t68 = sin(pkin(8));
t96 = t68 * t69;
t71 = cos(pkin(8));
t95 = t69 * t71;
t74 = sin(qJ(1));
t94 = t69 * t74;
t76 = cos(qJ(1));
t93 = t69 * t76;
t92 = t74 * t68;
t91 = t74 * t71;
t90 = t76 * t68;
t89 = t76 * t71;
t88 = t74 ^ 2 + t76 ^ 2;
t87 = t74 * qJ(2);
t67 = sin(pkin(9));
t70 = cos(pkin(9));
t52 = -t72 * t67 + t70 * t95;
t73 = sin(qJ(5));
t75 = cos(qJ(5));
t40 = -t52 * t73 + t75 * t96;
t41 = t52 * t75 + t73 * t96;
t51 = t67 * t95 + t72 * t70;
t24 = Icges(6,5) * t41 + Icges(6,6) * t40 + Icges(6,3) * t51;
t25 = Icges(6,4) * t41 + Icges(6,2) * t40 + Icges(6,6) * t51;
t26 = Icges(6,1) * t41 + Icges(6,4) * t40 + Icges(6,5) * t51;
t85 = t51 * t24 + t40 * t25 + t41 * t26;
t55 = -t72 * t91 + t90;
t43 = t55 * t70 - t67 * t94;
t53 = t72 * t92 + t89;
t33 = -t43 * t73 - t53 * t75;
t34 = t43 * t75 - t53 * t73;
t42 = t55 * t67 + t70 * t94;
t17 = t34 * rSges(6,1) + t33 * rSges(6,2) + t42 * rSges(6,3);
t84 = -pkin(2) * t72 - pkin(1);
t83 = t99 + t86;
t57 = t72 * t89 + t92;
t45 = t57 * t70 + t67 * t93;
t56 = t72 * t90 - t91;
t35 = -t45 * t73 + t56 * t75;
t36 = t45 * t75 + t56 * t73;
t82 = -t36 * rSges(6,1) - t35 * rSges(6,2);
t81 = -rSges(3,1) * t72 + rSges(3,2) * t69 - pkin(1);
t80 = -qJ(3) * t69 + t84;
t79 = (-rSges(4,3) - qJ(3)) * t69 + t84;
t63 = t76 * qJ(2);
t78 = t55 * pkin(3) + t80 * t74 + t63;
t77 = -t57 * pkin(3) - t56 * qJ(4) + t80 * t76 - t87;
t64 = t69 ^ 2;
t59 = -t76 * rSges(2,1) + t74 * rSges(2,2);
t58 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t47 = (-rSges(3,3) - qJ(2)) * t74 + t81 * t76;
t46 = t76 * rSges(3,3) + t81 * t74 + t63;
t44 = t57 * t67 - t70 * t93;
t29 = -t57 * rSges(4,1) + t56 * rSges(4,2) + t79 * t76 - t87;
t28 = t55 * rSges(4,1) + t53 * rSges(4,2) + t79 * t74 + t63;
t27 = t41 * rSges(6,1) + t40 * rSges(6,2) + t51 * rSges(6,3);
t20 = -t45 * rSges(5,1) + t44 * rSges(5,2) - t56 * rSges(5,3) + t77;
t19 = t43 * rSges(5,1) - t42 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t53 + t78;
t18 = t44 * rSges(6,3) - t82;
t16 = Icges(6,1) * t36 + Icges(6,4) * t35 + Icges(6,5) * t44;
t15 = Icges(6,1) * t34 + Icges(6,4) * t33 + Icges(6,5) * t42;
t14 = Icges(6,4) * t36 + Icges(6,2) * t35 + Icges(6,6) * t44;
t13 = Icges(6,4) * t34 + Icges(6,2) * t33 + Icges(6,6) * t42;
t12 = Icges(6,5) * t36 + Icges(6,6) * t35 + Icges(6,3) * t44;
t11 = Icges(6,5) * t34 + Icges(6,6) * t33 + Icges(6,3) * t42;
t10 = -t45 * pkin(4) + (-rSges(6,3) - pkin(6)) * t44 + t77 + t82;
t9 = t43 * pkin(4) + t42 * pkin(6) - t53 * qJ(4) + t17 + t78;
t8 = -t51 * t18 + t44 * t27;
t7 = t51 * t17 - t42 * t27;
t6 = -t44 * t17 + t42 * t18;
t5 = t85 * t51;
t4 = t44 * t24 + t35 * t25 + t36 * t26;
t3 = t42 * t24 + t33 * t25 + t34 * t26;
t2 = t51 * t12 + t40 * t14 + t41 * t16;
t1 = t51 * t11 + t40 * t13 + t41 * t15;
t21 = [Icges(5,1) * t52 ^ 2 + Icges(2,3) + (-0.2e1 * Icges(5,4) * t52 + Icges(5,2) * t51) * t51 + (Icges(3,2) + Icges(4,3)) * t101 + m(2) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + m(6) * (t10 ^ 2 + t9 ^ 2) + ((Icges(4,1) * t71 ^ 2 + Icges(3,1)) * t69 + (-Icges(4,5) * t71 + Icges(3,4)) * t103 + (0.2e1 * Icges(5,5) * t52 + Icges(4,6) * t103 - 0.2e1 * Icges(5,6) * t51 + (-0.2e1 * Icges(4,4) * t71 + (Icges(4,2) + Icges(5,3)) * t68) * t69) * t68) * t69 + t85; m(3) * (t74 * t46 + t76 * t47) + m(4) * (t74 * t28 + t76 * t29) + m(5) * (t74 * t19 + t76 * t20) + m(6) * (t76 * t10 + t74 * t9); 0.2e1 * (m(3) / 0.2e1 + t83) * t88; ((t28 * t76 - t29 * t74) * t99 + (t19 * t76 - t20 * t74) * t98 + (-t10 * t74 + t76 * t9) * t97) * t100; 0; 0.2e1 * t83 * (t88 * t64 + t101); m(5) * (t56 * t19 - t53 * t20) + m(6) * (-t53 * t10 + t56 * t9); (-t53 * t76 + t56 * t74) * t102; t86 * (t53 * t74 + t56 * t76 - t72 * t68) * t100; (t64 * t68 ^ 2 + t53 ^ 2 + t56 ^ 2) * t102; t5 + m(6) * (t8 * t10 + t7 * t9) + (t2 / 0.2e1 + t4 / 0.2e1) * t44 + (t1 / 0.2e1 + t3 / 0.2e1) * t42; m(6) * (t7 * t74 + t8 * t76); m(6) * (-t6 * t72 + (t7 * t76 - t74 * t8) * t69); m(6) * (-t8 * t53 + t7 * t56 + t6 * t96); m(6) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) + t51 * (t1 * t42 + t2 * t44 + t5) + t42 * (t3 * t51 + (t42 * t11 + t33 * t13 + t34 * t15) * t42 + (t42 * t12 + t33 * t14 + t34 * t16) * t44) + t44 * (t4 * t51 + (t44 * t11 + t35 * t13 + t36 * t15) * t42 + (t44 * t12 + t35 * t14 + t36 * t16) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
