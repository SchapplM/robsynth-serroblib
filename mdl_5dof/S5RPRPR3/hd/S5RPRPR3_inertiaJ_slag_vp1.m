% Calculate joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:16
% DurationCPUTime: 0.43s
% Computational Cost: add. (1697->140), mult. (1430->209), div. (0->0), fcn. (1440->10), ass. (0->72)
t65 = sin(qJ(1));
t87 = t65 * pkin(1);
t67 = cos(qJ(1));
t86 = t67 * pkin(1);
t62 = sin(pkin(9));
t63 = cos(pkin(9));
t64 = sin(qJ(5));
t66 = cos(qJ(5));
t38 = -Icges(6,3) * t63 + (Icges(6,5) * t66 - Icges(6,6) * t64) * t62;
t39 = -Icges(6,6) * t63 + (Icges(6,4) * t66 - Icges(6,2) * t64) * t62;
t40 = -Icges(6,5) * t63 + (Icges(6,1) * t66 - Icges(6,4) * t64) * t62;
t16 = -t63 * t38 + (-t39 * t64 + t40 * t66) * t62;
t85 = t16 * t63;
t61 = qJ(1) + pkin(8);
t60 = qJ(3) + t61;
t56 = sin(t60);
t84 = t56 * t62;
t57 = cos(t60);
t83 = t57 * t62;
t82 = t63 * t64;
t81 = t63 * t66;
t34 = t56 * t82 + t57 * t66;
t35 = -t56 * t81 + t57 * t64;
t80 = t35 * rSges(6,1) + t34 * rSges(6,2);
t79 = Icges(6,5) * t62;
t78 = Icges(6,6) * t62;
t77 = Icges(6,3) * t62;
t74 = -rSges(5,1) * t63 - pkin(3);
t45 = -t57 * rSges(4,1) + t56 * rSges(4,2);
t58 = sin(t61);
t73 = -pkin(2) * t58 - t87;
t59 = cos(t61);
t72 = -pkin(2) * t59 - t86;
t44 = -t56 * rSges(4,1) - t57 * rSges(4,2);
t36 = t56 * t66 - t57 * t82;
t37 = t56 * t64 + t57 * t81;
t71 = -t37 * rSges(6,1) - t36 * rSges(6,2);
t17 = Icges(6,5) * t35 + Icges(6,6) * t34 - t56 * t77;
t19 = Icges(6,4) * t35 + Icges(6,2) * t34 - t56 * t78;
t21 = Icges(6,1) * t35 + Icges(6,4) * t34 - t56 * t79;
t3 = -t63 * t17 + (-t19 * t64 + t21 * t66) * t62;
t18 = Icges(6,5) * t37 + Icges(6,6) * t36 + t57 * t77;
t20 = Icges(6,4) * t37 + Icges(6,2) * t36 + t57 * t78;
t22 = Icges(6,1) * t37 + Icges(6,4) * t36 + t57 * t79;
t4 = -t63 * t18 + (-t20 * t64 + t22 * t66) * t62;
t8 = t34 * t39 + t35 * t40 - t38 * t84;
t9 = t36 * t39 + t37 * t40 + t38 * t83;
t70 = -t85 - (t3 + t8) * t84 / 0.2e1 + (t4 + t9) * t83 / 0.2e1;
t69 = -pkin(4) * t63 - pkin(3) + (-rSges(6,3) - pkin(7)) * t62;
t53 = t57 * qJ(4);
t27 = rSges(5,2) * t84 + t57 * rSges(5,3) + t74 * t56 + t53;
t68 = Icges(5,2) * t63 ^ 2 + Icges(4,3) + t16 + (Icges(5,1) * t62 + 0.2e1 * Icges(5,4) * t63) * t62;
t28 = rSges(5,2) * t83 + t74 * t57 + (-rSges(5,3) - qJ(4)) * t56;
t12 = t69 * t56 + t53 + t80;
t13 = -t56 * qJ(4) + t69 * t57 + t71;
t52 = -t67 * rSges(2,1) + t65 * rSges(2,2);
t51 = -t65 * rSges(2,1) - t67 * rSges(2,2);
t43 = -t59 * rSges(3,1) + t58 * rSges(3,2) - t86;
t42 = -t58 * rSges(3,1) - t59 * rSges(3,2) - t87;
t41 = -t63 * rSges(6,3) + (rSges(6,1) * t66 - rSges(6,2) * t64) * t62;
t30 = t45 + t72;
t29 = t44 + t73;
t26 = t28 + t72;
t25 = t27 + t73;
t24 = rSges(6,3) * t83 - t71;
t23 = -rSges(6,3) * t84 + t80;
t15 = t63 * t24 + t41 * t83;
t14 = -t63 * t23 + t41 * t84;
t11 = t13 + t72;
t10 = t12 + t73;
t5 = (-t23 * t57 - t24 * t56) * t62;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t68; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t44 * t29 + t45 * t30) + m(5) * (t27 * t25 + t28 * t26) + m(6) * (t12 * t10 + t13 * t11) + t68; 0; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t68; m(5) * (t56 * t25 + t57 * t26) + m(6) * (t56 * t10 + t57 * t11); 0; m(6) * (t56 * t12 + t57 * t13) + m(5) * (t56 * t27 + t57 * t28); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t56 ^ 2 + t57 ^ 2); m(6) * (t14 * t10 + t15 * t11) + t70; m(6) * t5; m(6) * (t14 * t12 + t15 * t13) + t70; m(6) * (t14 * t56 + t15 * t57); m(6) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) - t63 * (-t85 + (-t3 * t56 + t4 * t57) * t62) - (-t8 * t63 - (-t17 * t84 + t34 * t19 + t35 * t21) * t84 + (-t18 * t84 + t34 * t20 + t35 * t22) * t83) * t84 + (-t9 * t63 - (t17 * t83 + t36 * t19 + t37 * t21) * t84 + (t18 * t83 + t36 * t20 + t37 * t22) * t83) * t83;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
