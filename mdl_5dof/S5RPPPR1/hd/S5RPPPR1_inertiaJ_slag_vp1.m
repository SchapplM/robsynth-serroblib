% Calculate joint inertia matrix for
% S5RPPPR1
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:47
% DurationCPUTime: 0.40s
% Computational Cost: add. (1236->146), mult. (1181->222), div. (0->0), fcn. (1194->10), ass. (0->71)
t81 = m(5) / 0.2e1;
t80 = m(6) / 0.2e1;
t79 = -m(5) - m(6);
t51 = sin(pkin(9));
t78 = pkin(4) * t51;
t56 = sin(qJ(1));
t77 = t56 * pkin(1);
t57 = cos(qJ(1));
t76 = t57 * pkin(1);
t49 = pkin(9) + qJ(5);
t44 = sin(t49);
t46 = cos(t49);
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t28 = -Icges(6,5) * t54 + (Icges(6,1) * t46 - Icges(6,4) * t44) * t52;
t23 = t52 * t46 * t28;
t26 = -Icges(6,3) * t54 + (Icges(6,5) * t46 - Icges(6,6) * t44) * t52;
t27 = -Icges(6,6) * t54 + (Icges(6,4) * t46 - Icges(6,2) * t44) * t52;
t74 = t44 * t27;
t75 = (-t54 * t26 - t52 * t74 + t23) * t54;
t50 = qJ(1) + pkin(7);
t45 = sin(t50);
t73 = t45 * t52;
t72 = t45 * t54;
t47 = cos(t50);
t71 = t47 * t52;
t70 = t47 * t54;
t30 = t44 * t72 + t47 * t46;
t31 = t47 * t44 - t46 * t72;
t69 = t31 * rSges(6,1) + t30 * rSges(6,2);
t68 = t45 ^ 2 + t47 ^ 2;
t67 = Icges(6,5) * t52;
t66 = Icges(6,6) * t52;
t65 = Icges(6,3) * t52;
t64 = t81 + t80;
t63 = t47 * qJ(3) - t77;
t53 = cos(pkin(9));
t62 = t51 * rSges(5,1) + t53 * rSges(5,2);
t32 = -t44 * t70 + t45 * t46;
t33 = t45 * t44 + t46 * t70;
t61 = -t33 * rSges(6,1) - t32 * rSges(6,2);
t60 = -rSges(4,1) * t54 + rSges(4,2) * t52 - pkin(2);
t59 = -(t53 * pkin(4) + pkin(3)) * t54 - pkin(2) + (-rSges(6,3) - pkin(6) - qJ(4)) * t52;
t58 = -pkin(2) + (-rSges(5,3) - qJ(4)) * t52 + (-rSges(5,1) * t53 + rSges(5,2) * t51 - pkin(3)) * t54;
t39 = -t57 * rSges(2,1) + t56 * rSges(2,2);
t38 = -t56 * rSges(2,1) - t57 * rSges(2,2);
t36 = -t47 * rSges(3,1) + t45 * rSges(3,2) - t76;
t35 = -t45 * rSges(3,1) - t47 * rSges(3,2) - t77;
t29 = -t54 * rSges(6,3) + (rSges(6,1) * t46 - rSges(6,2) * t44) * t52;
t22 = -t76 + (-rSges(4,3) - qJ(3)) * t45 + t60 * t47;
t21 = t47 * rSges(4,3) + t60 * t45 + t63;
t20 = rSges(6,3) * t71 - t61;
t19 = -rSges(6,3) * t73 + t69;
t18 = Icges(6,1) * t33 + Icges(6,4) * t32 + t47 * t67;
t17 = Icges(6,1) * t31 + Icges(6,4) * t30 - t45 * t67;
t16 = Icges(6,4) * t33 + Icges(6,2) * t32 + t47 * t66;
t15 = Icges(6,4) * t31 + Icges(6,2) * t30 - t45 * t66;
t14 = Icges(6,5) * t33 + Icges(6,6) * t32 + t47 * t65;
t13 = Icges(6,5) * t31 + Icges(6,6) * t30 - t45 * t65;
t12 = -t76 + (-qJ(3) - t62) * t45 + t58 * t47;
t11 = t58 * t45 + t62 * t47 + t63;
t9 = t54 * t20 + t29 * t71;
t8 = -t54 * t19 + t29 * t73;
t7 = -t76 + (-qJ(3) - t78) * t45 + t59 * t47 + t61;
t6 = t59 * t45 + t47 * t78 + t63 + t69;
t5 = t26 * t71 + t32 * t27 + t33 * t28;
t4 = -t26 * t73 + t30 * t27 + t31 * t28;
t3 = (-t19 * t47 - t20 * t45) * t52;
t2 = -t54 * t14 + (-t16 * t44 + t18 * t46) * t52;
t1 = -t54 * t13 + (-t15 * t44 + t17 * t46) * t52;
t10 = [Icges(2,3) + Icges(3,3) + t23 + (-t26 + (Icges(4,2) + Icges(5,3)) * t54) * t54 + m(2) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + (-t74 + (Icges(5,1) * t53 ^ 2 + Icges(4,1) + (-0.2e1 * Icges(5,4) * t53 + Icges(5,2) * t51) * t51) * t52 + 0.2e1 * (-Icges(5,5) * t53 + Icges(5,6) * t51 + Icges(4,4)) * t54) * t52; 0; m(3) + m(4) - t79; m(4) * (t45 * t21 + t47 * t22) + m(5) * (t45 * t11 + t47 * t12) + m(6) * (t45 * t6 + t47 * t7); 0; 0.2e1 * (m(4) / 0.2e1 + t64) * t68; 0.2e1 * ((t11 * t47 - t12 * t45) * t81 + (-t45 * t7 + t47 * t6) * t80) * t52; t79 * t54; 0; 0.2e1 * t64 * (t68 * t52 ^ 2 + t54 ^ 2); -t75 + m(6) * (t8 * t6 + t9 * t7) + ((t2 / 0.2e1 + t5 / 0.2e1) * t47 + (-t1 / 0.2e1 - t4 / 0.2e1) * t45) * t52; m(6) * t3; m(6) * (t8 * t45 + t9 * t47); m(6) * (-t3 * t54 + (-t45 * t9 + t47 * t8) * t52); m(6) * (t3 ^ 2 + t8 ^ 2 + t9 ^ 2) - t54 * (-t75 + (-t1 * t45 + t2 * t47) * t52) - (-t4 * t54 - (-t13 * t73 + t30 * t15 + t31 * t17) * t73 + (-t14 * t73 + t30 * t16 + t31 * t18) * t71) * t73 + (-t5 * t54 - (t13 * t71 + t32 * t15 + t33 * t17) * t73 + (t14 * t71 + t32 * t16 + t33 * t18) * t71) * t71;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
