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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:06
% EndTime: 2020-01-03 11:20:07
% DurationCPUTime: 0.42s
% Computational Cost: add. (1236->146), mult. (1181->225), div. (0->0), fcn. (1194->10), ass. (0->72)
t84 = m(5) / 0.2e1;
t83 = m(6) / 0.2e1;
t82 = -m(5) - m(6);
t56 = sin(pkin(9));
t81 = pkin(4) * t56;
t54 = pkin(9) + qJ(5);
t47 = sin(t54);
t49 = cos(t54);
t57 = sin(pkin(8));
t59 = cos(pkin(8));
t28 = -Icges(6,5) * t59 + (Icges(6,1) * t49 - Icges(6,4) * t47) * t57;
t23 = t57 * t49 * t28;
t26 = -Icges(6,3) * t59 + (Icges(6,5) * t49 - Icges(6,6) * t47) * t57;
t27 = -Icges(6,6) * t59 + (Icges(6,4) * t49 - Icges(6,2) * t47) * t57;
t78 = t47 * t27;
t80 = (-t59 * t26 - t57 * t78 + t23) * t59;
t58 = cos(pkin(9));
t79 = (t58 * pkin(4) + pkin(3)) * t59;
t55 = qJ(1) + pkin(7);
t48 = sin(t55);
t77 = t48 * t57;
t76 = t48 * t59;
t50 = cos(t55);
t75 = t50 * t57;
t74 = t50 * t59;
t61 = sin(qJ(1));
t51 = t61 * pkin(1);
t73 = t48 * pkin(2) + t51;
t72 = t48 ^ 2 + t50 ^ 2;
t71 = Icges(6,5) * t57;
t70 = Icges(6,6) * t57;
t69 = Icges(6,3) * t57;
t68 = t84 + t83;
t30 = -t47 * t76 - t50 * t49;
t31 = -t50 * t47 + t49 * t76;
t19 = t31 * rSges(6,1) + t30 * rSges(6,2) + rSges(6,3) * t77;
t62 = cos(qJ(1));
t52 = t62 * pkin(1);
t67 = t50 * pkin(2) + t48 * qJ(3) + t52;
t66 = rSges(4,1) * t59 - rSges(4,2) * t57;
t65 = t56 * rSges(5,1) + t58 * rSges(5,2);
t32 = t47 * t74 - t48 * t49;
t33 = -t48 * t47 - t49 * t74;
t64 = -t33 * rSges(6,1) - t32 * rSges(6,2);
t63 = (rSges(5,3) + qJ(4)) * t57 + (rSges(5,1) * t58 - rSges(5,2) * t56 + pkin(3)) * t59;
t60 = -pkin(6) - qJ(4);
t40 = t62 * rSges(2,1) - t61 * rSges(2,2);
t39 = t61 * rSges(2,1) + t62 * rSges(2,2);
t36 = t50 * rSges(3,1) - t48 * rSges(3,2) + t52;
t35 = t48 * rSges(3,1) + t50 * rSges(3,2) + t51;
t29 = -t59 * rSges(6,3) + (rSges(6,1) * t49 - rSges(6,2) * t47) * t57;
t22 = t48 * rSges(4,3) + t66 * t50 + t67;
t21 = (-rSges(4,3) - qJ(3)) * t50 + t66 * t48 + t73;
t20 = -rSges(6,3) * t75 - t64;
t18 = Icges(6,1) * t33 + Icges(6,4) * t32 - t50 * t71;
t17 = Icges(6,1) * t31 + Icges(6,4) * t30 + t48 * t71;
t16 = Icges(6,4) * t33 + Icges(6,2) * t32 - t50 * t70;
t15 = Icges(6,4) * t31 + Icges(6,2) * t30 + t48 * t70;
t14 = Icges(6,5) * t33 + Icges(6,6) * t32 - t50 * t69;
t13 = Icges(6,5) * t31 + Icges(6,6) * t30 + t48 * t69;
t12 = t65 * t48 + t63 * t50 + t67;
t11 = (-qJ(3) - t65) * t50 + t63 * t48 + t73;
t9 = t59 * t20 - t29 * t75;
t8 = -t59 * t19 - t29 * t77;
t7 = t48 * t81 + (t79 + (rSges(6,3) - t60) * t57) * t50 + t64 + t67;
t6 = (-qJ(3) - t81) * t50 + (-t57 * t60 + t79) * t48 + t19 + t73;
t5 = -t26 * t75 + t32 * t27 + t33 * t28;
t4 = t26 * t77 + t30 * t27 + t31 * t28;
t3 = (t19 * t50 + t20 * t48) * t57;
t2 = -t59 * t14 + (-t16 * t47 + t18 * t49) * t57;
t1 = -t59 * t13 + (-t15 * t47 + t17 * t49) * t57;
t10 = [Icges(2,3) + Icges(3,3) + t23 + (-t26 + (Icges(4,2) + Icges(5,3)) * t59) * t59 + m(2) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + (-t78 + (Icges(5,1) * t58 ^ 2 + Icges(4,1) + (-0.2e1 * Icges(5,4) * t58 + Icges(5,2) * t56) * t56) * t57 + 0.2e1 * (-Icges(5,5) * t58 + Icges(5,6) * t56 + Icges(4,4)) * t59) * t57; 0; m(3) + m(4) - t82; m(4) * (-t48 * t21 - t50 * t22) + m(5) * (-t48 * t11 - t50 * t12) + m(6) * (-t48 * t6 - t50 * t7); 0; 0.2e1 * (m(4) / 0.2e1 + t68) * t72; 0.2e1 * ((-t11 * t50 + t12 * t48) * t84 + (t48 * t7 - t50 * t6) * t83) * t57; t82 * t59; 0; 0.2e1 * t68 * (t72 * t57 ^ 2 + t59 ^ 2); -t80 + m(6) * (t8 * t6 + t9 * t7) + ((-t2 / 0.2e1 - t5 / 0.2e1) * t50 + (t1 / 0.2e1 + t4 / 0.2e1) * t48) * t57; m(6) * t3; m(6) * (-t8 * t48 - t9 * t50); m(6) * (-t3 * t59 + (t48 * t9 - t50 * t8) * t57); m(6) * (t3 ^ 2 + t8 ^ 2 + t9 ^ 2) - t59 * (-t80 + (t1 * t48 - t2 * t50) * t57) + (-t4 * t59 + (t13 * t77 + t30 * t15 + t31 * t17) * t77 - (t14 * t77 + t30 * t16 + t31 * t18) * t75) * t77 - (-t5 * t59 + (-t13 * t75 + t32 * t15 + t33 * t17) * t77 - (-t14 * t75 + t32 * t16 + t33 * t18) * t75) * t75;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
