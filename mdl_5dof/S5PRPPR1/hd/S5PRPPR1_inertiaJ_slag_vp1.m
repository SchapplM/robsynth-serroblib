% Calculate joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:21:48
% DurationCPUTime: 0.44s
% Computational Cost: add. (1212->138), mult. (1151->216), div. (0->0), fcn. (1170->8), ass. (0->66)
t52 = sin(pkin(9));
t53 = sin(pkin(8));
t54 = cos(pkin(9));
t55 = cos(pkin(8));
t79 = (rSges(5,3) + qJ(4)) * t53 + (rSges(5,1) * t54 - rSges(5,2) * t52 + pkin(3)) * t55;
t78 = m(5) / 0.2e1;
t77 = m(6) / 0.2e1;
t76 = -m(5) - m(6);
t75 = pkin(4) * t52;
t50 = pkin(9) + qJ(5);
t45 = sin(t50);
t47 = cos(t50);
t28 = -Icges(6,5) * t55 + (Icges(6,1) * t47 - Icges(6,4) * t45) * t53;
t23 = t53 * t47 * t28;
t26 = -Icges(6,3) * t55 + (Icges(6,5) * t47 - Icges(6,6) * t45) * t53;
t27 = -Icges(6,6) * t55 + (Icges(6,4) * t47 - Icges(6,2) * t45) * t53;
t72 = t45 * t27;
t74 = (-t55 * t26 - t53 * t72 + t23) * t55;
t73 = (t54 * pkin(4) + pkin(3)) * t55;
t51 = pkin(7) + qJ(2);
t46 = sin(t51);
t71 = t46 * t53;
t70 = t46 * t55;
t48 = cos(t51);
t69 = t48 * t53;
t68 = t48 * t55;
t67 = t48 * pkin(2) + t46 * qJ(3);
t66 = t46 ^ 2 + t48 ^ 2;
t65 = Icges(6,5) * t53;
t64 = Icges(6,6) * t53;
t63 = Icges(6,3) * t53;
t61 = t78 + t77;
t32 = -t45 * t68 + t46 * t47;
t33 = t46 * t45 + t47 * t68;
t20 = t33 * rSges(6,1) + t32 * rSges(6,2) + rSges(6,3) * t69;
t60 = rSges(4,1) * t55 - rSges(4,2) * t53;
t59 = t52 * rSges(5,1) + t54 * rSges(5,2);
t30 = -t45 * t70 - t48 * t47;
t31 = -t48 * t45 + t47 * t70;
t58 = -t31 * rSges(6,1) - t30 * rSges(6,2);
t56 = -pkin(6) - qJ(4);
t40 = t48 * qJ(3);
t36 = t48 * rSges(3,1) - t46 * rSges(3,2);
t35 = -t46 * rSges(3,1) - t48 * rSges(3,2);
t29 = -t55 * rSges(6,3) + (rSges(6,1) * t47 - rSges(6,2) * t45) * t53;
t22 = t46 * rSges(4,3) + t60 * t48 + t67;
t21 = t48 * rSges(4,3) + t40 + (-pkin(2) - t60) * t46;
t19 = rSges(6,3) * t71 - t58;
t18 = Icges(6,1) * t33 + Icges(6,4) * t32 + t48 * t65;
t17 = Icges(6,1) * t31 + Icges(6,4) * t30 + t46 * t65;
t16 = Icges(6,4) * t33 + Icges(6,2) * t32 + t48 * t64;
t15 = Icges(6,4) * t31 + Icges(6,2) * t30 + t46 * t64;
t14 = Icges(6,5) * t33 + Icges(6,6) * t32 + t48 * t63;
t13 = Icges(6,5) * t31 + Icges(6,6) * t30 + t46 * t63;
t12 = t59 * t46 + t79 * t48 + t67;
t11 = t40 + t59 * t48 + (-pkin(2) - t79) * t46;
t9 = -t55 * t20 - t29 * t69;
t8 = t55 * t19 + t29 * t71;
t7 = t46 * t75 + (-t53 * t56 + t73) * t48 + t20 + t67;
t6 = t48 * t75 + t40 + (-t73 - pkin(2) + (-rSges(6,3) + t56) * t53) * t46 + t58;
t5 = t26 * t69 + t32 * t27 + t33 * t28;
t4 = t26 * t71 + t30 * t27 + t31 * t28;
t3 = (t19 * t48 - t20 * t46) * t53;
t2 = -t55 * t14 + (-t16 * t45 + t18 * t47) * t53;
t1 = -t55 * t13 + (-t15 * t45 + t17 * t47) * t53;
t10 = [m(2) + m(3) + m(4) - t76; 0; Icges(3,3) + t23 + (-t26 + (Icges(4,2) + Icges(5,3)) * t55) * t55 + m(6) * (t6 ^ 2 + t7 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t21 ^ 2 + t22 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + (-t72 + (Icges(5,1) * t54 ^ 2 + Icges(4,1) + (-0.2e1 * Icges(5,4) * t54 + Icges(5,2) * t52) * t52) * t53 + 0.2e1 * (-Icges(5,5) * t54 + Icges(5,6) * t52 + Icges(4,4)) * t55) * t53; 0; m(6) * (t46 * t6 - t48 * t7) + m(5) * (t46 * t11 - t48 * t12) + m(4) * (t46 * t21 - t48 * t22); 0.2e1 * (m(4) / 0.2e1 + t61) * t66; t76 * t55; 0.2e1 * ((t46 * t7 + t48 * t6) * t77 + (t11 * t48 + t12 * t46) * t78) * t53; 0; 0.2e1 * t61 * (t66 * t53 ^ 2 + t55 ^ 2); m(6) * t3; -t74 + m(6) * (t8 * t6 + t9 * t7) + ((t2 / 0.2e1 + t5 / 0.2e1) * t48 + (t1 / 0.2e1 + t4 / 0.2e1) * t46) * t53; m(6) * (t8 * t46 - t9 * t48); m(6) * (-t3 * t55 + (t46 * t9 + t48 * t8) * t53); m(6) * (t3 ^ 2 + t8 ^ 2 + t9 ^ 2) + ((t14 * t69 + t32 * t16 + t33 * t18) * t69 + (t13 * t69 + t32 * t15 + t33 * t17) * t71 - t5 * t55) * t69 + ((t14 * t71 + t30 * t16 + t31 * t18) * t69 + (t13 * t71 + t30 * t15 + t31 * t17) * t71 - t4 * t55) * t71 - t55 * (-t74 + (t1 * t46 + t2 * t48) * t53);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
