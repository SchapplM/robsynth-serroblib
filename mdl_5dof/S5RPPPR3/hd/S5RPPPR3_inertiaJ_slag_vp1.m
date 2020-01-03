% Calculate joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (713->116), mult. (868->171), div. (0->0), fcn. (867->8), ass. (0->54)
t46 = cos(pkin(8));
t69 = t46 ^ 2;
t68 = m(5) / 0.2e1;
t65 = -m(5) - m(6);
t48 = sin(qJ(1));
t64 = t48 * pkin(1);
t63 = -rSges(6,3) - pkin(6);
t44 = qJ(1) + pkin(7);
t41 = cos(t44);
t62 = t41 * t46;
t45 = sin(pkin(8));
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t29 = t45 * t49 - t46 * t47;
t23 = t29 * t41;
t52 = t45 * t47 + t46 * t49;
t24 = t52 * t41;
t61 = t24 * rSges(6,1) + t23 * rSges(6,2);
t40 = sin(t44);
t60 = t40 ^ 2 + t41 ^ 2;
t59 = qJ(4) * t45;
t58 = t68 + m(6) / 0.2e1;
t50 = cos(qJ(1));
t42 = t50 * pkin(1);
t57 = t41 * pkin(2) + t40 * qJ(3) + t42;
t56 = t41 * qJ(3) - t64;
t55 = pkin(3) * t62 + t41 * t59 + t57;
t54 = rSges(4,1) * t46 - rSges(4,2) * t45;
t21 = t29 * t40;
t22 = t52 * t40;
t53 = t22 * rSges(6,1) + t21 * rSges(6,2);
t2 = t63 * t41 + (-t59 - pkin(2) + (-pkin(3) - pkin(4)) * t46) * t40 - t53 + t56;
t3 = pkin(4) * t62 + t63 * t40 + t55 + t61;
t51 = m(6) * (t2 * t41 + t3 * t40);
t33 = t50 * rSges(2,1) - t48 * rSges(2,2);
t32 = -t48 * rSges(2,1) - t50 * rSges(2,2);
t27 = t41 * rSges(3,1) - t40 * rSges(3,2) + t42;
t26 = -t40 * rSges(3,1) - t41 * rSges(3,2) - t64;
t18 = t29 * rSges(6,1) - rSges(6,2) * t52;
t17 = Icges(6,1) * t29 - Icges(6,4) * t52;
t16 = Icges(6,4) * t29 - Icges(6,2) * t52;
t15 = Icges(6,5) * t29 - Icges(6,6) * t52;
t13 = t40 * rSges(4,3) + t54 * t41 + t57;
t12 = t41 * rSges(4,3) + (-pkin(2) - t54) * t40 + t56;
t11 = Icges(6,1) * t24 + Icges(6,4) * t23 - Icges(6,5) * t40;
t10 = Icges(6,1) * t22 + Icges(6,4) * t21 + Icges(6,5) * t41;
t9 = Icges(6,4) * t24 + Icges(6,2) * t23 - Icges(6,6) * t40;
t8 = Icges(6,4) * t22 + Icges(6,2) * t21 + Icges(6,6) * t41;
t7 = Icges(6,5) * t24 + Icges(6,6) * t23 - Icges(6,3) * t40;
t6 = Icges(6,5) * t22 + Icges(6,6) * t21 + Icges(6,3) * t41;
t5 = t40 * rSges(5,2) + (rSges(5,1) * t46 + rSges(5,3) * t45) * t41 + t55;
t4 = t41 * rSges(5,2) + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t46 + (-rSges(5,3) - qJ(4)) * t45) * t40 + t56;
t1 = -t40 * t53 - t41 * t61;
t14 = [-t52 * t16 + t29 * t17 + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t69 + ((Icges(4,1) + Icges(5,1)) * t45 + 0.2e1 * (Icges(4,4) - Icges(5,5)) * t46) * t45 + m(2) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2); 0; m(3) + m(4) - t65; m(4) * (t40 * t12 - t41 * t13) + m(5) * (t40 * t4 - t41 * t5) + m(6) * (t40 * t2 - t41 * t3); 0; 0.2e1 * (m(4) / 0.2e1 + t58) * t60; 0.2e1 * ((t4 * t41 + t40 * t5) * t68 + t51 / 0.2e1) * t45; t65 * t46; 0; 0.2e1 * t58 * (t60 * t45 ^ 2 + t69); t18 * t51 - (t29 * t11 - t40 * t15 + t23 * t16 + t24 * t17 - t52 * t9) * t40 / 0.2e1 + (t29 * t10 + t41 * t15 + t21 * t16 + t22 * t17 - t52 * t8) * t41 / 0.2e1; m(6) * t1; 0; m(6) * (t60 * t45 * t18 - t1 * t46); m(6) * (t60 * t18 ^ 2 + t1 ^ 2) - t40 * (-(t24 * t11 + t23 * t9 - t40 * t7) * t40 + (t24 * t10 + t23 * t8 - t40 * t6) * t41) + t41 * (-(t22 * t11 + t21 * t9 + t41 * t7) * t40 + (t22 * t10 + t21 * t8 + t41 * t6) * t41);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
