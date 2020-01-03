% Calculate joint inertia matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:47
% DurationCPUTime: 0.30s
% Computational Cost: add. (689->108), mult. (838->162), div. (0->0), fcn. (843->6), ass. (0->48)
t43 = cos(pkin(8));
t62 = t43 ^ 2;
t61 = m(5) / 0.2e1;
t58 = -m(5) - m(6);
t57 = -rSges(6,3) - pkin(6);
t41 = pkin(7) + qJ(2);
t39 = cos(t41);
t56 = t39 * t43;
t42 = sin(pkin(8));
t44 = sin(qJ(5));
t45 = cos(qJ(5));
t29 = t42 * t45 - t43 * t44;
t23 = t29 * t39;
t47 = t42 * t44 + t43 * t45;
t24 = t47 * t39;
t55 = t24 * rSges(6,1) + t23 * rSges(6,2);
t38 = sin(t41);
t54 = t39 * pkin(2) + t38 * qJ(3);
t53 = t38 ^ 2 + t39 ^ 2;
t52 = qJ(4) * t42;
t51 = t61 + m(6) / 0.2e1;
t50 = pkin(3) * t56 + t39 * t52 + t54;
t49 = rSges(4,1) * t43 - rSges(4,2) * t42;
t21 = t29 * t38;
t22 = t47 * t38;
t48 = t22 * rSges(6,1) + t21 * rSges(6,2);
t34 = t39 * qJ(3);
t2 = t34 + t57 * t39 + (-t52 - pkin(2) + (-pkin(3) - pkin(4)) * t43) * t38 - t48;
t3 = pkin(4) * t56 + t57 * t38 + t50 + t55;
t46 = m(6) * (t2 * t39 + t3 * t38);
t27 = t39 * rSges(3,1) - t38 * rSges(3,2);
t26 = -t38 * rSges(3,1) - t39 * rSges(3,2);
t18 = t29 * rSges(6,1) - rSges(6,2) * t47;
t17 = Icges(6,1) * t29 - Icges(6,4) * t47;
t16 = Icges(6,4) * t29 - Icges(6,2) * t47;
t15 = Icges(6,5) * t29 - Icges(6,6) * t47;
t13 = t38 * rSges(4,3) + t49 * t39 + t54;
t12 = t39 * rSges(4,3) + t34 + (-pkin(2) - t49) * t38;
t11 = Icges(6,1) * t24 + Icges(6,4) * t23 - Icges(6,5) * t38;
t10 = Icges(6,1) * t22 + Icges(6,4) * t21 + Icges(6,5) * t39;
t9 = Icges(6,4) * t24 + Icges(6,2) * t23 - Icges(6,6) * t38;
t8 = Icges(6,4) * t22 + Icges(6,2) * t21 + Icges(6,6) * t39;
t7 = Icges(6,5) * t24 + Icges(6,6) * t23 - Icges(6,3) * t38;
t6 = Icges(6,5) * t22 + Icges(6,6) * t21 + Icges(6,3) * t39;
t5 = t38 * rSges(5,2) + (rSges(5,1) * t43 + rSges(5,3) * t42) * t39 + t50;
t4 = t39 * rSges(5,2) + t34 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t43 + (-rSges(5,3) - qJ(4)) * t42) * t38;
t1 = -t38 * t48 - t39 * t55;
t14 = [m(2) + m(3) + m(4) - t58; 0; -t47 * t16 + t29 * t17 + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t62 + ((Icges(4,1) + Icges(5,1)) * t42 + 0.2e1 * (Icges(4,4) - Icges(5,5)) * t43) * t42 + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2); 0; m(4) * (t38 * t12 - t39 * t13) + m(5) * (t38 * t4 - t39 * t5) + m(6) * (t38 * t2 - t39 * t3); 0.2e1 * (m(4) / 0.2e1 + t51) * t53; t58 * t43; 0.2e1 * ((t38 * t5 + t39 * t4) * t61 + t46 / 0.2e1) * t42; 0; 0.2e1 * t51 * (t53 * t42 ^ 2 + t62); m(6) * t1; t18 * t46 - (t29 * t11 - t38 * t15 + t23 * t16 + t24 * t17 - t47 * t9) * t38 / 0.2e1 + (t29 * t10 + t39 * t15 + t21 * t16 + t22 * t17 - t47 * t8) * t39 / 0.2e1; 0; m(6) * (t53 * t42 * t18 - t1 * t43); m(6) * (t53 * t18 ^ 2 + t1 ^ 2) - t38 * (-(t24 * t11 + t23 * t9 - t38 * t7) * t38 + (t24 * t10 + t23 * t8 - t38 * t6) * t39) + t39 * (-(t22 * t11 + t21 * t9 + t39 * t7) * t38 + (t22 * t10 + t21 * t8 + t39 * t6) * t39);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
