% Calculate joint inertia matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:27
% EndTime: 2019-12-31 16:18:28
% DurationCPUTime: 0.41s
% Computational Cost: add. (1082->100), mult. (1478->179), div. (0->0), fcn. (1601->6), ass. (0->59)
t50 = cos(pkin(6));
t47 = t50 ^ 2;
t49 = sin(pkin(6));
t63 = t49 ^ 2 + t47;
t72 = t49 / 0.2e1;
t48 = pkin(7) + qJ(3);
t44 = sin(t48);
t71 = t44 * t49;
t70 = t44 * t50;
t45 = cos(t48);
t51 = sin(qJ(4));
t52 = cos(qJ(4));
t26 = -Icges(5,3) * t45 + (Icges(5,5) * t52 - Icges(5,6) * t51) * t44;
t69 = t45 * t26;
t68 = t49 * t51;
t67 = t49 * t52;
t66 = t50 * t51;
t65 = t50 * t52;
t29 = -t45 * rSges(5,3) + (rSges(5,1) * t52 - rSges(5,2) * t51) * t44;
t64 = -t44 * pkin(3) + t45 * pkin(5) - t29;
t62 = Icges(5,5) * t44;
t61 = Icges(5,6) * t44;
t60 = Icges(5,3) * t44;
t59 = pkin(3) * t45 + pkin(5) * t44;
t53 = Icges(4,5) * t45 - Icges(4,6) * t44;
t41 = t44 * rSges(4,1) + t45 * rSges(4,2);
t39 = t45 * t65 + t68;
t38 = -t45 * t66 + t67;
t37 = t45 * t67 - t66;
t36 = -t45 * t68 - t65;
t30 = -Icges(4,3) * t50 + t53 * t49;
t28 = -Icges(5,5) * t45 + (Icges(5,1) * t52 - Icges(5,4) * t51) * t44;
t27 = -Icges(5,6) * t45 + (Icges(5,4) * t52 - Icges(5,2) * t51) * t44;
t25 = t64 * t50;
t24 = t64 * t49;
t23 = t39 * rSges(5,1) + t38 * rSges(5,2) + rSges(5,3) * t70;
t22 = t37 * rSges(5,1) + t36 * rSges(5,2) + rSges(5,3) * t71;
t21 = Icges(5,1) * t39 + Icges(5,4) * t38 + t50 * t62;
t20 = Icges(5,1) * t37 + Icges(5,4) * t36 + t49 * t62;
t19 = Icges(5,4) * t39 + Icges(5,2) * t38 + t50 * t61;
t18 = Icges(5,4) * t37 + Icges(5,2) * t36 + t49 * t61;
t17 = Icges(5,5) * t39 + Icges(5,6) * t38 + t50 * t60;
t16 = Icges(5,5) * t37 + Icges(5,6) * t36 + t49 * t60;
t15 = t63 * (rSges(4,1) * t45 - rSges(4,2) * t44);
t14 = -t45 * t23 - t29 * t70;
t13 = t45 * t22 + t29 * t71;
t12 = (t22 * t50 - t23 * t49) * t44;
t11 = (t59 * t50 + t23) * t50 + (t59 * t49 + t22) * t49;
t10 = -t45 * t17 + (-t19 * t51 + t21 * t52) * t44;
t9 = -t45 * t16 + (-t18 * t51 + t20 * t52) * t44;
t8 = t17 * t70 + t38 * t19 + t39 * t21;
t7 = t16 * t70 + t38 * t18 + t39 * t20;
t6 = t17 * t71 + t36 * t19 + t37 * t21;
t5 = t16 * t71 + t36 * t18 + t37 * t20;
t4 = t8 * t49 - t7 * t50;
t3 = t6 * t49 - t5 * t50;
t2 = -(t38 * t27 + t39 * t28) * t45 + (t7 * t49 + (t8 - t69) * t50) * t44;
t1 = -(t36 * t27 + t37 * t28) * t45 + (t6 * t50 + (t5 - t69) * t49) * t44;
t31 = [m(2) + m(3) + m(4) + m(5); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t63; m(4) * t15 + m(5) * t11; m(5) * (-t24 * t50 + t25 * t49); m(4) * (t63 * t41 ^ 2 + t15 ^ 2) + m(5) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + (-t47 * t30 - t3) * t50 + (-t49 * t30 * t50 + t4 + t63 * (Icges(4,3) * t49 + t53 * t50)) * t49; m(5) * t12; m(5) * (t13 * t49 - t14 * t50); m(5) * (t12 * t11 + t13 * t25 + t14 * t24) + t2 * t72 - t50 * t1 / 0.2e1 - t45 * (t10 * t49 - t9 * t50) / 0.2e1 + (t50 * t4 / 0.2e1 + t3 * t72) * t44; m(5) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) + t2 * t70 + t1 * t71 - t45 * (t45 ^ 2 * t26 + (t10 * t50 + t9 * t49 - (-t27 * t51 + t28 * t52) * t45) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t31(1), t31(2), t31(4), t31(7); t31(2), t31(3), t31(5), t31(8); t31(4), t31(5), t31(6), t31(9); t31(7), t31(8), t31(9), t31(10);];
Mq = res;
