% Calculate joint inertia matrix for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:42
% DurationCPUTime: 0.43s
% Computational Cost: add. (589->100), mult. (1490->183), div. (0->0), fcn. (1607->6), ass. (0->56)
t47 = cos(pkin(6));
t45 = t47 ^ 2;
t46 = sin(pkin(6));
t62 = t46 ^ 2 + t45;
t69 = t47 / 0.2e1;
t51 = cos(qJ(3));
t68 = t46 * t51;
t67 = t47 * t51;
t48 = sin(qJ(4));
t49 = sin(qJ(3));
t66 = t48 * t49;
t50 = cos(qJ(4));
t32 = Icges(5,3) * t49 + (Icges(5,5) * t50 - Icges(5,6) * t48) * t51;
t65 = t49 * t32;
t64 = t49 * t50;
t35 = t49 * rSges(5,3) + (rSges(5,1) * t50 - rSges(5,2) * t48) * t51;
t63 = t51 * pkin(3) + t49 * pkin(5) + t35;
t61 = Icges(5,5) * t51;
t60 = Icges(5,6) * t51;
t59 = Icges(5,3) * t51;
t58 = -pkin(3) * t49 + pkin(5) * t51;
t52 = Icges(4,5) * t49 + Icges(4,6) * t51;
t41 = t51 * rSges(4,1) - t49 * rSges(4,2);
t39 = t46 * t48 - t47 * t64;
t38 = t46 * t50 + t47 * t66;
t37 = t46 * t64 + t47 * t48;
t36 = -t46 * t66 + t47 * t50;
t34 = Icges(5,5) * t49 + (Icges(5,1) * t50 - Icges(5,4) * t48) * t51;
t33 = Icges(5,6) * t49 + (Icges(5,4) * t50 - Icges(5,2) * t48) * t51;
t26 = Icges(4,3) * t47 + t52 * t46;
t25 = t63 * t47;
t24 = t63 * t46;
t23 = t39 * rSges(5,1) + t38 * rSges(5,2) + rSges(5,3) * t67;
t22 = t37 * rSges(5,1) + t36 * rSges(5,2) - rSges(5,3) * t68;
t21 = Icges(5,1) * t39 + Icges(5,4) * t38 + t47 * t61;
t20 = Icges(5,1) * t37 + Icges(5,4) * t36 - t46 * t61;
t19 = Icges(5,4) * t39 + Icges(5,2) * t38 + t47 * t60;
t18 = Icges(5,4) * t37 + Icges(5,2) * t36 - t46 * t60;
t17 = Icges(5,5) * t39 + Icges(5,6) * t38 + t47 * t59;
t16 = Icges(5,5) * t37 + Icges(5,6) * t36 - t46 * t59;
t15 = t62 * (-rSges(4,1) * t49 - rSges(4,2) * t51);
t14 = -t49 * t23 + t35 * t67;
t13 = t49 * t22 + t35 * t68;
t12 = (-t22 * t47 - t23 * t46) * t51;
t11 = (t58 * t47 + t23) * t47 + (t58 * t46 - t22) * t46;
t10 = t49 * t17 + (-t19 * t48 + t21 * t50) * t51;
t9 = t49 * t16 + (-t18 * t48 + t20 * t50) * t51;
t8 = t17 * t67 + t38 * t19 + t39 * t21;
t7 = t16 * t67 + t38 * t18 + t39 * t20;
t6 = -t17 * t68 + t36 * t19 + t37 * t21;
t5 = -t16 * t68 + t36 * t18 + t37 * t20;
t4 = t8 * t46 + t7 * t47;
t3 = t6 * t46 + t5 * t47;
t2 = (t38 * t33 + t39 * t34) * t49 + (-t7 * t46 + (t8 + t65) * t47) * t51;
t1 = (t36 * t33 + t37 * t34) * t49 + (t6 * t47 + (-t5 - t65) * t46) * t51;
t27 = [m(2) + m(3) + m(4) + m(5); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t62; m(4) * t15 + m(5) * t11; m(5) * (t24 * t46 + t25 * t47) + m(4) * t62 * t41; m(4) * (t62 * t41 ^ 2 + t15 ^ 2) + m(5) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + (t45 * t26 + t3) * t47 + (t46 * t26 * t47 + t4 + t62 * (Icges(4,3) * t46 - t52 * t47)) * t46; m(5) * t12; m(5) * (-t13 * t47 + t14 * t46); m(5) * (t12 * t11 - t13 * t25 + t14 * t24) + t1 * t69 + t46 * t2 / 0.2e1 + t49 * (t10 * t46 + t9 * t47) / 0.2e1 + (-t46 * t3 / 0.2e1 + t4 * t69) * t51; m(5) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) - t1 * t68 + t2 * t67 + t49 * (t49 ^ 2 * t32 + (-t9 * t46 + t10 * t47 + (-t33 * t48 + t34 * t50) * t49) * t51);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t27(1), t27(2), t27(4), t27(7); t27(2), t27(3), t27(5), t27(8); t27(4), t27(5), t27(6), t27(9); t27(7), t27(8), t27(9), t27(10);];
Mq = res;
