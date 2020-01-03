% Calculate joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.29s
% Computational Cost: add. (339->103), mult. (808->160), div. (0->0), fcn. (813->6), ass. (0->46)
t40 = cos(pkin(6));
t60 = t40 ^ 2;
t59 = m(4) / 0.2e1;
t56 = -rSges(5,3) - pkin(5);
t44 = cos(qJ(1));
t55 = t40 * t44;
t39 = sin(pkin(6));
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t27 = t39 * t43 - t40 * t41;
t23 = t27 * t44;
t46 = t39 * t41 + t40 * t43;
t24 = t46 * t44;
t54 = t24 * rSges(5,1) + t23 * rSges(5,2);
t42 = sin(qJ(1));
t53 = t44 * pkin(1) + t42 * qJ(2);
t52 = t42 ^ 2 + t44 ^ 2;
t51 = qJ(3) * t39;
t50 = t59 + m(5) / 0.2e1;
t49 = pkin(2) * t55 + t44 * t51 + t53;
t48 = rSges(3,1) * t40 - rSges(3,2) * t39;
t21 = t27 * t42;
t22 = t46 * t42;
t47 = t22 * rSges(5,1) + t21 * rSges(5,2);
t34 = t44 * qJ(2);
t2 = t34 + t56 * t44 + (-t51 - pkin(1) + (-pkin(2) - pkin(3)) * t40) * t42 - t47;
t3 = pkin(3) * t55 + t56 * t42 + t49 + t54;
t45 = m(5) * (t2 * t44 + t3 * t42);
t29 = t44 * rSges(2,1) - t42 * rSges(2,2);
t28 = -t42 * rSges(2,1) - t44 * rSges(2,2);
t18 = t42 * rSges(3,3) + t48 * t44 + t53;
t17 = t44 * rSges(3,3) + t34 + (-pkin(1) - t48) * t42;
t16 = t27 * rSges(5,1) - rSges(5,2) * t46;
t15 = Icges(5,1) * t27 - Icges(5,4) * t46;
t14 = Icges(5,4) * t27 - Icges(5,2) * t46;
t13 = Icges(5,5) * t27 - Icges(5,6) * t46;
t11 = t42 * rSges(4,2) + (rSges(4,1) * t40 + rSges(4,3) * t39) * t44 + t49;
t10 = t44 * rSges(4,2) + t34 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t40 + (-rSges(4,3) - qJ(3)) * t39) * t42;
t9 = Icges(5,1) * t24 + Icges(5,4) * t23 - Icges(5,5) * t42;
t8 = Icges(5,1) * t22 + Icges(5,4) * t21 + Icges(5,5) * t44;
t7 = Icges(5,4) * t24 + Icges(5,2) * t23 - Icges(5,6) * t42;
t6 = Icges(5,4) * t22 + Icges(5,2) * t21 + Icges(5,6) * t44;
t5 = Icges(5,5) * t24 + Icges(5,6) * t23 - Icges(5,3) * t42;
t4 = Icges(5,5) * t22 + Icges(5,6) * t21 + Icges(5,3) * t44;
t1 = -t42 * t47 - t44 * t54;
t12 = [-t46 * t14 + t27 * t15 + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t60 + ((Icges(3,1) + Icges(4,1)) * t39 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t40) * t39 + m(2) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2); m(3) * (t42 * t17 - t44 * t18) + m(4) * (t42 * t10 - t44 * t11) + m(5) * (t42 * t2 - t44 * t3); 0.2e1 * (m(3) / 0.2e1 + t50) * t52; 0.2e1 * ((t10 * t44 + t11 * t42) * t59 + t45 / 0.2e1) * t39; 0; 0.2e1 * t50 * (t52 * t39 ^ 2 + t60); t16 * t45 - (-t42 * t13 + t23 * t14 + t24 * t15 + t27 * t9 - t46 * t7) * t42 / 0.2e1 + (t44 * t13 + t21 * t14 + t22 * t15 + t27 * t8 - t46 * t6) * t44 / 0.2e1; 0; m(5) * (t52 * t39 * t16 - t1 * t40); m(5) * (t52 * t16 ^ 2 + t1 ^ 2) - t42 * (-(t23 * t7 + t24 * t9 - t42 * t5) * t42 + (t23 * t6 + t24 * t8 - t42 * t4) * t44) + t44 * (-(t21 * t7 + t22 * t9 + t44 * t5) * t42 + (t21 * t6 + t22 * t8 + t44 * t4) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;
