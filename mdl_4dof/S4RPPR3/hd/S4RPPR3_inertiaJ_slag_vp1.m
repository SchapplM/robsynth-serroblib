% Calculate joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:51
% DurationCPUTime: 0.32s
% Computational Cost: add. (422->64), mult. (352->93), div. (0->0), fcn. (288->8), ass. (0->38)
t33 = qJ(1) + pkin(6);
t28 = sin(t33);
t25 = t28 ^ 2;
t30 = cos(t33);
t26 = t30 ^ 2;
t49 = t25 + t26;
t54 = t28 * t30;
t37 = sin(qJ(1));
t53 = t37 * pkin(1);
t32 = pkin(7) + qJ(4);
t29 = cos(t32);
t52 = rSges(5,1) * t29;
t27 = sin(t32);
t51 = rSges(5,2) * t27;
t50 = t28 * rSges(5,3) + t30 * t52;
t46 = rSges(4,3) + qJ(3);
t43 = -t51 + t52;
t40 = Icges(5,5) * t29 - Icges(5,6) * t27;
t34 = sin(pkin(7));
t35 = cos(pkin(7));
t39 = rSges(4,1) * t35 - rSges(4,2) * t34 + pkin(2);
t38 = cos(qJ(1));
t36 = -pkin(5) - qJ(3);
t31 = t38 * pkin(1);
t24 = t35 * pkin(3) + pkin(2);
t22 = t38 * rSges(2,1) - t37 * rSges(2,2);
t21 = -t37 * rSges(2,1) - t38 * rSges(2,2);
t18 = t27 * rSges(5,1) + t29 * rSges(5,2);
t13 = t30 * rSges(3,1) - t28 * rSges(3,2) + t31;
t12 = -t28 * rSges(3,1) - t30 * rSges(3,2) - t53;
t7 = Icges(5,3) * t28 + t40 * t30;
t6 = -Icges(5,3) * t30 + t40 * t28;
t5 = t46 * t28 + t39 * t30 + t31;
t4 = -t39 * t28 + t46 * t30 - t53;
t3 = -t28 * t36 + t31 + (t24 - t51) * t30 + t50;
t2 = -t53 + (rSges(5,3) - t36) * t30 + (-t24 - t43) * t28;
t1 = t30 * (-t30 * t51 + t50) + (-t30 * rSges(5,3) + t43 * t28) * t28;
t8 = [Icges(4,2) * t35 ^ 2 + Icges(5,2) * t29 ^ 2 + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t34 + 0.2e1 * Icges(4,4) * t35) * t34 + m(2) * (t21 ^ 2 + t22 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,1) * t27 + 0.2e1 * Icges(5,4) * t29) * t27; 0; m(3) + m(4) + m(5); m(4) * (t28 * t4 - t30 * t5) + m(5) * (t28 * t2 - t30 * t3); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t49; m(5) * (-t2 * t30 - t28 * t3) * t18 + t49 * (Icges(5,5) * t27 + Icges(5,6) * t29); m(5) * t1; 0; m(5) * (t49 * t18 ^ 2 + t1 ^ 2) + t28 * (t25 * t7 - t6 * t54) - t30 * (t26 * t6 - t7 * t54);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
