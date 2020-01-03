% Calculate joint inertia matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.36s
% Computational Cost: add. (577->75), mult. (586->111), div. (0->0), fcn. (466->6), ass. (0->45)
t51 = qJ(1) + qJ(2);
t48 = sin(t51);
t46 = t48 ^ 2;
t49 = cos(t51);
t47 = t49 ^ 2;
t54 = cos(qJ(4));
t80 = Icges(5,5) * t54;
t52 = sin(qJ(4));
t79 = Icges(5,6) * t52;
t27 = -t79 + t80;
t78 = rSges(5,1) * t52 + rSges(5,2) * t54;
t77 = t48 * t49;
t53 = sin(qJ(1));
t74 = t53 * pkin(1);
t71 = t78 * t49;
t70 = t49 * pkin(2) + t48 * qJ(3);
t69 = t46 + t47;
t66 = t49 * rSges(5,3) + t78 * t48;
t65 = t27 * t47 + (t80 / 0.2e1 - t79 / 0.2e1 + t27 / 0.2e1) * t46;
t23 = t49 * rSges(3,1) - t48 * rSges(3,2);
t22 = -t48 * rSges(3,1) - t49 * rSges(3,2);
t38 = t49 * qJ(3);
t12 = t49 * rSges(4,3) + t38 + (rSges(4,2) - pkin(2)) * t48;
t59 = Icges(5,5) * t52 + Icges(5,6) * t54;
t13 = -t49 * rSges(4,2) + t48 * rSges(4,3) + t70;
t9 = t49 * pkin(6) + t66 + t70;
t58 = Icges(5,1) * t54 ^ 2 + Icges(4,1) + Icges(3,3) + (-0.2e1 * Icges(5,4) * t54 + Icges(5,2) * t52) * t52;
t8 = t38 + (-rSges(5,3) - pkin(2) - pkin(6)) * t48 + t71;
t6 = t8 - t74;
t55 = cos(qJ(1));
t50 = t55 * pkin(1);
t7 = t50 + t9;
t57 = m(5) * (t48 * t6 - t49 * t7);
t56 = m(5) * (t48 * t8 - t49 * t9);
t32 = t55 * rSges(2,1) - t53 * rSges(2,2);
t31 = t54 * rSges(5,1) - t52 * rSges(5,2);
t30 = -t53 * rSges(2,1) - t55 * rSges(2,2);
t21 = t23 + t50;
t20 = t22 - t74;
t15 = Icges(5,3) * t48 - t59 * t49;
t14 = Icges(5,3) * t49 + t59 * t48;
t11 = t13 + t50;
t10 = t12 - t74;
t3 = t49 * (t48 * rSges(5,3) - t71) - t48 * t66;
t1 = [Icges(2,3) + m(2) * (t30 ^ 2 + t32 ^ 2) + m(3) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t58; m(3) * (t22 * t20 + t23 * t21) + m(4) * (t12 * t10 + t13 * t11) + m(5) * (t8 * t6 + t9 * t7) + t58; m(3) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t58; m(4) * (t48 * t10 - t49 * t11) + t57; m(4) * (t48 * t12 - t49 * t13) + t56; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t69; t31 * t57 + t65; t31 * t56 + t65; m(5) * t69 * t31; m(5) * (t69 * t31 ^ 2 + t3 ^ 2) + t49 * (t47 * t14 + t15 * t77) + t48 * (t14 * t77 + t46 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
