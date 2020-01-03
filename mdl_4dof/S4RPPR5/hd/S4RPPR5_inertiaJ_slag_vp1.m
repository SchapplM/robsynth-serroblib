% Calculate joint inertia matrix for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:44
% DurationCPUTime: 0.30s
% Computational Cost: add. (313->71), mult. (646->107), div. (0->0), fcn. (706->6), ass. (0->37)
t34 = sin(qJ(4));
t60 = Icges(5,5) * t34;
t36 = cos(qJ(4));
t59 = Icges(5,6) * t36;
t58 = -t60 / 0.2e1 - t59 / 0.2e1;
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t46 = sin(pkin(6));
t47 = cos(pkin(6));
t16 = -t35 * t46 - t37 * t47;
t17 = -t35 * t47 + t37 * t46;
t57 = t17 * t16;
t56 = t16 ^ 2;
t55 = t17 ^ 2;
t23 = -t34 * rSges(5,1) - t36 * rSges(5,2);
t54 = m(5) * t23;
t53 = rSges(5,1) * t36;
t52 = rSges(5,2) * t34;
t51 = t17 * rSges(5,3) - t16 * t53;
t50 = t37 * pkin(1) + t35 * qJ(2);
t45 = t37 * pkin(2) + t50;
t31 = t37 * qJ(2);
t44 = t31 + (-pkin(1) - pkin(2)) * t35;
t41 = t52 - t53;
t38 = -Icges(5,5) * t36 + Icges(5,6) * t34;
t25 = t37 * rSges(2,1) - t35 * rSges(2,2);
t24 = -t35 * rSges(2,1) - t37 * rSges(2,2);
t13 = t37 * rSges(3,1) + t35 * rSges(3,3) + t50;
t12 = t37 * rSges(3,3) + t31 + (-rSges(3,1) - pkin(1)) * t35;
t11 = -t16 * rSges(4,1) - t17 * rSges(4,2) + t45;
t10 = t17 * rSges(4,1) - t16 * rSges(4,2) + t44;
t5 = Icges(5,3) * t17 + t38 * t16;
t4 = -Icges(5,3) * t16 + t38 * t17;
t3 = t17 * pkin(5) + (-pkin(3) + t52) * t16 + t45 + t51;
t2 = (rSges(5,3) + pkin(5)) * t16 + (pkin(3) - t41) * t17 + t44;
t1 = t16 * (t16 * t52 + t51) + (-t16 * rSges(5,3) + t41 * t17) * t17;
t6 = [t36 ^ 2 * Icges(5,2) + Icges(3,2) + Icges(2,3) + Icges(4,3) + m(2) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,1) * t34 + 0.2e1 * Icges(5,4) * t36) * t34; m(3) * (t35 * t12 - t37 * t13) + m(4) * (t35 * t10 - t37 * t11) + m(5) * (t35 * t2 - t37 * t3); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t35 ^ 2 + t37 ^ 2); 0; 0; m(4) + m(5); (t55 / 0.2e1 + t56 / 0.2e1) * (-t59 - t60) + (t58 * t17 - t3 * t54) * t17 + (t58 * t16 - t2 * t54) * t16; (-t16 * t35 + t17 * t37) * t54; -m(5) * t1; m(5) * (t1 ^ 2 + (t55 + t56) * t23 ^ 2) + t17 * (-t4 * t57 + t55 * t5) - t16 * (t56 * t4 - t5 * t57);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
