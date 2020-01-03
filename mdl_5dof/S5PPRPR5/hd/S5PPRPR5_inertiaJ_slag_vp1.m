% Calculate joint inertia matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:24
% EndTime: 2019-12-31 17:33:25
% DurationCPUTime: 0.33s
% Computational Cost: add. (374->62), mult. (804->97), div. (0->0), fcn. (898->6), ass. (0->34)
t35 = sin(pkin(7));
t36 = cos(pkin(7));
t54 = sin(qJ(3));
t55 = cos(qJ(3));
t22 = -t35 * t54 - t36 * t55;
t21 = t22 ^ 2;
t23 = -t35 * t55 + t36 * t54;
t57 = t23 ^ 2;
t49 = t21 + t57;
t37 = sin(qJ(5));
t38 = cos(qJ(5));
t61 = rSges(6,1) * t37 + rSges(6,2) * t38;
t48 = m(5) / 0.2e1 + m(6) / 0.2e1;
t60 = 0.2e1 * t48;
t59 = t23 * t22;
t58 = t61 * t23;
t30 = -t38 * rSges(6,1) + t37 * rSges(6,2);
t56 = m(6) * t30;
t47 = -t23 * rSges(6,3) - t61 * t22;
t20 = t23 * pkin(3);
t2 = -t23 * pkin(6) - t22 * qJ(4) - t20 + t47;
t18 = t23 * qJ(4);
t3 = -t18 - t58 + (rSges(6,3) + pkin(3) + pkin(6)) * t22;
t46 = t23 * t2 - t22 * t3;
t42 = t22 * t36 + t23 * t35;
t39 = Icges(6,5) * t37 + Icges(6,6) * t38;
t14 = t22 * rSges(4,1) + t23 * rSges(4,2);
t13 = -t23 * rSges(4,1) + t22 * rSges(4,2);
t7 = -Icges(6,3) * t23 - t22 * t39;
t6 = -Icges(6,3) * t22 + t23 * t39;
t5 = -t23 * rSges(5,3) - t18 + (-rSges(5,2) + pkin(3)) * t22;
t4 = t23 * rSges(5,2) - t20 + (-rSges(5,3) - qJ(4)) * t22;
t1 = -t22 * t47 + (-t22 * rSges(6,3) + t58) * t23;
t8 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t48) * (t35 ^ 2 + t36 ^ 2); 0; m(4) * (t13 * t35 - t14 * t36) + m(5) * (t4 * t35 - t5 * t36) + m(6) * (t2 * t35 - t3 * t36); t38 ^ 2 * Icges(6,1) + Icges(5,1) + Icges(4,3) + m(4) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Icges(6,4) * t38 + Icges(6,2) * t37) * t37; 0; t42 * t60; m(5) * (-t22 * t5 + t23 * t4) + m(6) * t46; t49 * t60; m(6) * t1; -t42 * t56; -t46 * t56 + t49 * (Icges(6,5) * t38 - Icges(6,6) * t37); -t49 * t56; m(6) * (t49 * t30 ^ 2 + t1 ^ 2) - t22 * (t21 * t6 + t7 * t59) - t23 * (t57 * t7 + t6 * t59);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
