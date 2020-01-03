% Calculate joint inertia matrix for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:21
% EndTime: 2019-12-31 16:17:22
% DurationCPUTime: 0.24s
% Computational Cost: add. (263->49), mult. (578->85), div. (0->0), fcn. (654->6), ass. (0->29)
t28 = sin(qJ(4));
t49 = Icges(5,5) * t28;
t29 = cos(qJ(4));
t48 = Icges(5,6) * t29;
t47 = t49 / 0.2e1 + t48 / 0.2e1;
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t14 = -t26 * t41 - t27 * t42;
t15 = -t26 * t42 + t27 * t41;
t46 = t15 * t14;
t45 = t14 ^ 2;
t44 = t15 ^ 2;
t21 = -t28 * rSges(5,1) - t29 * rSges(5,2);
t43 = m(5) * t21;
t40 = rSges(5,1) * t29;
t39 = rSges(5,2) * t28;
t38 = -t14 * rSges(5,3) - t15 * t40;
t33 = t39 - t40;
t30 = -Icges(5,5) * t29 + Icges(5,6) * t28;
t11 = t14 * rSges(4,1) + t15 * rSges(4,2);
t10 = -t15 * rSges(4,1) + t14 * rSges(4,2);
t5 = Icges(5,3) * t15 + t30 * t14;
t4 = -Icges(5,3) * t14 + t30 * t15;
t3 = (-rSges(5,3) - pkin(5)) * t15 + (pkin(3) - t33) * t14;
t2 = -t14 * pkin(5) + (-pkin(3) + t39) * t15 + t38;
t1 = t15 * (t15 * t39 + t38) + (t15 * rSges(5,3) + t33 * t14) * t14;
t6 = [m(2) + m(3) + m(4) + m(5); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * (t26 ^ 2 + t27 ^ 2); 0; m(4) * (t10 * t26 - t11 * t27) + m(5) * (t2 * t26 - t3 * t27); t29 ^ 2 * Icges(5,2) + Icges(4,3) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,1) * t28 + 0.2e1 * Icges(5,4) * t29) * t28; m(5) * t1; (-t14 * t26 + t15 * t27) * t43; (-t44 / 0.2e1 - t45 / 0.2e1) * (-t48 - t49) + (t47 * t15 - t3 * t43) * t15 + (t47 * t14 - t2 * t43) * t14; m(5) * (t1 ^ 2 + (t44 + t45) * t21 ^ 2) + t15 * (-t4 * t46 + t44 * t5) - t14 * (t45 * t4 - t5 * t46);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
