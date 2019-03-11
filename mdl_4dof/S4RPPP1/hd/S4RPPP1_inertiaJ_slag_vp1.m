% Calculate joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:24
% EndTime: 2019-03-08 18:26:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (280->83), mult. (699->121), div. (0->0), fcn. (697->6), ass. (0->41)
t27 = cos(pkin(4));
t50 = t27 ^ 2;
t25 = sin(pkin(4));
t49 = 0.2e1 * t25;
t48 = m(3) / 0.2e1;
t47 = m(4) / 0.2e1;
t46 = m(5) / 0.2e1;
t26 = cos(pkin(6));
t45 = t25 * t26;
t28 = sin(qJ(1));
t44 = t25 * t28;
t29 = cos(qJ(1));
t43 = t25 * t29;
t24 = sin(pkin(6));
t42 = t28 * t24;
t41 = t28 * t26;
t40 = t29 * t24;
t39 = t29 * t26;
t37 = qJ(2) * t25;
t38 = t29 * pkin(1) + t28 * t37;
t36 = rSges(5,3) + qJ(4);
t35 = 0.2e1 * t27;
t34 = t47 + t46;
t15 = -t27 * t42 + t39;
t33 = t15 * pkin(2) + t38;
t32 = t25 * (rSges(5,1) + pkin(3));
t31 = -t28 * pkin(1) + t29 * t37;
t12 = -t27 * t39 + t42;
t30 = -t12 * qJ(3) + t31;
t23 = t25 ^ 2;
t17 = t29 * rSges(2,1) - t28 * rSges(2,2);
t16 = -t28 * rSges(2,1) - t29 * rSges(2,2);
t14 = t27 * t41 + t40;
t13 = t27 * t40 + t41;
t6 = t15 * rSges(3,1) - t14 * rSges(3,2) + rSges(3,3) * t44 + t38;
t5 = -t13 * rSges(3,1) + t12 * rSges(3,2) + rSges(3,3) * t43 + t31;
t4 = rSges(4,1) * t44 - t15 * rSges(4,2) + (rSges(4,3) + qJ(3)) * t14 + t33;
t3 = rSges(4,1) * t43 - t12 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t13 + t30;
t2 = t28 * t32 + t36 * t15 + (rSges(5,2) + qJ(3)) * t14 + t33;
t1 = -t12 * rSges(5,2) + t29 * t32 + (-pkin(2) - t36) * t13 + t30;
t7 = [Icges(2,3) + (Icges(3,3) + Icges(4,1) + Icges(5,1)) * t50 + m(2) * (t16 ^ 2 + t17 ^ 2) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + (((Icges(3,2) + Icges(4,3) + Icges(5,2)) * t45 + (-Icges(5,4) - Icges(4,5) + Icges(3,6)) * t35) * t26 + ((Icges(3,1) + Icges(4,2) + Icges(5,3)) * t25 * t24 + (Icges(3,5) - Icges(4,4) + Icges(5,5)) * t35 + 0.2e1 * (Icges(3,4) + Icges(4,6) - Icges(5,6)) * t45) * t24) * t25; ((t28 * t5 - t29 * t6) * t48 + (t28 * t3 - t29 * t4) * t47 + (t1 * t28 - t2 * t29) * t46) * t49; 0.2e1 * (t48 + t34) * (t50 + (t28 ^ 2 + t29 ^ 2) * t23); m(4) * (t12 * t4 + t14 * t3) + m(5) * (t14 * t1 + t12 * t2); t34 * (-t12 * t29 + t14 * t28 - t27 * t26) * t49; 0.2e1 * t34 * (t23 * t26 ^ 2 + t12 ^ 2 + t14 ^ 2); m(5) * (t15 * t1 + t13 * t2); m(5) * (-t13 * t29 + t15 * t28 + t27 * t24) * t25; m(5) * (-t23 * t24 * t26 + t13 * t12 + t15 * t14); m(5) * (t23 * t24 ^ 2 + t13 ^ 2 + t15 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
