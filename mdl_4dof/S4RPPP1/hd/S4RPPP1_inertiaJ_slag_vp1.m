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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:35
% EndTime: 2018-11-14 13:45:35
% DurationCPUTime: 0.23s
% Computational Cost: add. (658->89), mult. (822->120), div. (0->0), fcn. (698->10), ass. (0->45)
t55 = m(5) / 0.2e1;
t56 = m(4) / 0.2e1;
t48 = t56 + t55;
t59 = 0.2e1 * t48;
t34 = cos(pkin(4));
t58 = t34 ^ 2;
t57 = m(3) / 0.2e1;
t32 = sin(pkin(4));
t35 = sin(qJ(1));
t54 = t32 * t35;
t36 = cos(qJ(1));
t53 = t32 * t36;
t51 = qJ(2) * t32;
t52 = t36 * pkin(1) + t35 * t51;
t50 = rSges(5,3) + qJ(4);
t49 = 0.2e1 * t34;
t47 = pkin(4) - pkin(6);
t46 = pkin(4) + pkin(6);
t33 = cos(pkin(6));
t26 = -sin(t47) / 0.2e1;
t39 = sin(t46) / 0.2e1;
t37 = t39 + t26;
t14 = t33 * t36 - t35 * t37;
t45 = t14 * pkin(2) + t52;
t44 = t32 * (rSges(5,1) + pkin(3));
t43 = -t35 * pkin(1) + t36 * t51;
t31 = sin(pkin(6));
t40 = cos(t47) / 0.2e1;
t41 = cos(t46);
t38 = t40 + t41 / 0.2e1;
t11 = t31 * t35 - t36 * t38;
t42 = -t11 * qJ(3) + t43;
t22 = rSges(2,1) * t36 - t35 * rSges(2,2);
t21 = -t35 * rSges(2,1) - rSges(2,2) * t36;
t20 = t40 - t41 / 0.2e1;
t18 = t39 - t26;
t13 = t36 * t31 + t35 * t38;
t12 = t35 * t33 + t36 * t37;
t7 = rSges(3,1) * t14 - rSges(3,2) * t13 + rSges(3,3) * t54 + t52;
t6 = -t12 * rSges(3,1) + t11 * rSges(3,2) + rSges(3,3) * t53 + t43;
t4 = rSges(4,1) * t54 - t14 * rSges(4,2) + (rSges(4,3) + qJ(3)) * t13 + t45;
t3 = rSges(4,1) * t53 - t11 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t12 + t42;
t2 = t35 * t44 + t50 * t14 + (rSges(5,2) + qJ(3)) * t13 + t45;
t1 = -t11 * rSges(5,2) + t36 * t44 + (-pkin(2) - t50) * t12 + t42;
t5 = [Icges(2,3) + (Icges(3,3) + Icges(4,1) + Icges(5,1)) * t58 + m(2) * (t21 ^ 2 + t22 ^ 2) + m(3) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + ((Icges(3,1) + Icges(4,2) + Icges(5,3)) * t20 + (-Icges(4,4) + Icges(3,5) + Icges(5,5)) * t49) * t20 + ((Icges(3,2) + Icges(4,3) + Icges(5,2)) * t18 + (-Icges(5,4) - Icges(4,5) + Icges(3,6)) * t49 + 0.2e1 * (Icges(3,4) + Icges(4,6) - Icges(5,6)) * t20) * t18; 0.2e1 * ((t35 * t6 - t36 * t7) * t57 + (t3 * t35 - t36 * t4) * t56 + (t1 * t35 - t2 * t36) * t55) * t32; 0.2e1 * (t57 + t48) * (t58 + (t35 ^ 2 + t36 ^ 2) * t32 ^ 2); m(4) * (t11 * t4 + t13 * t3) + m(5) * (t1 * t13 + t11 * t2); (-t18 * t34 + (-t11 * t36 + t13 * t35) * t32) * t59; (t11 ^ 2 + t13 ^ 2 + t18 ^ 2) * t59; m(5) * (t1 * t14 + t12 * t2); m(5) * (t20 * t34 + (-t12 * t36 + t14 * t35) * t32); m(5) * (t11 * t12 + t13 * t14 - t18 * t20); m(5) * (t12 ^ 2 + t14 ^ 2 + t20 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
