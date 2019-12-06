% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:53
% EndTime: 2019-12-05 17:55:54
% DurationCPUTime: 0.41s
% Computational Cost: add. (234->102), mult. (294->140), div. (0->0), fcn. (285->10), ass. (0->54)
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t64 = -pkin(2) * t34 - pkin(1) + (-rSges(4,3) - pkin(6)) * t33;
t63 = -m(5) - m(6);
t32 = qJ(3) + pkin(9);
t28 = qJ(5) + t32;
t24 = cos(t28);
t39 = cos(qJ(1));
t49 = t39 * t24;
t23 = sin(t28);
t37 = sin(qJ(1));
t56 = t37 * t23;
t5 = t34 * t56 + t49;
t50 = t39 * t23;
t55 = t37 * t24;
t6 = t34 * t55 - t50;
t62 = t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = t34 * t50 - t55;
t8 = -t34 * t49 - t56;
t61 = -t7 * rSges(6,1) + t8 * rSges(6,2);
t60 = g(1) * t33;
t59 = g(2) * t39;
t36 = sin(qJ(3));
t58 = t36 * pkin(3);
t26 = sin(t32);
t19 = pkin(4) * t26 + t58;
t57 = t19 * t34;
t54 = t37 * t26;
t27 = cos(t32);
t53 = t37 * t27;
t52 = t37 * t36;
t38 = cos(qJ(3));
t51 = t37 * t38;
t48 = t39 * t26;
t47 = t39 * t27;
t46 = t39 * t36;
t45 = t39 * t38;
t35 = -qJ(4) - pkin(6);
t30 = t38 * pkin(3);
t20 = pkin(4) * t27 + t30;
t43 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t42 = -rSges(3,1) * t34 + rSges(3,2) * t33 - pkin(1);
t16 = t34 * t46 - t51;
t14 = t34 * t52 + t45;
t41 = -(t30 + pkin(2)) * t34 - pkin(1) + (-rSges(5,3) + t35) * t33;
t40 = -(pkin(2) + t20) * t34 - pkin(1) + (-rSges(6,3) - pkin(7) + t35) * t33;
t29 = t39 * qJ(2);
t17 = -t34 * t45 - t52;
t15 = t34 * t51 - t46;
t12 = -t34 * t47 - t54;
t11 = t34 * t48 - t53;
t10 = t34 * t53 - t48;
t9 = t34 * t54 + t47;
t1 = [-m(2) * (g(2) * (-t39 * rSges(2,1) + t37 * rSges(2,2)) + g(3) * (-t37 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(3) * t29 + (g(3) * rSges(3,3) + g(2) * t42) * t39 + (g(2) * (-rSges(3,3) - qJ(2)) + g(3) * t42) * t37) - m(4) * (g(2) * (t17 * rSges(4,1) + t16 * rSges(4,2)) + g(3) * (-t15 * rSges(4,1) + t14 * rSges(4,2) + t29) + t64 * t59 + (-g(2) * qJ(2) + g(3) * t64) * t37) - m(5) * (g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2)) + g(3) * (-t10 * rSges(5,1) + t9 * rSges(5,2) + t29) + (g(2) * t41 + g(3) * t58) * t39 + (g(2) * (-qJ(2) - t58) + g(3) * t41) * t37) - m(6) * (g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2)) + g(3) * (-t6 * rSges(6,1) + t5 * rSges(6,2) + t29) + (g(2) * t40 + g(3) * t19) * t39 + (g(2) * (-qJ(2) - t19) + g(3) * t40) * t37), (-m(3) - m(4) + t63) * (g(3) * t37 + t59), -m(4) * (g(2) * (t14 * rSges(4,1) + t15 * rSges(4,2)) + g(3) * (-t16 * rSges(4,1) + t17 * rSges(4,2))) - m(5) * (g(2) * (t9 * rSges(5,1) + t10 * rSges(5,2) + t14 * pkin(3)) + g(3) * (-t11 * rSges(5,1) + t12 * rSges(5,2) - t16 * pkin(3))) - m(6) * (g(2) * (t39 * t20 + t37 * t57 + t62) + g(3) * (t37 * t20 - t39 * t57 + t61)) + (-m(4) * (-rSges(4,1) * t36 - rSges(4,2) * t38) - m(5) * (-rSges(5,1) * t26 - rSges(5,2) * t27 - t58) - m(6) * (-t19 + t43)) * t60, t63 * (-g(1) * t34 + (-g(2) * t37 + g(3) * t39) * t33), -m(6) * (g(2) * t62 + g(3) * t61 + t43 * t60)];
taug = t1(:);
