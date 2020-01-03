% Calculate Gravitation load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:02
% EndTime: 2019-12-31 19:43:04
% DurationCPUTime: 0.60s
% Computational Cost: add. (208->121), mult. (469->164), div. (0->0), fcn. (483->8), ass. (0->47)
t31 = cos(qJ(1));
t28 = sin(qJ(1));
t59 = g(2) * t28;
t63 = g(1) * t31 + t59;
t62 = -m(5) - m(6);
t61 = g(1) * t28;
t27 = sin(qJ(2));
t58 = g(3) * t27;
t30 = cos(qJ(2));
t21 = t30 * pkin(2);
t57 = -pkin(7) - rSges(6,3);
t24 = sin(pkin(8));
t56 = t24 * t30;
t25 = cos(pkin(8));
t55 = t25 * t30;
t54 = t27 * t31;
t53 = t28 * t30;
t52 = t30 * t31;
t51 = t31 * t24;
t19 = t27 * qJ(3);
t50 = t19 + t21;
t49 = t31 * pkin(1) + t28 * pkin(6);
t48 = qJ(3) * t30;
t47 = rSges(5,3) + qJ(4);
t46 = -pkin(1) - t21;
t45 = t57 * t31;
t44 = pkin(3) * t55 + qJ(4) * t56 + t50;
t43 = pkin(2) * t52 + t31 * t19 + t49;
t22 = t31 * pkin(6);
t7 = t24 * t53 + t25 * t31;
t8 = t25 * t53 - t51;
t42 = -t8 * pkin(3) - t7 * qJ(4) + t22;
t10 = t28 * t24 + t25 * t52;
t41 = t10 * pkin(3) + t43;
t26 = sin(qJ(5));
t29 = cos(qJ(5));
t40 = t26 * t8 - t29 * t7;
t39 = -t26 * t7 - t29 * t8;
t38 = rSges(3,1) * t30 - rSges(3,2) * t27;
t36 = t24 * t29 - t25 * t26;
t35 = t24 * t26 + t25 * t29;
t15 = t31 * t48;
t12 = t28 * t48;
t9 = -t28 * t25 + t30 * t51;
t3 = t10 * t29 + t26 * t9;
t2 = -t10 * t26 + t29 * t9;
t1 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t28 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t31 + t22) + g(2) * (rSges(3,1) * t52 - rSges(3,2) * t54 + t49) + (g(1) * (-pkin(1) - t38) + g(2) * rSges(3,3)) * t28) - m(4) * (g(1) * (-rSges(4,1) * t8 + rSges(4,2) * t7 + t22) + g(2) * (t10 * rSges(4,1) - t9 * rSges(4,2) + rSges(4,3) * t54 + t43) + ((-rSges(4,3) - qJ(3)) * t27 + t46) * t61) - m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,3) * t7 + t42) + g(2) * (t10 * rSges(5,1) + rSges(5,2) * t54 + t47 * t9 + t41) + ((-rSges(5,2) - qJ(3)) * t27 + t46) * t61) - m(6) * (g(1) * (t39 * rSges(6,1) + t40 * rSges(6,2) - t8 * pkin(4) + t42) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t10 * pkin(4) + t9 * qJ(4) + t27 * t45 + t41) + ((-qJ(3) - t57) * t27 + t46) * t61), -m(3) * (g(3) * t38 + t63 * (-rSges(3,1) * t27 - rSges(3,2) * t30)) - m(4) * (g(1) * (rSges(4,3) * t52 + t15) + g(2) * (rSges(4,3) * t53 + t12) + g(3) * (rSges(4,1) * t55 - rSges(4,2) * t56 + t50) + (g(3) * rSges(4,3) + t63 * (-rSges(4,1) * t25 + rSges(4,2) * t24 - pkin(2))) * t27) - m(5) * (g(1) * (rSges(5,2) * t52 + t15) + g(2) * (rSges(5,2) * t53 + t12) + g(3) * (rSges(5,1) * t55 + rSges(5,3) * t56 + t44) + (g(3) * rSges(5,2) + t63 * (-pkin(2) + (-rSges(5,1) - pkin(3)) * t25 - t47 * t24)) * t27) - m(6) * (g(1) * t15 + g(2) * t12 + g(3) * t44 + (g(3) * (t35 * rSges(6,1) + t36 * rSges(6,2) + t25 * pkin(4)) + g(1) * t45 + t57 * t59) * t30 + (g(3) * t57 + t63 * (-pkin(2) + (-t26 * rSges(6,1) - t29 * rSges(6,2) - qJ(4)) * t24 + (-t29 * rSges(6,1) + t26 * rSges(6,2) - pkin(3) - pkin(4)) * t25)) * t27), (-m(4) + t62) * (-g(3) * t30 + t63 * t27), t62 * (g(1) * t9 + g(2) * t7 + t24 * t58), -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (-t40 * rSges(6,1) + t39 * rSges(6,2)) + (t36 * rSges(6,1) - t35 * rSges(6,2)) * t58)];
taug = t1(:);
