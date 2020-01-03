% Calculate Gravitation load on the joints for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:41
% EndTime: 2019-12-31 19:40:43
% DurationCPUTime: 0.45s
% Computational Cost: add. (148->101), mult. (294->137), div. (0->0), fcn. (262->6), ass. (0->47)
t49 = rSges(6,3) + pkin(7);
t23 = cos(qJ(1));
t20 = sin(qJ(1));
t53 = g(2) * t20;
t57 = g(1) * t23 + t53;
t56 = -m(5) - m(6);
t55 = -pkin(2) - pkin(3);
t22 = cos(qJ(2));
t52 = g(3) * t22;
t19 = sin(qJ(2));
t51 = t19 * pkin(4);
t15 = t22 * pkin(2);
t48 = t19 * rSges(5,1);
t47 = t19 * t23;
t18 = sin(qJ(5));
t46 = t20 * t18;
t21 = cos(qJ(5));
t45 = t20 * t21;
t44 = t20 * t22;
t43 = t22 * rSges(4,1);
t42 = t22 * rSges(5,2);
t41 = t22 * t23;
t40 = t23 * t18;
t39 = t23 * t21;
t12 = t19 * qJ(3);
t38 = t12 + t15;
t37 = t23 * pkin(1) + t20 * pkin(6);
t36 = qJ(3) * t22;
t35 = -rSges(5,3) - qJ(4);
t34 = rSges(5,2) + t55;
t33 = t22 * pkin(3) + t38;
t32 = -t49 + t55;
t31 = -pkin(1) - t12;
t30 = pkin(2) * t41 + t23 * t12 + t37;
t29 = g(1) * t34;
t28 = pkin(3) * t41 + t30;
t27 = g(1) * t32;
t26 = t22 * rSges(3,1) - t19 * rSges(3,2);
t24 = rSges(6,1) * t21 - rSges(6,2) * t18 + pkin(4);
t16 = t23 * pkin(6);
t9 = t23 * t36;
t7 = t20 * t36;
t5 = t19 * t39 - t46;
t4 = -t19 * t40 - t45;
t3 = -t19 * t45 - t40;
t2 = t19 * t46 - t39;
t1 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) - t20 * rSges(2,2))) - m(3) * (g(1) * (t23 * rSges(3,3) + t16) + g(2) * (rSges(3,1) * t41 - rSges(3,2) * t47 + t37) + (g(1) * (-pkin(1) - t26) + g(2) * rSges(3,3)) * t20) - m(4) * (g(1) * (t23 * rSges(4,2) + t16) + g(2) * (rSges(4,1) * t41 + rSges(4,3) * t47 + t30) + (g(1) * (-t19 * rSges(4,3) - t15 + t31 - t43) + g(2) * rSges(4,2)) * t20) - m(5) * (g(1) * t16 + g(2) * t28 + (g(1) * t35 + g(2) * (-t42 + t48)) * t23 + (g(1) * (t31 - t48) + g(2) * t35 + t22 * t29) * t20) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t16) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t28) + (-g(1) * qJ(4) + g(2) * (t49 * t22 + t51)) * t23 + (g(1) * (t31 - t51) - g(2) * qJ(4) + t22 * t27) * t20), -m(3) * (g(3) * t26 + t57 * (-rSges(3,1) * t19 - rSges(3,2) * t22)) - m(4) * (g(1) * (rSges(4,3) * t41 + t9) + g(2) * (rSges(4,3) * t44 + t7) + g(3) * (t38 + t43) + (g(3) * rSges(4,3) + t57 * (-rSges(4,1) - pkin(2))) * t19) - m(5) * (g(1) * (rSges(5,1) * t41 + t9) + g(2) * (rSges(5,1) * t44 + t7) + g(3) * (t33 - t42) + (g(3) * rSges(5,1) + t23 * t29 + t34 * t53) * t19) - m(6) * (g(1) * t9 + g(2) * t7 + g(3) * t33 + (g(3) * t49 + t57 * t24) * t22 + (g(3) * t24 + t23 * t27 + t32 * t53) * t19), (-m(4) + t56) * (t57 * t19 - t52), t56 * (-g(1) * t20 + g(2) * t23), -m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2)) + g(2) * (-t2 * rSges(6,1) + t3 * rSges(6,2)) + (rSges(6,1) * t18 + rSges(6,2) * t21) * t52)];
taug = t1(:);
