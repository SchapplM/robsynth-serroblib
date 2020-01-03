% Calculate Gravitation load on the joints for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:57
% EndTime: 2019-12-31 18:31:58
% DurationCPUTime: 0.35s
% Computational Cost: add. (200->77), mult. (239->99), div. (0->0), fcn. (210->8), ass. (0->40)
t17 = pkin(8) + qJ(3);
t15 = sin(t17);
t12 = t15 * qJ(4);
t16 = cos(t17);
t13 = t16 * pkin(3);
t37 = t12 + t13;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t50 = g(1) * t24 + g(2) * t22;
t49 = -m(5) - m(6);
t46 = g(3) * t16;
t45 = rSges(6,3) + pkin(7);
t20 = -pkin(6) - qJ(2);
t44 = pkin(4) - t20;
t21 = sin(qJ(5));
t43 = t22 * t21;
t23 = cos(qJ(5));
t42 = t22 * t23;
t41 = t24 * t21;
t40 = t24 * t23;
t39 = rSges(5,1) - t20;
t38 = rSges(4,3) - t20;
t35 = rSges(3,3) + qJ(2);
t34 = -pkin(3) - t45;
t19 = cos(pkin(8));
t14 = pkin(2) * t19 + pkin(1);
t9 = t24 * t14;
t33 = t24 * t37 + t9;
t32 = t45 * t16;
t31 = -t14 - t12;
t30 = t50 * qJ(4) * t16;
t29 = rSges(4,1) * t16 - rSges(4,2) * t15;
t28 = rSges(6,1) * t21 + rSges(6,2) * t23;
t27 = -t16 * rSges(5,2) + t15 * rSges(5,3);
t26 = rSges(3,1) * t19 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t5 = -t15 * t43 + t40;
t4 = t15 * t42 + t41;
t3 = t15 * t41 + t42;
t2 = t15 * t40 - t43;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t22 - rSges(2,2) * t24) + g(2) * (rSges(2,1) * t24 - rSges(2,2) * t22)) - m(3) * ((g(1) * t35 + g(2) * t26) * t24 + (-g(1) * t26 + g(2) * t35) * t22) - m(4) * (g(2) * t9 + (g(1) * t38 + g(2) * t29) * t24 + (g(1) * (-t14 - t29) + g(2) * t38) * t22) - m(5) * (g(2) * t33 + (g(1) * t39 + g(2) * t27) * t24 + (g(1) * (-t27 + t31 - t13) + g(2) * t39) * t22) - m(6) * ((t3 * rSges(6,1) + t2 * rSges(6,2) + t22 * t44 + t24 * t32 + t33) * g(2) + (t5 * rSges(6,1) - t4 * rSges(6,2) + t44 * t24 + (t16 * t34 + t31) * t22) * g(1)), (-m(3) - m(4) + t49) * (g(1) * t22 - g(2) * t24), -m(4) * g(3) * t29 - m(5) * (g(3) * (t27 + t37) + t30) - m(6) * (g(3) * (t15 * t28 + t32 + t37) + t30) + t50 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * t28) * t16 + (m(4) * rSges(4,1) - m(5) * (rSges(5,2) - pkin(3)) - m(6) * t34) * t15), t49 * (t15 * t50 - t46), -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5) + (-rSges(6,1) * t23 + rSges(6,2) * t21) * t46)];
taug = t1(:);
