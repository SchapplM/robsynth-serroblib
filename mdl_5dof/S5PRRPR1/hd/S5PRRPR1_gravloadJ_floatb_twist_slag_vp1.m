% Calculate Gravitation load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:35
% DurationCPUTime: 0.18s
% Computational Cost: add. (213->45), mult. (124->55), div. (0->0), fcn. (92->8), ass. (0->25)
t46 = pkin(7) + qJ(4) + rSges(6,3);
t24 = pkin(9) + qJ(5);
t19 = sin(t24);
t21 = cos(t24);
t45 = rSges(6,1) * t21 - rSges(6,2) * t19;
t44 = rSges(5,3) + qJ(4);
t27 = cos(pkin(9));
t43 = rSges(5,2) * sin(pkin(9)) - pkin(3) - rSges(5,1) * t27;
t42 = -pkin(4) * t27 - pkin(3) - t45;
t41 = -m(5) - m(6);
t25 = pkin(8) + qJ(2);
t20 = sin(t25);
t40 = pkin(2) * t20;
t23 = qJ(3) + t25;
t16 = sin(t23);
t17 = cos(t23);
t35 = t17 * rSges(4,1) - rSges(4,2) * t16;
t34 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t32 = t44 * t16 - t43 * t17;
t31 = t43 * t16 + t44 * t17;
t30 = t46 * t16 - t42 * t17;
t29 = t42 * t16 + t46 * t17;
t22 = cos(t25);
t15 = pkin(2) * t22;
t1 = [(-m(2) - m(3) - m(4) + t41) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t20 - rSges(3,2) * t22) + g(2) * (rSges(3,1) * t22 - rSges(3,2) * t20)) - m(4) * (g(1) * (t34 - t40) + g(2) * (t15 + t35)) - m(5) * (g(1) * (t31 - t40) + g(2) * (t15 + t32)) - m(6) * (g(1) * (t29 - t40) + g(2) * (t15 + t30)), -m(4) * (g(1) * t34 + g(2) * t35) - m(5) * (g(1) * t31 + g(2) * t32) - m(6) * (g(1) * t29 + g(2) * t30), t41 * (g(1) * t16 - g(2) * t17), -m(6) * (g(3) * t45 + (g(1) * t17 + g(2) * t16) * (-rSges(6,1) * t19 - rSges(6,2) * t21))];
taug = t1(:);
