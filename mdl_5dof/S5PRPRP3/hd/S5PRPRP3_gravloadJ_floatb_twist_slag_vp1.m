% Calculate Gravitation load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:39
% EndTime: 2019-12-05 15:32:41
% DurationCPUTime: 0.31s
% Computational Cost: add. (146->51), mult. (196->73), div. (0->0), fcn. (171->8), ass. (0->27)
t35 = rSges(6,1) + pkin(4);
t11 = sin(pkin(7));
t12 = cos(pkin(7));
t34 = g(1) * t12 + g(2) * t11;
t15 = sin(qJ(2));
t33 = pkin(2) * t15;
t30 = rSges(5,3) + pkin(6);
t14 = sin(qJ(4));
t29 = t11 * t14;
t16 = cos(qJ(4));
t28 = t11 * t16;
t27 = t12 * t14;
t26 = t12 * t16;
t25 = rSges(6,3) + qJ(5) + pkin(6);
t24 = -m(4) - m(5) - m(6);
t10 = qJ(2) + pkin(8);
t8 = cos(t10);
t3 = -t8 * t27 + t28;
t1 = -t8 * t29 - t26;
t22 = rSges(5,1) * t16 - rSges(5,2) * t14 + pkin(3);
t21 = -rSges(6,2) * t14 + t35 * t16 + pkin(3);
t17 = cos(qJ(2));
t9 = t17 * pkin(2);
t7 = sin(t10);
t4 = -t8 * t26 - t29;
t2 = -t8 * t28 + t27;
t5 = [(-m(2) - m(3) + t24) * g(3), -m(3) * (g(3) * (rSges(3,1) * t17 - t15 * rSges(3,2)) + t34 * (-rSges(3,1) * t15 - rSges(3,2) * t17)) - m(4) * (g(3) * (rSges(4,1) * t8 - rSges(4,2) * t7 + t9) + t34 * (-rSges(4,1) * t7 - rSges(4,2) * t8 - t33)) - m(5) * (g(3) * (t22 * t8 + t30 * t7 + t9) + t34 * (-t22 * t7 + t30 * t8 - t33)) - m(6) * (g(3) * (t21 * t8 + t25 * t7 + t9) + t34 * (-t21 * t7 + t25 * t8 - t33)), t24 * (g(1) * t11 - g(2) * t12), -m(5) * (g(1) * (rSges(5,1) * t3 + rSges(5,2) * t4) + g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2)) - m(6) * (g(1) * (rSges(6,2) * t4 + t35 * t3) + g(2) * (rSges(6,2) * t2 + t35 * t1)) + (-m(5) * (-rSges(5,1) * t14 - rSges(5,2) * t16) - m(6) * (-rSges(6,2) * t16 - t35 * t14)) * g(3) * t7, -m(6) * (-g(3) * t8 + t34 * t7)];
taug = t5(:);
