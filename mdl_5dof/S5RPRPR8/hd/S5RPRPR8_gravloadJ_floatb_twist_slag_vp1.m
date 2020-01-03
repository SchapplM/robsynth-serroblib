% Calculate Gravitation load on the joints for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:18
% EndTime: 2019-12-31 18:21:19
% DurationCPUTime: 0.42s
% Computational Cost: add. (238->72), mult. (242->102), div. (0->0), fcn. (217->10), ass. (0->36)
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t17 = sin(pkin(9));
t18 = cos(pkin(9));
t28 = rSges(5,1) * t18 - rSges(5,2) * t17 + pkin(3);
t32 = rSges(5,3) + qJ(4);
t49 = t20 * t32 + t28 * t22;
t46 = rSges(6,3) + pkin(7) + qJ(4);
t45 = t22 * rSges(4,1) - t20 * rSges(4,2);
t16 = qJ(1) + pkin(8);
t11 = sin(t16);
t13 = cos(t16);
t44 = g(1) * t13 + g(2) * t11;
t15 = pkin(9) + qJ(5);
t10 = sin(t15);
t12 = cos(t15);
t9 = t18 * pkin(4) + pkin(3);
t43 = m(5) * t28 + m(6) * (rSges(6,1) * t12 - rSges(6,2) * t10 + t9) + m(4) * rSges(4,1);
t42 = -m(5) - m(6);
t40 = pkin(4) * t17;
t21 = sin(qJ(1));
t37 = t21 * pkin(1);
t36 = t11 * t22;
t35 = t13 * t22;
t23 = cos(qJ(1));
t14 = t23 * pkin(1);
t31 = t13 * pkin(2) + t11 * pkin(6) + t14;
t30 = t13 * pkin(6) - t37;
t29 = t17 * rSges(5,1) + t18 * rSges(5,2);
t26 = t46 * t20 + t22 * t9;
t25 = m(4) * rSges(4,2) - m(5) * t32 - m(6) * t46;
t4 = t11 * t10 + t12 * t35;
t3 = -t10 * t35 + t11 * t12;
t2 = t13 * t10 - t12 * t36;
t1 = t10 * t36 + t13 * t12;
t5 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (-t11 * rSges(3,1) - t13 * rSges(3,2) - t37) + g(2) * (t13 * rSges(3,1) - t11 * rSges(3,2) + t14)) - m(4) * (g(1) * (t13 * rSges(4,3) + t30) + g(2) * (t45 * t13 + t31) + (g(1) * (-pkin(2) - t45) + g(2) * rSges(4,3)) * t11) - m(5) * (g(1) * t30 + g(2) * t31 + (g(1) * t29 + t49 * g(2)) * t13 + (g(2) * t29 + (-pkin(2) - t49) * g(1)) * t11) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t30) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t31) + (g(1) * t40 + g(2) * t26) * t13 + (g(1) * (-pkin(2) - t26) + g(2) * t40) * t11), (-m(3) - m(4) + t42) * g(3), (t25 * t20 - t43 * t22) * g(3) + t44 * (t43 * t20 + t25 * t22), t42 * (-g(3) * t22 + t44 * t20), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (-rSges(6,1) * t10 - rSges(6,2) * t12) * t20)];
taug = t5(:);
