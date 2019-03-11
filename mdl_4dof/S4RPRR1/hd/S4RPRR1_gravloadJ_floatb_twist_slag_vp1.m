% Calculate Gravitation load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:44
% EndTime: 2019-03-08 18:31:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (110->35), mult. (66->44), div. (0->0), fcn. (42->8), ass. (0->22)
t15 = sin(qJ(1));
t25 = t15 * pkin(1);
t14 = qJ(1) + pkin(7);
t11 = cos(t14);
t16 = cos(qJ(1));
t13 = t16 * pkin(1);
t24 = pkin(2) * t11 + t13;
t12 = qJ(3) + t14;
t7 = sin(t12);
t8 = cos(t12);
t23 = t8 * rSges(4,1) - t7 * rSges(4,2);
t9 = qJ(4) + t12;
t4 = sin(t9);
t5 = cos(t9);
t22 = t5 * rSges(5,1) - t4 * rSges(5,2);
t21 = pkin(3) * t8 + t22;
t20 = -t7 * rSges(4,1) - t8 * rSges(4,2);
t19 = -t4 * rSges(5,1) - t5 * rSges(5,2);
t10 = sin(t14);
t18 = -pkin(2) * t10 - t25;
t17 = -pkin(3) * t7 + t19;
t1 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - t16 * rSges(2,2)) + g(2) * (t16 * rSges(2,1) - t15 * rSges(2,2))) - m(3) * (g(1) * (-t10 * rSges(3,1) - t11 * rSges(3,2) - t25) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t13)) - m(4) * (g(1) * (t18 + t20) + g(2) * (t23 + t24)) - m(5) * (g(1) * (t17 + t18) + g(2) * (t21 + t24)) (-m(3) - m(4) - m(5)) * g(3), -m(4) * (g(1) * t20 + g(2) * t23) - m(5) * (g(1) * t17 + g(2) * t21) -m(5) * (g(1) * t19 + g(2) * t22)];
taug  = t1(:);
