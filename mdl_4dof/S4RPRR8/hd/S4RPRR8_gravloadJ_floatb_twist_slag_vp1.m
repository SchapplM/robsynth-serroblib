% Calculate Gravitation load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (80->47), mult. (116->63), div. (0->0), fcn. (89->6), ass. (0->22)
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t28 = g(1) * t11 - g(2) * t13;
t27 = t13 * pkin(1) + t11 * qJ(2);
t9 = qJ(3) + qJ(4);
t5 = cos(t9);
t26 = rSges(5,1) * t5;
t4 = sin(t9);
t25 = rSges(5,2) * t4;
t12 = cos(qJ(3));
t24 = pkin(3) * t12;
t21 = rSges(4,3) + pkin(5);
t20 = rSges(5,3) + pkin(6) + pkin(5);
t19 = -rSges(5,1) * t4 - rSges(5,2) * t5;
t10 = sin(qJ(3));
t17 = rSges(4,1) * t10 + rSges(4,2) * t12;
t7 = t13 * qJ(2);
t16 = g(1) * t7 + g(2) * t27;
t15 = pkin(3) * t10 - t19;
t3 = t13 * t25;
t2 = t11 * t26;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t11 - rSges(2,2) * t13) + g(2) * (rSges(2,1) * t13 - rSges(2,2) * t11)) - m(3) * (g(1) * (rSges(3,3) * t13 + t7 + (rSges(3,2) - pkin(1)) * t11) + g(2) * (-rSges(3,2) * t13 + rSges(3,3) * t11 + t27)) - m(4) * ((g(1) * t17 + g(2) * t21) * t13 + (g(1) * (-pkin(1) - t21) + g(2) * t17) * t11 + t16) - m(5) * ((g(1) * t15 + g(2) * t20) * t13 + (g(1) * (-pkin(1) - t20) + g(2) * t15) * t11 + t16), (-m(3) - m(4) - m(5)) * t28, -m(4) * (-g(3) * t17 + t28 * (rSges(4,1) * t12 - rSges(4,2) * t10)) - m(5) * (g(1) * (t2 + (t24 - t25) * t11) + g(2) * (t3 + (-t24 - t26) * t13) - g(3) * t15), -m(5) * (g(1) * (-t11 * t25 + t2) + g(2) * (-t13 * t26 + t3) + g(3) * t19)];
taug = t1(:);
