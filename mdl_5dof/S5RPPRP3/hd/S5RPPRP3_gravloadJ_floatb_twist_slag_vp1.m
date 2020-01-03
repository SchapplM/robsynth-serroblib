% Calculate Gravitation load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:45
% DurationCPUTime: 0.24s
% Computational Cost: add. (131->49), mult. (124->63), div. (0->0), fcn. (93->6), ass. (0->20)
t8 = qJ(1) + pkin(7);
t5 = sin(t8);
t6 = cos(t8);
t30 = -g(1) * t5 + g(2) * t6;
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t29 = rSges(6,1) + pkin(4);
t15 = rSges(6,2) * t12 + t29 * t10;
t11 = sin(qJ(1));
t25 = pkin(1) * t11;
t24 = rSges(5,3) + pkin(6);
t23 = rSges(6,3) + qJ(5) + pkin(6);
t21 = -m(4) - m(5) - m(6);
t13 = cos(qJ(1));
t7 = t13 * pkin(1);
t20 = t6 * pkin(2) + t5 * qJ(3) + t7;
t19 = t6 * qJ(3) - t25;
t17 = rSges(5,1) * t10 + rSges(5,2) * t12;
t14 = g(1) * t19 + g(2) * t20;
t1 = [-m(2) * (g(1) * (-t11 * rSges(2,1) - rSges(2,2) * t13) + g(2) * (rSges(2,1) * t13 - t11 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t5 - rSges(3,2) * t6 - t25) + g(2) * (rSges(3,1) * t6 - rSges(3,2) * t5 + t7)) - m(4) * (g(1) * (rSges(4,3) * t6 + (rSges(4,2) - pkin(2)) * t5 + t19) + g(2) * (-rSges(4,2) * t6 + rSges(4,3) * t5 + t20)) - m(5) * ((g(1) * t17 + g(2) * t24) * t6 + (g(1) * (-pkin(2) - t24) + g(2) * t17) * t5 + t14) - m(6) * ((g(1) * t15 + g(2) * t23) * t6 + (g(1) * (-pkin(2) - t23) + g(2) * t15) * t5 + t14), (-m(3) + t21) * g(3), -t21 * t30, (m(5) * t17 + m(6) * t15) * g(3) + t30 * (m(5) * (rSges(5,1) * t12 - rSges(5,2) * t10) + m(6) * (-rSges(6,2) * t10 + t29 * t12)), -m(6) * (g(1) * t6 + g(2) * t5)];
taug = t1(:);
