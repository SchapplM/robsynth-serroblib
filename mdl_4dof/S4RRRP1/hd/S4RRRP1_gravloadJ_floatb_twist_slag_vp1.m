% Calculate Gravitation load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:53
% EndTime: 2019-03-08 18:35:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (130->37), mult. (84->48), div. (0->0), fcn. (56->6), ass. (0->24)
t27 = rSges(5,1) + pkin(3);
t12 = qJ(1) + qJ(2);
t8 = sin(t12);
t26 = pkin(2) * t8;
t13 = sin(qJ(1));
t25 = t13 * pkin(1);
t9 = cos(t12);
t24 = t9 * rSges(3,1) - rSges(3,2) * t8;
t10 = qJ(3) + t12;
t6 = sin(t10);
t7 = cos(t10);
t23 = t7 * rSges(4,1) - rSges(4,2) * t6;
t5 = pkin(2) * t9;
t22 = t23 + t5;
t21 = -t6 * rSges(5,2) + t27 * t7;
t20 = -rSges(3,1) * t8 - rSges(3,2) * t9;
t19 = -rSges(4,1) * t6 - rSges(4,2) * t7;
t18 = t21 + t5;
t17 = -t7 * rSges(5,2) - t27 * t6;
t16 = t19 - t26;
t15 = t17 - t26;
t14 = cos(qJ(1));
t11 = t14 * pkin(1);
t1 = [-m(2) * (g(1) * (-t13 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t13 * rSges(2,2))) - m(3) * (g(1) * (t20 - t25) + g(2) * (t11 + t24)) - m(4) * (g(1) * (t16 - t25) + g(2) * (t11 + t22)) - m(5) * (g(1) * (t15 - t25) + g(2) * (t11 + t18)) -m(3) * (g(1) * t20 + g(2) * t24) - m(4) * (g(1) * t16 + g(2) * t22) - m(5) * (g(1) * t15 + g(2) * t18) -m(4) * (g(1) * t19 + g(2) * t23) - m(5) * (g(1) * t17 + g(2) * t21) -m(5) * g(3)];
taug  = t1(:);
