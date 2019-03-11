% Calculate Gravitation load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.09s
% Computational Cost: add. (95->27), mult. (54->33), div. (0->0), fcn. (34->4), ass. (0->15)
t21 = rSges(5,1) + pkin(3);
t20 = qJ(4) + rSges(5,3);
t14 = pkin(6) + qJ(2);
t11 = sin(t14);
t19 = pkin(2) * t11;
t13 = qJ(3) + t14;
t10 = cos(t13);
t9 = sin(t13);
t18 = t21 * t10 + t20 * t9;
t17 = t10 * rSges(4,1) - rSges(4,2) * t9;
t16 = -rSges(4,1) * t9 - rSges(4,2) * t10;
t15 = t20 * t10 - t21 * t9;
t12 = cos(t14);
t8 = pkin(2) * t12;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t11 - rSges(3,2) * t12) + g(2) * (rSges(3,1) * t12 - rSges(3,2) * t11)) - m(4) * (g(1) * (t16 - t19) + g(2) * (t17 + t8)) - m(5) * (g(1) * (t15 - t19) + g(2) * (t8 + t18)) -m(4) * (g(1) * t16 + g(2) * t17) - m(5) * (g(1) * t15 + g(2) * t18) -m(5) * (g(1) * t9 - g(2) * t10)];
taug  = t1(:);
