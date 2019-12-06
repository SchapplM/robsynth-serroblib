% Calculate Gravitation load on the joints for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:59
% EndTime: 2019-12-05 15:03:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (109->37), mult. (131->53), div. (0->0), fcn. (108->6), ass. (0->20)
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t29 = g(1) * t11 + g(2) * t10;
t9 = pkin(8) + qJ(3);
t8 = cos(t9);
t28 = g(3) * t8;
t27 = -m(5) - m(6);
t7 = sin(t9);
t26 = t8 * pkin(3) + t7 * qJ(4);
t23 = rSges(6,3) + pkin(6);
t12 = sin(qJ(5));
t21 = t10 * t12;
t13 = cos(qJ(5));
t20 = t10 * t13;
t19 = t11 * t12;
t18 = t11 * t13;
t17 = -m(3) - m(4) + t27;
t16 = t29 * qJ(4) * t8;
t15 = rSges(6,1) * t12 + rSges(6,2) * t13;
t1 = [(-m(2) + t17) * g(3), t17 * (g(1) * t10 - g(2) * t11), -m(4) * g(3) * (rSges(4,1) * t8 - t7 * rSges(4,2)) - m(5) * (g(3) * (-rSges(5,2) * t8 + rSges(5,3) * t7 + t26) + t16) - m(6) * (g(3) * (t15 * t7 + t23 * t8 + t26) + t16) + t29 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * t15) * t8 + (m(4) * rSges(4,1) - m(5) * (rSges(5,2) - pkin(3)) - m(6) * (-pkin(3) - t23)) * t7), t27 * (t29 * t7 - t28), -m(6) * (g(1) * ((t7 * t18 - t21) * rSges(6,1) + (-t7 * t19 - t20) * rSges(6,2)) + g(2) * ((t7 * t20 + t19) * rSges(6,1) + (-t7 * t21 + t18) * rSges(6,2)) + (-rSges(6,1) * t13 + rSges(6,2) * t12) * t28)];
taug = t1(:);
