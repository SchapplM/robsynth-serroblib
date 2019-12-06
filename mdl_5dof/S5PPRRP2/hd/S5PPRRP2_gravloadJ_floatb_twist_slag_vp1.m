% Calculate Gravitation load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:08:34
% DurationCPUTime: 0.30s
% Computational Cost: add. (149->44), mult. (196->69), div. (0->0), fcn. (180->6), ass. (0->24)
t36 = -m(5) - m(6);
t15 = sin(pkin(7));
t16 = cos(pkin(7));
t34 = g(1) * t16 + g(2) * t15;
t35 = rSges(6,1) + pkin(4);
t24 = rSges(6,3) + qJ(5);
t14 = pkin(8) + qJ(3);
t12 = sin(t14);
t30 = g(3) * t12;
t17 = sin(qJ(4));
t28 = t15 * t17;
t18 = cos(qJ(4));
t27 = t15 * t18;
t26 = t16 * t17;
t25 = t16 * t18;
t23 = -m(3) - m(4) + t36;
t22 = rSges(5,1) * t18 - rSges(5,2) * t17;
t20 = t24 * t17 + t35 * t18;
t13 = cos(t14);
t4 = t13 * t25 + t28;
t3 = t13 * t26 - t27;
t2 = t13 * t27 - t26;
t1 = t13 * t28 + t25;
t5 = [(-m(2) + t23) * g(3), t23 * (g(1) * t15 - g(2) * t16), (-m(4) * (g(3) * rSges(4,1) - t34 * rSges(4,2)) - m(5) * (t34 * rSges(5,3) + g(3) * t22) - m(6) * (t34 * rSges(6,2) + g(3) * t20)) * t13 + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * g(3) + t34 * (m(4) * rSges(4,1) - m(5) * (-pkin(3) - t22) - m(6) * (-pkin(3) - t20))) * t12 + t36 * (g(3) * (t13 * pkin(3) + t12 * pkin(6)) + t34 * pkin(6) * t13), -m(5) * (g(1) * (-rSges(5,1) * t3 - t4 * rSges(5,2)) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (t24 * t4 - t3 * t35) + g(2) * (-t1 * t35 + t24 * t2)) + (-m(5) * (-rSges(5,1) * t17 - rSges(5,2) * t18) - m(6) * (-t35 * t17 + t24 * t18)) * t30, -m(6) * (g(1) * t3 + g(2) * t1 + t17 * t30)];
taug = t5(:);
