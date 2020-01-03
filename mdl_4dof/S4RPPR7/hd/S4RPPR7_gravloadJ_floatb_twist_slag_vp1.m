% Calculate Gravitation load on the joints for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (63->37), mult. (89->49), div. (0->0), fcn. (68->6), ass. (0->17)
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t25 = g(1) * t12 - g(2) * t13;
t24 = -m(4) - m(5);
t23 = t13 * pkin(1) + t12 * qJ(2);
t20 = rSges(5,3) + pkin(5) + qJ(3);
t19 = rSges(4,3) + qJ(3);
t8 = pkin(6) + qJ(4);
t3 = sin(t8);
t4 = cos(t8);
t17 = -t3 * rSges(5,1) - rSges(5,2) * t4;
t9 = sin(pkin(6));
t16 = rSges(4,1) * t9 + rSges(4,2) * cos(pkin(6));
t6 = t13 * qJ(2);
t15 = g(1) * t6 + g(2) * t23;
t14 = pkin(3) * t9 - t17;
t1 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t13) + g(2) * (rSges(2,1) * t13 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t13 + t6 + (rSges(3,2) - pkin(1)) * t12) + g(2) * (-rSges(3,2) * t13 + t12 * rSges(3,3) + t23)) - m(4) * ((g(1) * t16 + g(2) * t19) * t13 + (g(1) * (-pkin(1) - t19) + g(2) * t16) * t12 + t15) - m(5) * ((g(1) * t14 + g(2) * t20) * t13 + (g(1) * (-pkin(1) - t20) + g(2) * t14) * t12 + t15), (-m(3) + t24) * t25, t24 * (g(1) * t13 + g(2) * t12), -m(5) * (g(3) * t17 + t25 * (rSges(5,1) * t4 - rSges(5,2) * t3))];
taug = t1(:);
