% Calculate Gravitation load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:19
% DurationCPUTime: 0.31s
% Computational Cost: add. (134->68), mult. (177->82), div. (0->0), fcn. (143->6), ass. (0->28)
t34 = rSges(6,3) + qJ(5);
t35 = rSges(6,1) + pkin(4);
t12 = qJ(3) + pkin(7);
t7 = sin(t12);
t8 = cos(t12);
t18 = -t34 * t8 + t35 * t7;
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t36 = m(4) * (rSges(4,1) * t16 - rSges(4,2) * t14);
t15 = sin(qJ(1));
t29 = g(1) * t15;
t30 = pkin(3) * t16;
t33 = t30 * t29;
t32 = -m(5) - m(6);
t31 = pkin(3) * t14;
t17 = cos(qJ(1));
t28 = g(2) * t17;
t27 = rSges(4,3) + pkin(6);
t26 = t17 * pkin(1) + t15 * qJ(2);
t25 = rSges(5,1) * t8 - rSges(5,2) * t7;
t24 = rSges(5,1) * t7 + rSges(5,2) * t8;
t23 = -t28 + t29;
t21 = rSges(4,1) * t14 + rSges(4,2) * t16;
t10 = t17 * qJ(2);
t13 = -qJ(4) - pkin(6);
t20 = g(1) * (t15 * t13 + t17 * t31 + t10) + g(2) * (t15 * t31 + t26);
t19 = t34 * t7 + t35 * t8;
t1 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - rSges(2,2) * t17) + g(2) * (rSges(2,1) * t17 - t15 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t17 + t10 + (rSges(3,2) - pkin(1)) * t15) + g(2) * (-rSges(3,2) * t17 + t15 * rSges(3,3) + t26)) - m(4) * (g(1) * t10 + g(2) * t26 + (g(1) * t21 + g(2) * t27) * t17 + (g(1) * (-pkin(1) - t27) + g(2) * t21) * t15) - m(5) * ((g(1) * t24 + g(2) * (rSges(5,3) - t13)) * t17 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t24) * t15 + t20) - m(6) * ((g(1) * t18 + g(2) * (rSges(6,2) - t13)) * t17 + (g(1) * (-rSges(6,2) - pkin(1)) + g(2) * t18) * t15 + t20), (-m(3) - m(4) + t32) * t23, m(4) * g(3) * t21 - m(5) * (t33 + g(3) * (-t24 - t31)) - m(6) * (t33 + g(3) * (-t18 - t31)) + (-m(5) * t25 - m(6) * t19 - t36) * t29 + (t36 - m(5) * (-t25 - t30) - m(6) * (-t19 - t30)) * t28, t32 * (g(1) * t17 + g(2) * t15), -m(6) * (g(3) * t7 - t23 * t8)];
taug = t1(:);
