% Calculate Gravitation load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (165->62), mult. (175->86), div. (0->0), fcn. (149->8), ass. (0->29)
t13 = qJ(1) + pkin(8);
t10 = sin(t13);
t11 = cos(t13);
t44 = -g(1) * t10 + g(2) * t11;
t15 = sin(qJ(4));
t18 = cos(qJ(4));
t31 = rSges(6,3) + pkin(7);
t42 = pkin(4) * t15 - t31 * t18;
t40 = t15 * rSges(5,1) + rSges(5,2) * t18;
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t20 = m(6) * (rSges(6,1) * t17 - rSges(6,2) * t14 + pkin(4)) + m(5) * rSges(5,1);
t39 = -m(5) * rSges(5,2) + m(6) * t31;
t38 = -pkin(2) - pkin(6);
t16 = sin(qJ(1));
t35 = pkin(1) * t16;
t29 = t14 * t15;
t27 = t15 * t17;
t26 = -m(4) - m(5) - m(6);
t19 = cos(qJ(1));
t12 = t19 * pkin(1);
t25 = t11 * pkin(2) + t10 * qJ(3) + t12;
t24 = t11 * qJ(3) - t35;
t23 = t11 * pkin(6) + t25;
t4 = -t10 * t14 + t11 * t27;
t3 = t10 * t17 + t11 * t29;
t2 = t10 * t27 + t11 * t14;
t1 = -t10 * t29 + t11 * t17;
t5 = [-m(2) * (g(1) * (-t16 * rSges(2,1) - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - t16 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t10 - rSges(3,2) * t11 - t35) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10 + t12)) - m(4) * (g(1) * (rSges(4,3) * t11 + (rSges(4,2) - pkin(2)) * t10 + t24) + g(2) * (-rSges(4,2) * t11 + rSges(4,3) * t10 + t25)) - m(5) * (g(1) * (t11 * t40 + t24) + g(2) * (rSges(5,3) * t11 + t23) + (g(1) * (-rSges(5,3) + t38) + g(2) * t40) * t10) - m(6) * ((rSges(6,1) * t2 + rSges(6,2) * t1 + t42 * t10 + t23) * g(2) + (rSges(6,1) * t4 - rSges(6,2) * t3 + t38 * t10 + t42 * t11 + t24) * g(1)), (-m(3) + t26) * g(3), -t26 * t44, (t20 * t15 - t18 * t39) * g(3) + t44 * (t15 * t39 + t20 * t18), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t14 - rSges(6,2) * t17) * t18)];
taug = t5(:);
