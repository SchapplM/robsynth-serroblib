% Calculate Gravitation load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (99->63), mult. (185->89), div. (0->0), fcn. (159->6), ass. (0->27)
t36 = rSges(6,3) + pkin(7);
t13 = sin(qJ(1));
t16 = cos(qJ(1));
t35 = g(1) * t16 + g(2) * t13;
t11 = sin(qJ(5));
t14 = cos(qJ(5));
t34 = m(5) * rSges(5,1) + m(6) * (rSges(6,1) * t14 - rSges(6,2) * t11 + pkin(4));
t30 = -rSges(5,3) - pkin(6);
t29 = t16 * pkin(1) + t13 * qJ(2);
t28 = t11 * t16;
t27 = t13 * t11;
t26 = t13 * t14;
t25 = t14 * t16;
t24 = -pkin(1) - qJ(3);
t23 = -m(4) - m(5) - m(6);
t22 = t16 * qJ(3) + t29;
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t21 = rSges(5,1) * t12 + rSges(5,2) * t15;
t19 = m(5) * rSges(5,2) - m(6) * t36;
t18 = pkin(4) * t12 - t36 * t15;
t9 = t16 * qJ(2);
t4 = t12 * t25 - t27;
t3 = -t12 * t28 - t26;
t2 = -t12 * t26 - t28;
t1 = t12 * t27 - t25;
t5 = [-m(2) * (g(1) * (-t13 * rSges(2,1) - rSges(2,2) * t16) + g(2) * (rSges(2,1) * t16 - t13 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t16 + t9 + (rSges(3,2) - pkin(1)) * t13) + g(2) * (-rSges(3,2) * t16 + t13 * rSges(3,3) + t29)) - m(4) * (g(1) * (rSges(4,2) * t16 + t9) + g(2) * (rSges(4,3) * t16 + t22) + (g(1) * (-rSges(4,3) + t24) + g(2) * rSges(4,2)) * t13) - m(5) * (g(1) * t9 + g(2) * t22 + (g(1) * t30 + g(2) * t21) * t16 + (g(1) * (-t21 + t24) + g(2) * t30) * t13) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t9) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t22) + (-g(1) * pkin(6) + g(2) * t18) * t16 + (g(1) * (-t18 + t24) - g(2) * pkin(6)) * t13), (-m(3) + t23) * (g(1) * t13 - g(2) * t16), t23 * t35, (t34 * t12 + t19 * t15) * g(3) + t35 * (t19 * t12 - t34 * t15), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t11 - rSges(6,2) * t14) * t15)];
taug = t5(:);
