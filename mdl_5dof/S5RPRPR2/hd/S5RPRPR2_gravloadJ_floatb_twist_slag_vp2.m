% Calculate Gravitation load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (223->37), mult. (143->37), div. (0->0), fcn. (104->10), ass. (0->23)
t18 = pkin(9) + qJ(5);
t13 = sin(t18);
t15 = cos(t18);
t42 = mrSges(6,1) * t15 - mrSges(6,2) * t13;
t21 = cos(pkin(9));
t41 = -sin(pkin(9)) * mrSges(5,2) + m(5) * pkin(3) + m(6) * (pkin(4) * t21 + pkin(3)) + t21 * mrSges(5,1) + mrSges(4,1) + t42;
t40 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3) - m(6) * (-pkin(7) - qJ(4));
t39 = m(5) + m(6);
t33 = m(4) + t39;
t19 = qJ(1) + pkin(8);
t32 = m(3) + t33;
t29 = pkin(2) * t33 + mrSges(3,1);
t28 = pkin(1) * t32 + mrSges(2,1);
t17 = qJ(3) + t19;
t10 = sin(t17);
t11 = cos(t17);
t26 = t41 * t10 - t40 * t11;
t25 = t40 * t10 + t41 * t11;
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t16 = cos(t19);
t14 = sin(t19);
t1 = [(mrSges(2,2) * t24 + mrSges(3,2) * t16 + t14 * t29 + t23 * t28 + t26) * g(3) + (-t23 * mrSges(2,2) - t14 * mrSges(3,2) + t16 * t29 + t24 * t28 + t25) * g(2), -t32 * g(1), g(2) * t25 + g(3) * t26, t39 * (-g(2) * t11 - g(3) * t10), -g(1) * t42 + (-g(2) * t10 + g(3) * t11) * (mrSges(6,1) * t13 + mrSges(6,2) * t15)];
taug = t1(:);
