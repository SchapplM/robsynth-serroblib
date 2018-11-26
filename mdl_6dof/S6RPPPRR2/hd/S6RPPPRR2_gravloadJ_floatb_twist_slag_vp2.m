% Calculate Gravitation load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:41
% EndTime: 2018-11-23 15:37:41
% DurationCPUTime: 0.31s
% Computational Cost: add. (259->66), mult. (228->77), div. (0->0), fcn. (183->10), ass. (0->36)
t17 = qJ(1) + pkin(9);
t12 = sin(t17);
t14 = cos(t17);
t49 = -g(1) * t12 + g(2) * t14;
t50 = -m(6) - m(7);
t48 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t16 = pkin(10) + qJ(5);
t11 = sin(t16);
t13 = cos(t16);
t18 = sin(pkin(10));
t30 = t11 * mrSges(6,1) + t13 * mrSges(6,2);
t47 = mrSges(3,2) - mrSges(4,3) - t18 * mrSges(5,1) - cos(pkin(10)) * mrSges(5,2) - t30 - m(7) * (pkin(5) * t11 - pkin(8) * t13) + t13 * mrSges(7,3);
t22 = sin(qJ(1));
t45 = pkin(1) * t22;
t44 = pkin(4) * t18;
t24 = cos(qJ(1));
t15 = t24 * pkin(1);
t21 = sin(qJ(6));
t41 = t12 * t21;
t23 = cos(qJ(6));
t40 = t12 * t23;
t39 = t14 * t21;
t38 = t14 * t23;
t37 = -m(5) + t50;
t36 = t14 * pkin(2) + t12 * qJ(3) + t15;
t35 = m(4) - t37;
t34 = m(7) * pkin(8) + mrSges(7,3);
t33 = t14 * qJ(3) - t45;
t29 = -pkin(2) * t12 + t33;
t26 = m(7) * pkin(5) + t23 * mrSges(7,1) - t21 * mrSges(7,2);
t20 = -pkin(7) - qJ(4);
t4 = t11 * t38 - t41;
t3 = t11 * t39 + t40;
t2 = t11 * t40 + t39;
t1 = -t11 * t41 + t38;
t5 = [(-m(3) * t15 - mrSges(2,1) * t24 - t2 * mrSges(7,1) + t22 * mrSges(2,2) - t1 * mrSges(7,2) + (-m(4) - m(5)) * t36 + t50 * (t12 * t44 - t14 * t20 + t36) + (-m(5) * qJ(4) - t48) * t14 + t47 * t12) * g(2) + (m(3) * t45 - m(4) * t29 - m(5) * t33 + t22 * mrSges(2,1) - t4 * mrSges(7,1) + mrSges(2,2) * t24 + t3 * mrSges(7,2) + t50 * (t12 * t20 + t14 * t44 + t29) + (-m(5) * (-pkin(2) - qJ(4)) + t48) * t12 + t47 * t14) * g(1) (-m(3) - t35) * g(3), t49 * t35 (g(1) * t14 + g(2) * t12) * t37 (t11 * t26 - t13 * t34 + t30) * g(3) + ((mrSges(6,1) + t26) * t13 + (-mrSges(6,2) + t34) * t11) * t49, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - g(3) * (-mrSges(7,1) * t21 - mrSges(7,2) * t23) * t13];
taug  = t5(:);
