% Calculate Gravitation load on the joints for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:16
% EndTime: 2019-12-05 18:11:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (253->54), mult. (211->51), div. (0->0), fcn. (158->10), ass. (0->28)
t20 = pkin(9) + qJ(3);
t14 = sin(t20);
t16 = qJ(4) + t20;
t10 = sin(t16);
t11 = cos(t16);
t13 = qJ(5) + t16;
t7 = sin(t13);
t8 = cos(t13);
t39 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t33 = -t11 * mrSges(5,1) + mrSges(5,2) * t10 - t39;
t51 = mrSges(4,2) * t14 + t33;
t24 = sin(qJ(1));
t25 = cos(qJ(1));
t49 = g(1) * t25 + g(2) * t24;
t22 = cos(pkin(9));
t12 = t22 * pkin(2) + pkin(1);
t15 = cos(t20);
t9 = pkin(3) * t15;
t3 = t9 + t12;
t6 = pkin(4) * t11;
t48 = m(4) * t12 + m(5) * t3 + mrSges(4,1) * t15 + mrSges(2,1) + m(3) * pkin(1) + t22 * mrSges(3,1) - sin(pkin(9)) * mrSges(3,2) + m(6) * (t6 + t3) - t51;
t23 = -pkin(6) - qJ(2);
t19 = -pkin(7) + t23;
t47 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) + m(6) * (-pkin(8) + t19) - mrSges(6,3) + m(5) * t19 - mrSges(5,3) + m(4) * t23 - mrSges(4,3);
t40 = m(5) * pkin(3) + mrSges(4,1);
t34 = mrSges(6,1) * t7 + mrSges(6,2) * t8;
t29 = mrSges(5,2) * t11 + t34;
t1 = [(t47 * t24 - t48 * t25) * g(2) + (t48 * t24 + t47 * t25) * g(1), (-g(1) * t24 + g(2) * t25) * (m(3) + m(4) + m(5) + m(6)), (-m(6) * (t6 + t9) - t40 * t15 + t51) * g(3) + t49 * (-m(6) * (-pkin(3) * t14 - pkin(4) * t10) + mrSges(5,1) * t10 + mrSges(4,2) * t15 + t40 * t14 + t29), (-m(6) * t6 + t33) * g(3) + t49 * ((m(6) * pkin(4) + mrSges(5,1)) * t10 + t29), -g(3) * t39 + t49 * t34];
taug = t1(:);
