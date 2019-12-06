% Calculate Gravitation load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (275->48), mult. (207->48), div. (0->0), fcn. (157->10), ass. (0->32)
t28 = pkin(9) + qJ(4);
t23 = qJ(5) + t28;
t19 = cos(t23);
t12 = t19 * mrSges(6,1);
t31 = cos(pkin(9));
t20 = t31 * pkin(3) + pkin(2);
t22 = cos(t28);
t21 = sin(t28);
t43 = t21 * mrSges(5,2);
t56 = -t43 - sin(pkin(9)) * mrSges(4,2) + m(4) * pkin(2) + m(5) * t20 + m(6) * (pkin(4) * t22 + t20) + t31 * mrSges(4,1) + t22 * mrSges(5,1) + mrSges(3,1) + t12;
t32 = -pkin(7) - qJ(3);
t55 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - m(6) * (-pkin(8) + t32) - m(5) * t32;
t41 = m(6) * pkin(4) + mrSges(5,1);
t51 = -mrSges(5,2) * t22 - t41 * t21;
t29 = qJ(1) + qJ(2);
t24 = sin(t29);
t18 = sin(t23);
t44 = t18 * t24;
t46 = mrSges(6,2) * t19;
t50 = mrSges(6,1) * t44 + t24 * t46;
t25 = cos(t29);
t49 = g(3) * t25;
t45 = t18 * mrSges(6,2);
t42 = -m(4) - m(5) - m(6);
t40 = t12 - t45;
t39 = -mrSges(6,1) * t18 - t46;
t38 = mrSges(2,1) + (m(3) - t42) * pkin(1);
t36 = -mrSges(6,2) * t44 + t56 * t24 - t55 * t25;
t35 = (-t45 + t56) * t25 + t55 * t24;
t34 = cos(qJ(1));
t33 = sin(qJ(1));
t1 = [(t34 * mrSges(2,2) + t38 * t33 + t36) * g(3) + (-t33 * mrSges(2,2) + t38 * t34 + t35) * g(2), t35 * g(2) + t36 * g(3), (g(2) * t25 + g(3) * t24) * t42, (t51 * t24 - t50) * g(2) + (-t41 * t22 - t40 + t43) * g(1) + (-t39 - t51) * t49, -g(1) * t40 - g(2) * t50 - t39 * t49];
taug = t1(:);
