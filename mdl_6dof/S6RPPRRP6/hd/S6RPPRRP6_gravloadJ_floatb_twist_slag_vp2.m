% Calculate Gravitation load on the joints for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:12
% EndTime: 2019-03-09 02:10:13
% DurationCPUTime: 0.62s
% Computational Cost: add. (177->73), mult. (371->94), div. (0->0), fcn. (330->6), ass. (0->41)
t67 = mrSges(6,1) + mrSges(7,1);
t66 = -mrSges(6,2) + mrSges(7,3);
t65 = mrSges(6,3) + mrSges(7,2);
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t64 = t66 * t18 + t67 * t21;
t22 = cos(qJ(4));
t63 = t65 * t22;
t62 = -m(4) - m(5);
t61 = -m(6) - m(7);
t19 = sin(qJ(4));
t28 = pkin(5) * t21 + qJ(6) * t18;
t60 = (-m(7) * t28 - t64) * t22 - t65 * t19;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t59 = -g(1) * t23 - g(2) * t20;
t31 = t19 * mrSges(5,1) + t22 * mrSges(5,2);
t58 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - t31;
t57 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t56 = m(7) * pkin(5) + t67;
t55 = m(7) * qJ(6) + t66;
t50 = g(3) * t22;
t49 = t19 * t23;
t48 = t20 * t18;
t47 = t20 * t21;
t46 = t20 * t22;
t43 = t22 * t23;
t42 = t23 * t18;
t41 = t23 * t21;
t40 = -pkin(1) - qJ(3);
t39 = t23 * pkin(1) + t20 * qJ(2);
t38 = t23 * qJ(3) + t39;
t37 = t61 + t62;
t15 = t23 * qJ(2);
t36 = -t23 * pkin(7) + t15;
t27 = -t20 * pkin(7) + t38;
t4 = t19 * t41 - t48;
t3 = t19 * t42 + t47;
t2 = t19 * t47 + t42;
t1 = t19 * t48 - t41;
t5 = [(-m(3) * t39 - m(4) * t38 - m(5) * t27 + t65 * t43 + t61 * (pkin(4) * t49 - pkin(8) * t43 + t27) - t56 * t4 - t55 * t3 + t58 * t23 + t57 * t20) * g(2) + (-m(5) * t36 + t61 * (pkin(8) * t46 + t36) + t56 * t2 + (-m(3) - m(4)) * t15 + t55 * t1 + t57 * t23 + (m(3) * pkin(1) + t62 * t40 + t61 * (-t19 * pkin(4) + t40) - t58 - t63) * t20) * g(1) (-g(1) * t20 + g(2) * t23) * (m(3) - t37) -t59 * t37, t59 * (mrSges(5,1) * t22 - mrSges(5,2) * t19) + (t61 * (t19 * t20 * pkin(8) + pkin(4) * t46) + t60 * t20) * g(2) + (t61 * (pkin(4) * t43 + pkin(8) * t49) + t60 * t23) * g(1) + (t31 + (m(6) * pkin(4) - m(7) * (-pkin(4) - t28) + t64) * t19 + t61 * t22 * pkin(8) - t63) * g(3) (t56 * t18 - t55 * t21) * t50 + (t56 * t1 - t55 * t2) * g(2) + (t56 * t3 - t55 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t18 * t50) * m(7)];
taug  = t5(:);
