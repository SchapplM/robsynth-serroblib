% Calculate Gravitation load on the joints for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:07
% EndTime: 2019-03-09 02:44:08
% DurationCPUTime: 0.56s
% Computational Cost: add. (293->82), mult. (332->90), div. (0->0), fcn. (276->8), ass. (0->41)
t72 = -mrSges(4,1) - mrSges(5,1);
t71 = mrSges(6,2) - mrSges(7,3);
t70 = mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t57 = -pkin(3) - pkin(4);
t69 = -m(6) * t57 - m(7) * (-pkin(8) + t57) - t71;
t19 = qJ(1) + pkin(9);
t13 = sin(t19);
t14 = cos(t19);
t62 = g(1) * t14 + g(2) * t13;
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t66 = t70 * t21 - t72 * t24;
t64 = m(6) + m(7);
t44 = m(5) + t64;
t60 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2) + mrSges(6,3);
t59 = t71 * t24 - t66;
t54 = g(3) * t24;
t22 = sin(qJ(1));
t53 = t22 * pkin(1);
t17 = t24 * pkin(3);
t52 = t24 * pkin(8);
t25 = cos(qJ(1));
t18 = t25 * pkin(1);
t51 = t14 * t24;
t20 = sin(qJ(6));
t49 = t21 * t20;
t23 = cos(qJ(6));
t48 = t21 * t23;
t15 = t21 * qJ(4);
t46 = t17 + t15;
t43 = t14 * pkin(2) + t13 * pkin(7) + t18;
t42 = t24 * pkin(4) + t46;
t41 = t14 * pkin(7) - t53;
t40 = -pkin(2) - t15;
t39 = pkin(3) * t51 + t14 * t15 + t43;
t30 = m(7) * pkin(5) + t23 * mrSges(7,1) - t20 * mrSges(7,2);
t4 = -t13 * t20 + t14 * t48;
t3 = -t13 * t23 - t14 * t49;
t2 = -t13 * t48 - t14 * t20;
t1 = t13 * t49 - t14 * t23;
t5 = [(-m(3) * t18 - m(4) * t43 - m(5) * t39 - t25 * mrSges(2,1) - t4 * mrSges(7,1) + t22 * mrSges(2,2) - t3 * mrSges(7,2) - t64 * (pkin(4) * t51 - t13 * qJ(5) + t39) + t60 * t13 + (-mrSges(3,1) - m(7) * (t21 * pkin(5) + t52) + t59) * t14) * g(2) + (m(3) * t53 + t22 * mrSges(2,1) - t2 * mrSges(7,1) + t25 * mrSges(2,2) - t1 * mrSges(7,2) + (-m(4) - m(5)) * t41 - t64 * (-t14 * qJ(5) + t41) + t60 * t14 + (mrSges(3,1) + m(4) * pkin(2) - m(5) * (t40 - t17) - m(6) * t40 - m(7) * (-pkin(2) + (-pkin(5) - qJ(4)) * t21) + t69 * t24 + t66) * t13) * g(1) (-m(3) - m(4) - t44) * g(3) (-m(5) * t46 - m(6) * t42 - m(7) * (t42 + t52) - t30 * t21 + t59) * g(3) + ((m(5) * pkin(3) + t69 - t72) * t21 + (-qJ(4) * t44 - t30 - t70) * t24) * t62 (-t62 * t21 + t54) * t44, t64 * (g(1) * t13 - g(2) * t14) -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - (mrSges(7,1) * t20 + mrSges(7,2) * t23) * t54];
taug  = t5(:);
