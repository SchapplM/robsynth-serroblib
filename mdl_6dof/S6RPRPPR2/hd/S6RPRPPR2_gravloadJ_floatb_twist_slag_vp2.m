% Calculate Gravitation load on the joints for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:17
% DurationCPUTime: 0.55s
% Computational Cost: add. (342->84), mult. (302->89), div. (0->0), fcn. (247->10), ass. (0->46)
t56 = m(6) + m(7);
t43 = m(5) + t56;
t70 = mrSges(5,1) - mrSges(6,2);
t69 = -mrSges(5,2) + mrSges(6,3);
t20 = qJ(1) + pkin(9);
t14 = sin(t20);
t16 = cos(t20);
t60 = g(1) * t16 + g(2) * t14;
t19 = qJ(3) + pkin(10);
t13 = sin(t19);
t15 = cos(t19);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t66 = -t26 * mrSges(4,1) + t23 * mrSges(4,2) - t69 * t13 - t70 * t15;
t64 = m(3) + m(4);
t10 = t13 * qJ(5);
t47 = t16 * t15;
t63 = pkin(4) * t47 + t16 * t10;
t58 = -m(4) * pkin(2) - mrSges(3,1) + t66;
t57 = -m(4) * pkin(7) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t52 = g(3) * t15;
t11 = t15 * pkin(4);
t24 = sin(qJ(1));
t51 = t24 * pkin(1);
t17 = t26 * pkin(3);
t27 = cos(qJ(1));
t18 = t27 * pkin(1);
t12 = t17 + pkin(2);
t50 = t16 * t12 + t18;
t22 = sin(qJ(6));
t49 = t14 * t22;
t25 = cos(qJ(6));
t48 = t14 * t25;
t46 = t16 * t22;
t45 = t16 * t25;
t42 = t11 + t10 + t17;
t41 = -t12 - t10;
t21 = -qJ(4) - pkin(7);
t40 = -t14 * t21 + t50;
t39 = m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3);
t33 = t22 * mrSges(7,1) + t25 * mrSges(7,2);
t4 = -t13 * t49 + t45;
t3 = t13 * t48 + t46;
t2 = t13 * t46 + t48;
t1 = t13 * t45 - t49;
t5 = [(-t27 * mrSges(2,1) + t24 * mrSges(2,2) - m(5) * t40 - m(6) * (t40 + t63) - m(7) * (pkin(8) * t47 + t50 + t63) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t47 - t64 * t18 + t58 * t16 + (-m(7) * (pkin(5) - t21) + t57) * t14) * g(2) + (t24 * mrSges(2,1) - t4 * mrSges(7,1) + t27 * mrSges(2,2) + t3 * mrSges(7,2) + t64 * t51 - t43 * (-t16 * t21 - t51) + (-m(7) * pkin(5) + t57) * t16 + (m(5) * t12 - m(6) * (t41 - t11) - m(7) * t41 - t39 * t15 - t58) * t14) * g(1) (-t43 - t64) * g(3) (-m(5) * t17 - m(6) * t42 - m(7) * (t15 * pkin(8) + t42) - t15 * mrSges(7,3) - t33 * t13 + t66) * g(3) + (mrSges(4,2) * t26 + (m(6) * pkin(4) - t39 + t70) * t13 + (-qJ(5) * t56 - t33 - t69) * t15 + (t43 * pkin(3) + mrSges(4,1)) * t23) * t60 (-g(1) * t14 + g(2) * t16) * t43 (-t13 * t60 + t52) * t56, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - (-mrSges(7,1) * t25 + mrSges(7,2) * t22) * t52];
taug  = t5(:);
