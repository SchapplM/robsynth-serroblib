% Calculate Gravitation load on the joints for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:17
% EndTime: 2019-03-09 03:21:19
% DurationCPUTime: 0.59s
% Computational Cost: add. (267->73), mult. (374->84), div. (0->0), fcn. (319->8), ass. (0->39)
t36 = -m(5) - m(6) - m(7);
t69 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t56 = -g(1) * t20 + g(2) * t23;
t67 = mrSges(6,1) + mrSges(7,1);
t60 = mrSges(6,2) + mrSges(7,2);
t15 = qJ(3) + pkin(9);
t10 = sin(t15);
t11 = cos(t15);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t65 = -t19 * mrSges(4,1) - t10 * mrSges(5,1) - t22 * mrSges(4,2) + t69 * t11;
t51 = m(7) * pkin(5);
t63 = -m(3) - m(4);
t18 = sin(qJ(5));
t21 = cos(qJ(5));
t9 = t21 * pkin(5) + pkin(4);
t59 = m(6) * pkin(4) + m(7) * t9 - t60 * t18 + t67 * t21;
t55 = t51 + t67;
t54 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t16 = -qJ(6) - pkin(8);
t42 = t11 * t16;
t46 = t11 * pkin(8);
t53 = mrSges(2,2) - mrSges(3,3) - m(6) * (t10 * pkin(4) - t46) - m(7) * (t10 * t9 + t42) + t65;
t45 = t19 * pkin(3);
t41 = t20 * t18;
t40 = t20 * t21;
t39 = t23 * t18;
t38 = t23 * t21;
t37 = t23 * pkin(1) + t20 * qJ(2);
t13 = t23 * qJ(2);
t35 = -t20 * pkin(1) + t13;
t3 = t10 * t39 + t40;
t1 = -t10 * t41 + t38;
t17 = -qJ(4) - pkin(7);
t4 = t10 * t38 - t41;
t2 = t10 * t40 + t39;
t5 = [(-t39 * t51 + t63 * t37 + t36 * (-t23 * t17 + t20 * t45 + t37) - t67 * t2 - t60 * t1 + (-m(4) * pkin(7) - t54) * t23 + t53 * t20) * g(2) + (t41 * t51 - m(3) * t35 - m(4) * t13 - t67 * t4 + t36 * (t20 * t17 + t23 * t45 + t35) + t60 * t3 + (-m(4) * (-pkin(1) - pkin(7)) + t54) * t20 + t53 * t23) * g(1), t56 * (-t36 - t63) (m(5) * t45 - m(6) * (-t45 + t46) - m(7) * (-t42 - t45) + t59 * t10 - t65) * g(3) - t56 * (mrSges(4,2) * t19 + (-mrSges(5,1) - t59) * t11 + (-m(6) * pkin(8) + m(7) * t16 - t69) * t10 + (t36 * pkin(3) - mrSges(4,1)) * t22) (g(1) * t23 + g(2) * t20) * t36 (t55 * t18 + t60 * t21) * g(3) * t11 + (-t55 * t3 - t60 * t4) * g(2) + (-t55 * t1 + t60 * t2) * g(1) (-g(3) * t10 - t56 * t11) * m(7)];
taug  = t5(:);
