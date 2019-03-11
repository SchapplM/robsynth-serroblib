% Calculate Gravitation load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:05
% EndTime: 2019-03-09 02:50:06
% DurationCPUTime: 0.64s
% Computational Cost: add. (333->93), mult. (382->94), div. (0->0), fcn. (323->10), ass. (0->48)
t70 = m(6) + m(7);
t44 = m(5) + t70;
t79 = m(4) + t44;
t78 = mrSges(4,1) - mrSges(5,2);
t77 = -mrSges(4,2) + mrSges(5,3);
t18 = pkin(10) + qJ(6);
t14 = sin(t18);
t16 = cos(t18);
t20 = sin(pkin(10));
t22 = cos(pkin(10));
t34 = t20 * mrSges(6,1) + t22 * mrSges(6,2);
t59 = pkin(5) * t20;
t76 = -m(7) * t59 - t14 * mrSges(7,1) - t16 * mrSges(7,2) - t34;
t24 = -pkin(8) - qJ(5);
t75 = -m(6) * (-pkin(3) - qJ(5)) + mrSges(6,3) - m(7) * (-pkin(3) + t24) + mrSges(7,3);
t26 = sin(qJ(1));
t27 = cos(qJ(1));
t66 = g(1) * t27 + g(2) * t26;
t71 = m(5) + m(7);
t19 = pkin(9) + qJ(3);
t15 = sin(t19);
t10 = t15 * qJ(4);
t17 = cos(t19);
t50 = t17 * t27;
t69 = pkin(3) * t50 + t27 * t10;
t47 = t17 * pkin(3) + t10;
t68 = m(5) * t47;
t67 = t77 * t15 + t78 * t17;
t23 = cos(pkin(9));
t63 = -mrSges(2,1) - m(3) * pkin(1) - t23 * mrSges(3,1) + sin(pkin(9)) * mrSges(3,2) - t67;
t25 = -pkin(7) - qJ(2);
t61 = -m(6) * (pkin(4) - t25) - m(7) * (pkin(5) * t22 + pkin(4)) - t22 * mrSges(6,1) + t20 * mrSges(6,2) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * qJ(2) - mrSges(3,3);
t56 = g(3) * t17;
t54 = t15 * t27;
t53 = t16 * t27;
t52 = t17 * mrSges(7,3);
t51 = t17 * t24;
t49 = t26 * t14;
t48 = t26 * t16;
t45 = qJ(5) * t17;
t13 = pkin(2) * t23 + pkin(1);
t8 = t27 * t13;
t43 = -t26 * t25 + t8;
t4 = -t15 * t49 + t53;
t3 = t14 * t27 + t15 * t48;
t2 = t14 * t54 + t48;
t1 = t15 * t53 - t49;
t5 = [(-m(4) * t43 - m(6) * (t8 + t69) - mrSges(6,3) * t50 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t34 * t54 - t71 * (t43 + t69) + t61 * t26 + (-m(6) * t45 - m(7) * (t15 * t59 - t51) - t52 + t63) * t27) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(4) + t71) * t25 + t61) * t27 + (t68 + t75 * t17 + (m(6) * qJ(4) + t34 - m(7) * (-qJ(4) - t59)) * t15 + t79 * t13 - t63) * t26) * g(1) (-g(1) * t26 + g(2) * t27) * (m(3) + t79) (-t68 - m(6) * (t45 + t47) - t17 * mrSges(6,3) - m(7) * (t47 - t51) - t52 + t76 * t15 - t67) * g(3) + ((m(5) * pkin(3) + t75 + t78) * t15 + (-qJ(4) * t44 + t76 - t77) * t17) * t66 (-t66 * t15 + t56) * t44 (-t15 * g(3) - t66 * t17) * t70, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - (-mrSges(7,1) * t16 + mrSges(7,2) * t14) * t56];
taug  = t5(:);
