% Calculate Gravitation load on the joints for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:31
% EndTime: 2019-03-09 08:41:32
% DurationCPUTime: 0.80s
% Computational Cost: add. (331->104), mult. (543->116), div. (0->0), fcn. (491->8), ass. (0->51)
t75 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t74 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t81 = -m(6) - m(7);
t84 = -mrSges(7,2) - mrSges(6,3);
t27 = pkin(9) + qJ(5);
t20 = sin(t27);
t21 = cos(t27);
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t80 = t29 * mrSges(5,2);
t46 = t28 * mrSges(5,1) + t80;
t83 = -t75 * t20 + t74 * t21 - t46;
t82 = -m(4) - m(5);
t31 = sin(qJ(2));
t22 = t31 * qJ(3);
t33 = cos(qJ(2));
t61 = t33 * pkin(2) + t22;
t79 = (-mrSges(3,1) + mrSges(4,2)) * t33 + (mrSges(3,2) - mrSges(4,3)) * t31;
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t78 = g(1) * t34 + g(2) * t32;
t49 = m(5) * (-pkin(2) - qJ(4)) - mrSges(5,3);
t77 = (-mrSges(4,3) + t83) * t33 + (-mrSges(4,2) - t49 + (m(4) - t81) * pkin(2) - t84) * t31;
t76 = t84 * t33 + t79;
t73 = -m(5) * pkin(3) - mrSges(5,1) * t29 + mrSges(5,2) * t28 - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t72 = pkin(4) * t28;
t69 = g(3) * t33;
t30 = -pkin(8) - qJ(4);
t68 = t30 * t33;
t67 = t31 * t32;
t66 = t31 * t34;
t65 = t32 * t21;
t62 = t33 * t34;
t60 = t34 * pkin(1) + t32 * pkin(7);
t59 = qJ(3) * t33;
t58 = qJ(4) * t33;
t57 = -m(5) + t81;
t56 = t33 * t72;
t14 = t31 * t72;
t53 = t28 * t66;
t50 = pkin(2) * t62 + t34 * t22 + t60;
t42 = -pkin(1) - t61;
t25 = t34 * pkin(7);
t19 = pkin(4) * t29 + pkin(3);
t16 = t34 * t59;
t15 = t32 * t59;
t4 = -t20 * t67 + t21 * t34;
t3 = t20 * t34 + t31 * t65;
t2 = t20 * t66 + t65;
t1 = t20 * t32 - t21 * t66;
t5 = [(-t66 * t80 - m(3) * t60 - t53 * mrSges(5,1) + t82 * t50 + t81 * (pkin(4) * t53 + t32 * t19 - t30 * t62 + t50) - t75 * t2 - t74 * t1 + (-mrSges(5,3) + t84) * t62 + (-m(5) * t58 - mrSges(2,1) + t79) * t34 + t73 * t32) * g(2) + (t81 * (t34 * t19 + t32 * t68 + t25) - t75 * t4 - t74 * t3 + (-m(3) + t82) * t25 + t73 * t34 + (mrSges(2,1) - m(4) * t42 - t49 * t33 - (-m(5) * qJ(3) - t46) * t31 + t81 * (t42 - t14) + (m(3) + m(5)) * pkin(1) - t76) * t32) * g(1), t78 * (mrSges(3,1) * t31 + mrSges(3,2) * t33) + (t81 * (t30 * t67 + t32 * t56 + t15) + t82 * t15 + t77 * t32) * g(2) + (t81 * (t30 * t66 + t34 * t56 + t16) + t82 * t16 + t77 * t34) * g(1) + (-m(4) * t61 - m(5) * (t58 + t61) - t33 * mrSges(5,3) + t81 * (t14 + t61 - t68) + t83 * t31 + t76) * g(3) (-t78 * t31 + t69) * (m(4) - t57) (t31 * g(3) + t78 * t33) * t57 (t74 * t20 + t75 * t21) * t69 + (-t75 * t3 + t74 * t4) * g(2) + (t75 * t1 - t74 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t21 * t69) * m(7)];
taug  = t5(:);
