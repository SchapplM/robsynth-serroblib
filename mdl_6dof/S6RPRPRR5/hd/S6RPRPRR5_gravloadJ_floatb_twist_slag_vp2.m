% Calculate Gravitation load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:59
% EndTime: 2019-03-09 03:48:01
% DurationCPUTime: 0.82s
% Computational Cost: add. (454->110), mult. (567->130), div. (0->0), fcn. (563->10), ass. (0->56)
t85 = -m(6) - m(7);
t30 = sin(qJ(6));
t32 = cos(qJ(6));
t84 = m(7) * pkin(5) + t32 * mrSges(7,1) - t30 * mrSges(7,2) + mrSges(6,1);
t31 = sin(qJ(1));
t26 = pkin(10) + qJ(3);
t25 = cos(t26);
t60 = qJ(4) * t25;
t13 = t31 * t60;
t24 = sin(t26);
t70 = -pkin(3) - pkin(4);
t81 = t24 * t70;
t83 = t31 * t81 + t13;
t33 = cos(qJ(1));
t15 = t33 * t60;
t82 = t33 * t81 + t15;
t80 = (mrSges(4,1) + mrSges(5,1)) * t25 + (-mrSges(4,2) + mrSges(5,3)) * t24;
t56 = m(7) * pkin(9) + mrSges(7,3);
t28 = cos(pkin(10));
t79 = -mrSges(2,1) - m(3) * pkin(1) - t28 * mrSges(3,1) + sin(pkin(10)) * mrSges(3,2) - t80;
t78 = -mrSges(6,2) + t56;
t29 = -pkin(7) - qJ(2);
t77 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) + t85 * (-pkin(8) - t29);
t20 = t24 * qJ(4);
t62 = t25 * t33;
t76 = pkin(3) * t62 + t33 * t20;
t65 = sin(qJ(5));
t66 = cos(qJ(5));
t8 = t24 * t66 - t25 * t65;
t75 = g(1) * t33 + g(2) * t31;
t3 = t8 * t31;
t7 = t24 * t65 + t25 * t66;
t4 = t7 * t31;
t73 = -t4 * mrSges(6,2) + t84 * t3;
t72 = -t8 * mrSges(6,2) - t84 * t7;
t52 = t33 * t65;
t53 = t33 * t66;
t5 = -t24 * t52 - t25 * t53;
t6 = -t24 * t53 + t25 * t52;
t71 = t5 * mrSges(6,2) - t84 * t6;
t22 = t25 * pkin(3);
t61 = t22 + t20;
t59 = m(5) - t85;
t58 = t25 * pkin(4) + t61;
t23 = pkin(2) * t28 + pkin(1);
t16 = t33 * t23;
t50 = -t29 * t31 + t16;
t49 = -t23 - t20;
t48 = pkin(4) * t62 + t16 + t76;
t42 = t4 * t30 - t32 * t33;
t41 = -t30 * t33 - t4 * t32;
t35 = t25 * mrSges(5,3) + (-m(5) * pkin(3) - mrSges(5,1)) * t24;
t34 = (t25 * t70 + t49) * t31;
t2 = -t30 * t31 - t32 * t5;
t1 = t30 * t5 - t31 * t32;
t9 = [(-m(4) * t50 - m(5) * (t50 + t76) - m(6) * t48 + t5 * mrSges(6,1) - m(7) * (-pkin(5) * t5 + t48) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t78 * t6 + t79 * t33 + t77 * t31) * g(2) + (-m(6) * t34 + t4 * mrSges(6,1) - t41 * mrSges(7,1) - t42 * mrSges(7,2) - (-t4 * pkin(5) + t34) * m(7) - t78 * t3 + (m(4) * t23 - m(5) * (t49 - t22) - t79) * t31 + ((m(4) + m(5)) * t29 + t77) * t33) * g(1) (-g(1) * t31 + g(2) * t33) * (m(3) + m(4) + t59) t75 * (mrSges(4,1) * t24 + mrSges(4,2) * t25) + (-m(5) * t13 - t35 * t31 - m(6) * t83 - m(7) * (-pkin(9) * t4 + t83) + t4 * mrSges(7,3) + t73) * g(2) + (-m(5) * t15 - t35 * t33 - m(6) * t82 - m(7) * (t5 * pkin(9) + t82) - t5 * mrSges(7,3) + t71) * g(1) + (-m(5) * t61 - m(6) * t58 - m(7) * (-pkin(9) * t8 + t58) + t8 * mrSges(7,3) + t72 - t80) * g(3) (g(3) * t25 - t75 * t24) * t59 (-t56 * t8 - t72) * g(3) + (-t4 * t56 - t73) * g(2) + (t5 * t56 - t71) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) - g(3) * (-t30 * mrSges(7,1) - t32 * mrSges(7,2)) * t8];
taug  = t9(:);
