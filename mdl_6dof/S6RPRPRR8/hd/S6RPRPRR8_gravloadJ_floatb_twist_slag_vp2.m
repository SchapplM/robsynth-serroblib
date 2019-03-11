% Calculate Gravitation load on the joints for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:56
% EndTime: 2019-03-09 03:57:57
% DurationCPUTime: 0.65s
% Computational Cost: add. (311->90), mult. (389->104), div. (0->0), fcn. (338->10), ass. (0->52)
t48 = -m(5) - m(6) - m(7);
t84 = -mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t82 = g(1) * t30 - g(2) * t33;
t25 = qJ(3) + pkin(10);
t18 = sin(t25);
t19 = cos(t25);
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t80 = -t29 * mrSges(4,1) - t18 * mrSges(5,1) - t32 * mrSges(4,2) + t84 * t19;
t69 = m(7) * pkin(5);
t78 = -m(3) - m(4);
t76 = mrSges(6,1) + t69;
t31 = cos(qJ(5));
t17 = t31 * pkin(5) + pkin(4);
t26 = qJ(5) + qJ(6);
t20 = sin(t26);
t21 = cos(t26);
t28 = sin(qJ(5));
t75 = m(6) * pkin(4) + m(7) * t17 + t31 * mrSges(6,1) + t21 * mrSges(7,1) - t28 * mrSges(6,2) - t20 * mrSges(7,2);
t72 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t34 = -pkin(9) - pkin(8);
t58 = t19 * t34;
t62 = t19 * pkin(8);
t71 = mrSges(2,2) - mrSges(3,3) - m(6) * (t18 * pkin(4) - t62) - m(7) * (t18 * t17 + t58) + t80;
t52 = t33 * t21;
t57 = t30 * t20;
t5 = -t18 * t57 + t52;
t53 = t33 * t20;
t56 = t30 * t21;
t6 = t18 * t56 + t53;
t68 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t18 * t53 + t56;
t8 = t18 * t52 - t57;
t67 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t63 = g(3) * t19;
t61 = t29 * pkin(3);
t55 = t30 * t28;
t54 = t30 * t31;
t51 = t33 * t28;
t50 = t33 * t31;
t49 = t33 * pkin(1) + t30 * qJ(2);
t23 = t33 * qJ(2);
t47 = -t30 * pkin(1) + t23;
t42 = -mrSges(7,1) * t20 - mrSges(7,2) * t21;
t11 = t18 * t51 + t54;
t9 = -t18 * t55 + t50;
t27 = -qJ(4) - pkin(7);
t12 = t18 * t50 - t55;
t10 = t18 * t54 + t51;
t1 = [(-t51 * t69 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t78 * t49 + t48 * (-t33 * t27 + t30 * t61 + t49) + (-m(4) * pkin(7) - t72) * t33 + t71 * t30) * g(2) + (t55 * t69 - m(3) * t47 - m(4) * t23 - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + t48 * (t30 * t27 + t33 * t61 + t47) + (-m(4) * (-pkin(1) - pkin(7)) + t72) * t30 + t71 * t33) * g(1), -t82 * (-t48 - t78) (m(5) * t61 - m(6) * (-t61 + t62) - m(7) * (-t58 - t61) + t75 * t18 - t80) * g(3) + t82 * (mrSges(4,2) * t29 + (-mrSges(5,1) - t75) * t19 + (-m(6) * pkin(8) + m(7) * t34 - t84) * t18 + (t48 * pkin(3) - mrSges(4,1)) * t32) (g(1) * t33 + g(2) * t30) * t48 (mrSges(6,2) * t31 + t76 * t28 - t42) * t63 + (-t12 * mrSges(6,2) - t76 * t11 - t67) * g(2) + (t10 * mrSges(6,2) - t76 * t9 - t68) * g(1), -g(1) * t68 - g(2) * t67 - t42 * t63];
taug  = t1(:);
