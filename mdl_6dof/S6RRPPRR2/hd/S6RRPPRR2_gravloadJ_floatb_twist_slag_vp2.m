% Calculate Gravitation load on the joints for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:51
% EndTime: 2019-03-09 08:49:52
% DurationCPUTime: 0.78s
% Computational Cost: add. (464->113), mult. (467->120), div. (0->0), fcn. (413->12), ass. (0->63)
t30 = sin(pkin(11));
t31 = cos(pkin(11));
t92 = -m(5) * pkin(3) - t31 * mrSges(5,1) + t30 * mrSges(5,2) - mrSges(4,1);
t82 = -m(6) - m(7);
t57 = m(5) - t82;
t91 = m(4) + t57;
t90 = mrSges(4,2) - mrSges(6,3) - mrSges(7,3);
t29 = qJ(2) + pkin(10);
t21 = sin(t29);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t85 = g(1) * t37 + g(2) * t35;
t89 = t85 * t21;
t18 = t31 * pkin(4) + pkin(3);
t28 = pkin(11) + qJ(5);
t22 = cos(t28);
t13 = pkin(5) * t22 + t18;
t24 = qJ(6) + t28;
t16 = sin(t24);
t17 = cos(t24);
t20 = sin(t28);
t88 = -m(6) * t18 - m(7) * t13 - t22 * mrSges(6,1) - t17 * mrSges(7,1) + t20 * mrSges(6,2) + t16 * mrSges(7,2);
t23 = cos(t29);
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t87 = -t36 * mrSges(3,1) + t34 * mrSges(3,2) + t90 * t21 + t92 * t23;
t86 = m(5) * qJ(4) + mrSges(5,3);
t84 = -m(3) * pkin(1) - mrSges(2,1) + t87;
t81 = m(7) * pkin(5) + mrSges(6,1);
t32 = -qJ(3) - pkin(7);
t80 = -m(3) * pkin(7) + m(5) * t32 - t30 * mrSges(5,1) - t31 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t63 = t35 * t16;
t5 = t17 * t37 + t23 * t63;
t62 = t35 * t17;
t6 = t16 * t37 - t23 * t62;
t76 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t59 = t37 * t23;
t7 = -t16 * t59 + t62;
t8 = t17 * t59 + t63;
t75 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t73 = pkin(4) * t30;
t72 = pkin(5) * t20;
t69 = g(3) * t21;
t26 = t36 * pkin(2);
t68 = t21 * mrSges(5,3);
t33 = -pkin(8) - qJ(4);
t27 = -pkin(9) + t33;
t65 = t21 * t27;
t64 = t21 * t33;
t61 = t35 * t20;
t60 = t35 * t22;
t58 = qJ(4) * t21;
t50 = -mrSges(7,1) * t16 - mrSges(7,2) * t17;
t49 = t13 * t23 - t65;
t48 = t18 * t23 - t64;
t11 = -t20 * t59 + t60;
t9 = t22 * t37 + t23 * t61;
t19 = t26 + pkin(1);
t15 = t37 * t19;
t14 = t72 + t73;
t12 = t22 * t59 + t61;
t10 = t20 * t37 - t23 * t60;
t1 = [(-m(5) * t15 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + (-m(4) + t82) * (-t35 * t32 + t15) + (-m(6) * t73 - m(7) * t14 + t80) * t35 + (-m(6) * t48 - m(7) * t49 - t86 * t21 + t84) * t37) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (m(4) * t32 - m(6) * (-t32 + t73) - m(7) * (t14 - t32) + t80) * t37 + (m(4) * t19 - m(5) * (-t19 - t58) + t68 - m(6) * (-t19 - t48) - m(7) * (-t19 - t49) - t84) * t35) * g(1) (-t88 - t92) * t89 + (-m(4) * t26 - m(5) * (t26 + t58) - t68 - m(6) * (t26 - t64) - m(7) * (t26 - t65) + t87 + t88 * t23) * g(3) + (mrSges(3,2) * t36 + (m(6) * t33 + m(7) * t27 - t86 + t90) * t23 + (t91 * pkin(2) + mrSges(3,1)) * t34) * t85 (-g(1) * t35 + g(2) * t37) * t91 (g(3) * t23 - t89) * t57 (m(7) * t72 + mrSges(6,1) * t20 + mrSges(6,2) * t22 - t50) * t69 + (-t10 * mrSges(6,2) + t81 * t9 - t76) * g(2) + (t12 * mrSges(6,2) - t81 * t11 - t75) * g(1), -g(1) * t75 - g(2) * t76 - t50 * t69];
taug  = t1(:);
