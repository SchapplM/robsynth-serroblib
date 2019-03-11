% Calculate Gravitation load on the joints for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:21
% EndTime: 2019-03-09 14:02:23
% DurationCPUTime: 0.81s
% Computational Cost: add. (592->130), mult. (597->149), div. (0->0), fcn. (549->12), ass. (0->68)
t105 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t46 = pkin(11) + qJ(4);
t37 = cos(t46);
t31 = pkin(4) * t37;
t48 = cos(pkin(11));
t34 = t48 * pkin(3) + pkin(2);
t25 = t31 + t34;
t38 = qJ(5) + t46;
t33 = cos(t38);
t27 = pkin(5) * t33;
t17 = t27 + t25;
t35 = qJ(6) + t38;
t28 = sin(t35);
t29 = cos(t35);
t32 = sin(t38);
t36 = sin(t46);
t104 = -m(5) * t34 - m(6) * t25 - m(7) * t17 - t37 * mrSges(5,1) - t33 * mrSges(6,1) - t29 * mrSges(7,1) + t36 * mrSges(5,2) + t32 * mrSges(6,2) + t28 * mrSges(7,2);
t49 = -pkin(8) - qJ(3);
t45 = -pkin(9) + t49;
t39 = -pkin(10) + t45;
t103 = m(5) * t49 + m(6) * t45 + m(7) * t39 - t105;
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t102 = g(1) * t53 + g(2) * t51;
t52 = cos(qJ(2));
t47 = sin(pkin(11));
t62 = m(4) * pkin(2) + t48 * mrSges(4,1) - t47 * mrSges(4,2);
t101 = t62 * t52;
t100 = m(6) * pkin(4) + mrSges(5,1);
t66 = -mrSges(7,1) * t28 - mrSges(7,2) * t29;
t99 = mrSges(6,1) * t32 + mrSges(6,2) * t33 - t66;
t78 = t52 * t53;
t15 = -t32 * t78 + t51 * t33;
t16 = t51 * t32 + t33 * t78;
t7 = -t28 * t78 + t51 * t29;
t8 = t51 * t28 + t29 * t78;
t88 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t98 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t88;
t79 = t51 * t52;
t13 = t32 * t79 + t33 * t53;
t14 = t32 * t53 - t33 * t79;
t5 = t28 * t79 + t29 * t53;
t6 = t28 * t53 - t29 * t79;
t89 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t97 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t89;
t96 = m(4) + m(5) + m(6) + m(7);
t95 = -m(3) - t96;
t50 = sin(qJ(2));
t70 = t52 * mrSges(3,1) - t50 * mrSges(3,2);
t94 = t105 * t50 + mrSges(2,1) + t70;
t30 = pkin(4) * t36;
t87 = pkin(5) * t32;
t23 = -t30 - t87;
t40 = t47 * pkin(3);
t92 = -m(5) * t40 - m(6) * (t30 + t40) - m(7) * (-t23 + t40) + mrSges(2,2) - mrSges(3,3) - mrSges(4,1) * t47 - t48 * mrSges(4,2);
t90 = m(7) * pkin(5);
t84 = g(3) * t50;
t71 = m(4) * qJ(3) + mrSges(4,3);
t65 = t17 * t52 - t50 * t39;
t64 = t25 * t52 - t45 * t50;
t63 = t34 * t52 - t49 * t50;
t20 = -t36 * t78 + t51 * t37;
t18 = t36 * t79 + t37 * t53;
t58 = t50 * t71 + t101;
t24 = t27 + t31;
t21 = t51 * t36 + t37 * t78;
t19 = t36 * t53 - t37 * t79;
t1 = [(-t21 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t20 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) + t95 * (t53 * pkin(1) + t51 * pkin(7)) + t92 * t51 + (-m(5) * t63 - m(6) * t64 - m(7) * t65 - t58 - t94) * t53) * g(2) + (-t19 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t18 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (m(3) * pkin(1) - m(4) * (-qJ(3) * t50 - pkin(1)) + t50 * mrSges(4,3) + t101 - m(5) * (-pkin(1) - t63) - m(6) * (-pkin(1) - t64) - m(7) * (-pkin(1) - t65) + t94) * t51 + (t95 * pkin(7) + t92) * t53) * g(1) (-t58 - t70) * g(3) + (t104 * g(3) + t102 * (mrSges(3,2) - t71 + t103)) * t52 + (t103 * g(3) + t102 * (mrSges(3,1) + t62 - t104)) * t50 (g(3) * t52 - t102 * t50) * t96 (m(6) * t30 - m(7) * t23 + mrSges(5,1) * t36 + mrSges(5,2) * t37 + t99) * t84 + (-t19 * mrSges(5,2) - m(7) * (t23 * t79 - t24 * t53) + t100 * t18 + t97) * g(2) + (t21 * mrSges(5,2) - m(7) * (t23 * t78 + t51 * t24) - t100 * t20 + t98) * g(1) (m(7) * t87 + t99) * t84 + (t13 * t90 + t97) * g(2) + (-t15 * t90 + t98) * g(1), -g(1) * t88 - g(2) * t89 - t66 * t84];
taug  = t1(:);
