% Calculate Gravitation load on the joints for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:06
% EndTime: 2019-12-05 18:55:09
% DurationCPUTime: 0.64s
% Computational Cost: add. (334->95), mult. (409->113), div. (0->0), fcn. (365->10), ass. (0->67)
t113 = mrSges(5,3) + mrSges(6,3);
t39 = qJ(4) + qJ(5);
t34 = sin(t39);
t41 = sin(qJ(4));
t114 = -t41 * mrSges(5,2) - t34 * mrSges(6,2);
t112 = m(5) + m(6);
t111 = pkin(5) * t112;
t40 = qJ(2) + qJ(3);
t35 = sin(t40);
t110 = t114 * t35;
t46 = cos(qJ(1));
t109 = t113 * t46;
t36 = cos(t39);
t44 = cos(qJ(4));
t108 = -t44 * mrSges(5,1) - t36 * mrSges(6,1);
t33 = t44 * pkin(3) + pkin(2);
t107 = m(5) * pkin(2) + m(6) * t33 - t108;
t106 = t113 * t35;
t31 = t35 * pkin(5);
t45 = cos(qJ(2));
t38 = t45 * pkin(1);
t56 = -t38 - t31;
t42 = sin(qJ(2));
t84 = t35 * mrSges(4,2);
t105 = -m(4) * t38 - t45 * mrSges(3,1) + t42 * mrSges(3,2) + t84;
t37 = cos(t40);
t43 = sin(qJ(1));
t104 = (t110 + (-t111 - t113) * t37) * t43;
t81 = t37 * t46;
t103 = -t109 * t37 + t110 * t46 - t81 * t111;
t94 = m(6) * pkin(3);
t102 = mrSges(5,1) + t94;
t90 = pkin(1) * t42;
t101 = t107 * t35 + t112 * t90;
t30 = t37 * mrSges(4,1);
t97 = -mrSges(2,1) - t30 + t105;
t96 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t95 = -t30 - t106 + (t108 - t114) * t37;
t71 = t46 * t36;
t79 = t43 * t34;
t5 = t37 * t79 + t71;
t72 = t46 * t34;
t78 = t43 * t36;
t6 = -t37 * t78 + t72;
t92 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t37 * t72 + t78;
t8 = t37 * t71 + t79;
t91 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t87 = g(3) * t35;
t32 = t37 * pkin(2);
t86 = mrSges(4,2) * t37;
t17 = t37 * t33;
t77 = t43 * t41;
t76 = t43 * t44;
t70 = t46 * t41;
t69 = t46 * t44;
t64 = t17 + t31;
t61 = t56 * t46;
t60 = t32 + t31;
t59 = t41 * t94;
t53 = -mrSges(6,1) * t34 - mrSges(6,2) * t36;
t11 = -t37 * t70 + t76;
t9 = t37 * t77 + t69;
t47 = t86 + (mrSges(4,1) + t107) * t35;
t12 = t37 * t69 + t77;
t10 = -t37 * t76 + t70;
t1 = [(-m(5) * (pkin(2) * t81 - t61) - t12 * mrSges(5,1) - t11 * mrSges(5,2) - m(6) * (pkin(3) * t77 + t33 * t81 - t61) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - t109 * t35 + t97 * t46 + t96 * t43) * g(2) + (-t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (-t59 + t96) * t46 + (-m(5) * (t56 - t32) - m(6) * (t56 - t17) - t97 + t106) * t43) * g(1), (t101 * t43 + t104) * g(2) + (t101 * t46 + t103) * g(1) + (-m(5) * (t38 + t60) - m(6) * (t38 + t64) + t95 + t105) * g(3) + (m(4) * t90 + mrSges(3,1) * t42 + mrSges(4,1) * t35 + mrSges(3,2) * t45 + t86) * (g(1) * t46 + g(2) * t43), (-m(5) * t60 - m(6) * t64 + t84 + t95) * g(3) + (t47 * t43 + t104) * g(2) + (t47 * t46 + t103) * g(1), (mrSges(5,1) * t41 + mrSges(5,2) * t44 - t53 + t59) * t87 + (-t10 * mrSges(5,2) + t102 * t9 - t92) * g(2) + (t12 * mrSges(5,2) - t102 * t11 - t91) * g(1), -g(1) * t91 - g(2) * t92 - t53 * t87];
taug = t1(:);
