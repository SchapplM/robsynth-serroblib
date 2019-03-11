% Calculate Gravitation load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:13
% EndTime: 2019-03-09 09:01:15
% DurationCPUTime: 1.01s
% Computational Cost: add. (679->110), mult. (1653->156), div. (0->0), fcn. (2020->12), ass. (0->60)
t98 = m(6) + m(7);
t112 = m(5) + t98;
t102 = qJ(4) * t112 - mrSges(4,2) + mrSges(5,3);
t47 = sin(qJ(6));
t51 = cos(qJ(6));
t105 = m(7) * pkin(5) + t51 * mrSges(7,1) - t47 * mrSges(7,2) + mrSges(6,1);
t103 = -t98 * pkin(9) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t100 = t47 * mrSges(7,1) + t51 * mrSges(7,2) - t103;
t70 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t46 = cos(pkin(6));
t50 = sin(qJ(1));
t53 = cos(qJ(2));
t90 = t50 * t53;
t49 = sin(qJ(2));
t54 = cos(qJ(1));
t93 = t49 * t54;
t33 = -t46 * t90 - t93;
t63 = t33 * pkin(2);
t88 = t53 * t54;
t91 = t50 * t49;
t108 = t46 * t88 - t91;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t99 = t105 * t48 + t70 * t52 + t102;
t82 = sin(pkin(11));
t83 = cos(pkin(11));
t36 = -t49 * t83 - t53 * t82;
t106 = t46 * t36;
t61 = -t49 * t82 + t53 * t83;
t21 = t106 * t50 + t61 * t54;
t107 = t106 * t54 - t50 * t61;
t97 = pkin(2) * t53;
t45 = sin(pkin(6));
t95 = t45 * t50;
t94 = t45 * t54;
t44 = pkin(1) + t97;
t39 = t54 * t44;
t85 = t21 * pkin(3) + t39;
t41 = t45 * t97;
t77 = pkin(4) * t95 + t85;
t76 = m(4) + t112;
t75 = -m(3) * pkin(1) - mrSges(2,1);
t69 = t108 * pkin(2);
t58 = t46 * t61;
t17 = t50 * t36 + t54 * t58;
t8 = t17 * t48 + t52 * t94;
t6 = -t17 * t52 + t48 * t94;
t56 = mrSges(2,2) + (-m(3) * pkin(8) - mrSges(5,1) - mrSges(3,3) - mrSges(4,3)) * t45 + t76 * (pkin(2) * t46 * t49 + (-pkin(8) - qJ(3)) * t45);
t34 = -t46 * t91 + t88;
t32 = -t46 * t93 - t90;
t29 = t36 * t45;
t28 = t61 * t45;
t23 = -t28 * t48 + t46 * t52;
t20 = t36 * t54 - t50 * t58;
t13 = t107 * pkin(3);
t4 = -t20 * t48 + t52 * t95;
t3 = t20 * t52 + t48 * t95;
t2 = t21 * t47 + t4 * t51;
t1 = t21 * t51 - t4 * t47;
t5 = [(-t34 * mrSges(3,1) - t33 * mrSges(3,2) - m(4) * t39 - m(5) * t85 - m(6) * t77 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t77) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t75 * t54 + t70 * t3 + t102 * t20 + t103 * t21 + t56 * t50) * g(2) + (-t32 * mrSges(3,1) + t108 * mrSges(3,2) - m(5) * t13 + t70 * t6 - t102 * t17 - t105 * t8 + (t76 * t44 - t75) * t50 - t100 * t107 + t56 * t54 + t98 * (-pkin(4) * t94 - t13)) * g(1) (-(mrSges(3,1) * t53 - mrSges(3,2) * t49) * t45 - m(4) * t41 - t112 * (t28 * pkin(3) + t41) + t99 * t29 - t100 * t28) * g(3) + (-m(4) * t69 - mrSges(3,1) * t108 - mrSges(3,2) * t32 - t112 * (t17 * pkin(3) + t69) - t100 * t17 + t99 * t107) * g(2) + (-m(4) * t63 - mrSges(3,1) * t33 + mrSges(3,2) * t34 - t100 * t20 - t99 * t21 - t112 * (t20 * pkin(3) + t63)) * g(1) (-t46 * g(3) + (-t50 * g(1) + t54 * g(2)) * t45) * t76, t112 * (g(1) * t20 + g(2) * t17 + g(3) * t28) (t70 * t23 - t105 * (-t28 * t52 - t46 * t48)) * g(3) + (-t105 * t6 - t70 * t8) * g(2) + (t105 * t3 + t4 * t70) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t107 * t51 + t47 * t8) * mrSges(7,1) + (t107 * t47 + t51 * t8) * mrSges(7,2)) - g(3) * ((-t23 * t47 - t29 * t51) * mrSges(7,1) + (-t23 * t51 + t29 * t47) * mrSges(7,2))];
taug  = t5(:);
