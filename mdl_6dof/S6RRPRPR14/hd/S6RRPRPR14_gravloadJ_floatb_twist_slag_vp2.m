% Calculate Gravitation load on the joints for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:35
% EndTime: 2019-03-09 11:32:39
% DurationCPUTime: 1.20s
% Computational Cost: add. (498->112), mult. (1193->144), div. (0->0), fcn. (1368->10), ass. (0->62)
t42 = sin(qJ(6));
t46 = cos(qJ(6));
t103 = m(6) + m(7);
t60 = -qJ(5) * t103 + mrSges(5,2) - mrSges(6,3);
t114 = -t42 * mrSges(7,1) - t46 * mrSges(7,2) + t60;
t62 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t80 = m(5) + t103;
t77 = m(4) + t80;
t112 = t77 * qJ(3);
t102 = mrSges(3,2) - mrSges(4,3);
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t111 = -t114 * t47 + t62 * t43 + t102;
t105 = pkin(4) * t103 - t62;
t99 = -t102 + t112;
t98 = mrSges(4,2) - mrSges(3,1) - mrSges(6,1) - mrSges(5,3);
t64 = -t46 * mrSges(7,1) + t42 * mrSges(7,2);
t97 = -m(7) * (-pkin(5) - pkin(9)) - t64 - t98;
t94 = t80 * pkin(9) - t98;
t93 = t111 - t112;
t92 = m(7) * pkin(5);
t90 = pkin(4) * t43;
t41 = sin(pkin(6));
t44 = sin(qJ(2));
t89 = t41 * t44;
t45 = sin(qJ(1));
t88 = t41 * t45;
t48 = cos(qJ(2));
t87 = t41 * t48;
t49 = cos(qJ(1));
t86 = t41 * t49;
t83 = pkin(2) * t87 + qJ(3) * t89;
t82 = t49 * pkin(1) + pkin(8) * t88;
t81 = cos(pkin(6));
t72 = t45 * t81;
t28 = -t44 * t72 + t48 * t49;
t79 = t28 * pkin(2) + t82;
t78 = pkin(9) * t87 + t83;
t75 = -pkin(1) * t45 + pkin(8) * t86;
t71 = t49 * t81;
t25 = t44 * t45 - t48 * t71;
t19 = t25 * pkin(2);
t74 = -pkin(9) * t25 - t19;
t27 = t49 * t44 + t48 * t72;
t21 = t27 * pkin(2);
t73 = -pkin(9) * t27 - t21;
t70 = pkin(3) * t88 + t79;
t26 = t44 * t71 + t45 * t48;
t68 = -t26 * pkin(2) + t75;
t66 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t41;
t61 = pkin(3) * t86 + t68;
t59 = -t25 * t43 + t47 * t86;
t9 = t25 * t47 + t43 * t86;
t56 = -t64 + t92;
t23 = t81 * t43 + t47 * t87;
t14 = t28 * t90;
t13 = t26 * t90;
t8 = t27 * t43 + t47 * t88;
t7 = -t27 * t47 + t43 * t88;
t2 = t28 * t46 + t42 * t7;
t1 = -t28 * t42 + t46 * t7;
t3 = [(-t49 * mrSges(2,1) - m(3) * t82 - m(4) * t79 - m(5) * t70 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t60 * t7 - t99 * t27 + t62 * t8 + t66 * t45 + (-t92 - t94) * t28 - t103 * (t8 * pkin(4) + t70)) * g(2) + (t45 * mrSges(2,1) - m(3) * t75 - m(4) * t68 - m(5) * t61 + t99 * t25 + t66 * t49 + t62 * t59 + t114 * t9 + (t56 + t94) * t26 + t103 * (-pkin(4) * t59 - t61)) * g(1) (-m(4) * t83 - m(5) * t78 - t103 * (t89 * t90 + t78) + ((-t56 + t98) * t48 + t111 * t44) * t41) * g(3) + (m(4) * t19 - m(5) * t74 - m(6) * (t13 + t74) - m(7) * (t13 - t19) + t93 * t26 + t97 * t25) * g(2) + (m(4) * t21 - m(5) * t73 - m(6) * (t14 + t73) - m(7) * (t14 - t21) + t93 * t28 + t97 * t27) * g(1) (-g(1) * t27 - g(2) * t25 + g(3) * t87) * t77 (t114 * (-t43 * t87 + t81 * t47) + t105 * t23) * g(3) + (-t105 * t9 - t114 * t59) * g(2) + (t105 * t7 + t114 * t8) * g(1), t103 * (-g(1) * t7 + g(2) * t9 - g(3) * t23) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t26 * t42 - t46 * t9) * mrSges(7,1) + (-t26 * t46 + t42 * t9) * mrSges(7,2)) - g(3) * ((t23 * t46 - t42 * t89) * mrSges(7,1) + (-t23 * t42 - t46 * t89) * mrSges(7,2))];
taug  = t3(:);
