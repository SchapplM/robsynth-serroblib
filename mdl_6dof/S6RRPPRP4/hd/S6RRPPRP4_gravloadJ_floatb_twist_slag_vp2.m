% Calculate Gravitation load on the joints for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:46:46
% EndTime: 2018-11-23 16:46:47
% DurationCPUTime: 0.96s
% Computational Cost: add. (365->98), mult. (877->121), div. (0->0), fcn. (928->8), ass. (0->50)
t107 = mrSges(4,1) + mrSges(5,1);
t105 = mrSges(4,2) - mrSges(5,3);
t114 = m(6) + m(7);
t43 = sin(pkin(9));
t44 = cos(pkin(9));
t113 = -t105 * t43 + t107 * t44;
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t104 = g(1) * t50 + g(2) * t47;
t112 = mrSges(4,3) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t46 = sin(qJ(2));
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t64 = t43 * t45 + t44 * t48;
t60 = t46 * t64;
t106 = mrSges(2,2) - mrSges(3,3);
t38 = t46 * qJ(3);
t49 = cos(qJ(2));
t80 = t49 * pkin(2) + t38;
t103 = m(5) + t114;
t102 = m(4) + t103;
t91 = t43 * t48;
t101 = t46 * (t44 * t45 - t91);
t72 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t99 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t69 = t49 * mrSges(3,1) - t46 * mrSges(3,2);
t97 = t112 * t46 + t69;
t94 = pkin(8) * t46;
t89 = t44 * t49;
t84 = t46 * t50;
t83 = t47 * t49;
t82 = t49 * t50;
t81 = t50 * t43;
t79 = t50 * pkin(1) + t47 * pkin(7);
t77 = qJ(4) * t43;
t75 = -pkin(2) - t77;
t74 = pkin(3) * t89 + t49 * t77 + t80;
t73 = pkin(2) * t82 + t50 * t38 + t79;
t22 = t43 * t83 + t44 * t50;
t23 = t44 * t83 - t81;
t41 = t50 * pkin(7);
t70 = -t23 * pkin(3) - qJ(4) * t22 + t41;
t3 = t22 * t48 - t23 * t45;
t65 = t22 * t45 + t23 * t48;
t24 = -t47 * t44 + t49 * t81;
t25 = t47 * t43 + t44 * t82;
t58 = t25 * pkin(3) + t24 * qJ(4) + t73;
t6 = t24 * t45 + t25 * t48;
t5 = -t24 * t48 + t25 * t45;
t1 = [(-m(3) * t79 - m(4) * t73 - m(5) * t58 - t114 * (t25 * pkin(4) - pkin(8) * t84 + t58) - t72 * t6 + (-mrSges(2,1) - t69) * t50 - t99 * t5 + t106 * t47 - t107 * t25 + t105 * t24 - t112 * t84) * g(2) + (-m(5) * t70 - t114 * (-t23 * pkin(4) + t47 * t94 + t70) + t106 * t50 + (-m(3) - m(4)) * t41 + t72 * t65 - t99 * t3 + t107 * t23 - t105 * t22 + (m(3) * pkin(1) + mrSges(2,1) - t102 * (-pkin(1) - t80) + t97) * t47) * g(1) (-m(4) * t80 - m(5) * t74 - t114 * (pkin(4) * t89 + t74 - t94) - t99 * t45 * t89 - t97 + (-t72 * t64 + t99 * t91 - t113) * t49) * g(3) + (t101 * t99 + t60 * t72 + (t114 * pkin(8) - qJ(3) * t102 + mrSges(3,2) - t112) * t49 + (mrSges(3,1) - t114 * ((-pkin(3) - pkin(4)) * t44 + t75) + m(4) * pkin(2) - m(5) * (-pkin(3) * t44 + t75) + t113) * t46) * t104 (t49 * g(3) - t104 * t46) * t102, t103 * (-g(3) * t43 * t46 - g(1) * t24 - g(2) * t22) (t101 * t72 - t60 * t99) * g(3) + (-t3 * t72 - t65 * t99) * g(2) + (t5 * t72 - t6 * t99) * g(1) (-g(1) * t5 + g(2) * t3 - g(3) * t101) * m(7)];
taug  = t1(:);
