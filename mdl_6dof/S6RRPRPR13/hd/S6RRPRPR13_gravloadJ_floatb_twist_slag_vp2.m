% Calculate Gravitation load on the joints for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:12
% EndTime: 2018-11-23 17:09:13
% DurationCPUTime: 0.94s
% Computational Cost: add. (1254->106), mult. (1471->132), div. (0->0), fcn. (1429->16), ass. (0->60)
t42 = pkin(11) + qJ(6);
t39 = sin(t42);
t40 = cos(t42);
t43 = sin(pkin(11));
t45 = cos(pkin(11));
t56 = -m(7) * (pkin(5) * t45 + pkin(4)) - m(6) * pkin(4) - mrSges(6,1) * t45 + mrSges(6,2) * t43 - mrSges(5,1);
t115 = -t40 * mrSges(7,1) + t39 * mrSges(7,2) + t56;
t111 = m(6) + m(7);
t88 = m(5) + t111;
t110 = m(4) + t88;
t104 = qJ(3) * t110 - mrSges(3,2) + mrSges(4,3);
t113 = -t45 * mrSges(6,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t69 = t39 * mrSges(7,1) + t40 * mrSges(7,2);
t112 = pkin(2) * t110 + t43 * mrSges(6,1) - t113 + (m(5) + m(6)) * pkin(9) + m(7) * (pkin(5) * t43 + pkin(9)) + t69;
t114 = m(6) * qJ(5) - m(7) * (-pkin(10) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t102 = -t88 * pkin(9) - (m(7) * pkin(5) + mrSges(6,1)) * t43 + t113;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t99 = t114 * t51 + t115 * t48 - t104;
t44 = sin(pkin(6));
t50 = sin(qJ(1));
t95 = t44 * t50;
t53 = cos(qJ(1));
t94 = t44 * t53;
t52 = cos(qJ(2));
t92 = t50 * t52;
t91 = t53 * t52;
t89 = t53 * pkin(1) + pkin(8) * t95;
t87 = pkin(6) - qJ(2);
t86 = pkin(6) + qJ(2);
t75 = sin(t86);
t71 = t75 / 0.2e1;
t76 = sin(t87);
t65 = t71 - t76 / 0.2e1;
t22 = -t50 * t65 + t91;
t85 = t22 * pkin(2) + t89;
t80 = -pkin(1) * t50 + pkin(8) * t94;
t77 = cos(t86);
t19 = t53 * t65 + t92;
t74 = -t19 * pkin(2) + t80;
t73 = cos(t87) / 0.2e1;
t72 = t76 / 0.2e1;
t70 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t44;
t66 = t72 - t75 / 0.2e1;
t49 = sin(qJ(2));
t58 = t73 + t77 / 0.2e1;
t18 = t49 * t50 - t53 * t58;
t7 = -t18 * t48 + t51 * t94;
t5 = t18 * t51 + t48 * t94;
t46 = cos(pkin(6));
t30 = t73 - t77 / 0.2e1;
t29 = t71 + t72;
t21 = t53 * t49 + t50 * t58;
t17 = -t29 * t48 + t46 * t51;
t16 = t29 * t51 + t46 * t48;
t4 = t21 * t48 + t51 * t95;
t3 = -t21 * t51 + t48 * t95;
t2 = t22 * t39 + t4 * t40;
t1 = t22 * t40 - t39 * t4;
t6 = [(-m(3) * t89 - m(4) * t85 - t53 * mrSges(2,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t102 * t22 - t104 * t21 - t114 * t3 + t56 * t4 + t70 * t50 - t88 * (pkin(3) * t95 + t85)) * g(2) + (t50 * mrSges(2,1) - m(3) * t80 - m(4) * t74 + t104 * t18 - t114 * t5 + t70 * t53 + t115 * t7 + (-t102 + t69) * t19 + t88 * (-pkin(3) * t94 - t74)) * g(1) (-t112 * t29 + t99 * t30) * g(3) + (t99 * (-t53 * t66 + t92) + t112 * t18) * g(2) + (t99 * (t50 * t66 + t91) + t112 * t21) * g(1) (-g(1) * t21 - g(2) * t18 + g(3) * t29) * t110 (-t114 * t17 - t115 * t16) * g(3) + (t114 * t7 + t115 * t5) * g(2) + (-t114 * t4 - t115 * t3) * g(1), t111 * (-g(1) * t3 + g(2) * t5 - g(3) * t16) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t40 + t39 * t7) * mrSges(7,1) + (-t19 * t39 + t40 * t7) * mrSges(7,2)) - g(3) * ((-t17 * t39 + t30 * t40) * mrSges(7,1) + (-t17 * t40 - t30 * t39) * mrSges(7,2))];
taug  = t6(:);
