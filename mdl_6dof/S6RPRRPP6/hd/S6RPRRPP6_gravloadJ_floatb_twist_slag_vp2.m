% Calculate Gravitation load on the joints for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:27
% EndTime: 2019-03-09 04:46:29
% DurationCPUTime: 0.80s
% Computational Cost: add. (315->99), mult. (494->123), div. (0->0), fcn. (451->8), ass. (0->55)
t28 = -qJ(5) - pkin(8);
t81 = m(6) + m(7);
t90 = t81 * t28;
t89 = mrSges(6,1) + mrSges(7,1);
t88 = -mrSges(6,2) + mrSges(7,3);
t87 = -m(4) - m(5);
t86 = mrSges(6,3) + mrSges(7,2);
t27 = qJ(4) + pkin(9);
t21 = sin(t27);
t22 = cos(t27);
t85 = t88 * t21 + t89 * t22;
t30 = sin(qJ(3));
t32 = cos(qJ(4));
t34 = cos(qJ(1));
t60 = t32 * t34;
t29 = sin(qJ(4));
t31 = sin(qJ(1));
t63 = t31 * t29;
t5 = -t30 * t63 + t60;
t62 = t31 * t32;
t66 = t30 * t34;
t7 = t29 * t66 + t62;
t77 = -g(1) * t31 + g(2) * t34;
t84 = -m(3) + t87;
t83 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t33 = cos(qJ(3));
t45 = t30 * mrSges(4,1) + t33 * mrSges(4,2);
t82 = mrSges(2,2) - mrSges(3,3) - t45 - m(5) * (pkin(3) * t30 - pkin(8) * t33) + t33 * mrSges(5,3);
t76 = m(7) * pkin(5) + t89;
t75 = m(7) * qJ(6) + t88;
t20 = pkin(4) * t32 + pkin(3);
t42 = pkin(5) * t22 + qJ(6) * t21;
t80 = -m(7) * (-t20 - t42) + m(6) * t20 + t85;
t78 = -t86 + t90;
t72 = -pkin(1) - pkin(7);
t71 = pkin(4) * t29;
t68 = g(3) * t33;
t67 = t29 * t34;
t65 = t31 * t21;
t64 = t31 * t22;
t61 = t31 * t33;
t59 = t33 * t34;
t58 = t34 * t22;
t56 = t34 * pkin(1) + t31 * qJ(2);
t52 = t34 * pkin(7) + t56;
t51 = m(5) * pkin(8) + mrSges(5,3);
t38 = m(5) * pkin(3) + t32 * mrSges(5,1) - t29 * mrSges(5,2);
t24 = t34 * qJ(2);
t8 = t30 * t60 - t63;
t6 = t30 * t62 + t67;
t4 = t30 * t58 - t65;
t3 = t21 * t66 + t64;
t2 = t21 * t34 + t30 * t64;
t1 = t30 * t65 - t58;
t9 = [(-m(3) * t56 - t6 * mrSges(5,1) - t5 * mrSges(5,2) + t86 * t61 + t87 * t52 - t81 * (t31 * t30 * t20 + pkin(4) * t67 + t28 * t61 + t52) - t76 * t2 - t75 * t1 + t83 * t34 + t82 * t31) * g(2) + (-t8 * mrSges(5,1) + t7 * mrSges(5,2) - t81 * ((-t71 + t72) * t31 + t20 * t66 + t28 * t59 + t24) + t86 * t59 - t76 * t4 - t75 * t3 + t84 * t24 + (m(3) * pkin(1) + t87 * t72 - t83) * t31 + t82 * t34) * g(1), t77 * (t81 - t84) ((t86 * t30 + t80 * t33) * t34 - t66 * t90) * g(2) + (((-m(7) * t42 - t85) * t33 + t78 * t30) * t31 - t81 * t20 * t61) * g(1) + (t45 + (-t51 + t78) * t33 + (t38 + t80) * t30) * g(3) + ((mrSges(4,1) + t38) * t33 + (-mrSges(4,2) + t51) * t30) * t77 (mrSges(5,1) * t29 + mrSges(5,2) * t32 + t76 * t21 - t75 * t22 + t81 * t71) * t68 + (-t7 * mrSges(5,1) - t8 * mrSges(5,2) - t76 * t3 + t75 * t4) * g(2) + (-t5 * mrSges(5,1) + t6 * mrSges(5,2) + t76 * t1 - t75 * t2) * g(1) + (-t5 * g(1) - t7 * g(2)) * t81 * pkin(4) (-t30 * g(3) - t33 * t77) * t81 (-g(1) * t1 + g(2) * t3 - t21 * t68) * m(7)];
taug  = t9(:);
