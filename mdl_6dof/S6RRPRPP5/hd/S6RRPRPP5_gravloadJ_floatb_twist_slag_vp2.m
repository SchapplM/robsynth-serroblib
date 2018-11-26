% Calculate Gravitation load on the joints for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2018-11-23 16:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:59:02
% EndTime: 2018-11-23 16:59:03
% DurationCPUTime: 0.88s
% Computational Cost: add. (280->95), mult. (621->114), div. (0->0), fcn. (584->6), ass. (0->48)
t81 = mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t75 = m(7) * pkin(5) + t81;
t79 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t85 = m(6) + m(7);
t90 = t85 * qJ(5) - t79;
t78 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t88 = -t75 * t27 + t90 * t30;
t87 = -m(7) * qJ(6) - t78 + (-m(5) - t85) * (-pkin(2) - pkin(8));
t86 = -m(4) - m(5);
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t84 = (mrSges(3,1) - mrSges(4,2)) * t31 + (-mrSges(3,2) + mrSges(4,3)) * t28;
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t83 = g(1) * t32 + g(2) * t29;
t82 = mrSges(2,1) + t84;
t80 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t77 = t85 - t86;
t76 = (-mrSges(4,3) + t88) * t31 + (m(4) * pkin(2) - mrSges(4,2) + t87) * t28;
t73 = pkin(4) * t27;
t70 = g(3) * t31;
t23 = t31 * pkin(2);
t68 = t27 * t29;
t67 = t27 * t32;
t66 = t29 * t30;
t65 = t30 * t32;
t64 = t31 * t32;
t19 = t28 * qJ(3);
t62 = t23 + t19;
t24 = t32 * pkin(7);
t61 = t32 * pkin(3) + t24;
t60 = t32 * pkin(1) + t29 * pkin(7);
t59 = qJ(3) * t31;
t58 = qJ(6) * t31;
t57 = t31 * t73;
t56 = t31 * pkin(8) + t62;
t51 = t28 * t73 + t56;
t50 = pkin(2) * t64 + t32 * t19 + t60;
t41 = t29 * pkin(3) + pkin(8) * t64 + t50;
t14 = t32 * t59;
t13 = t29 * t59;
t8 = -t28 * t68 + t65;
t7 = t28 * t66 + t67;
t6 = t28 * t67 + t66;
t5 = -t28 * t65 + t68;
t1 = [(-m(3) * t60 - m(4) * t50 - m(5) * t41 - t85 * (t6 * pkin(4) + t5 * qJ(5) + t41) + t78 * t64 - t75 * t6 + t79 * t5 + (m(7) * t58 - t82) * t32 + t80 * t29) * g(2) + (-m(5) * t61 - t85 * (t8 * pkin(4) + qJ(5) * t7 + t61) + (-m(3) - m(4)) * t24 - t75 * t8 + t79 * t7 + t80 * t32 + (m(3) * pkin(1) + m(4) * t23 - t77 * (-pkin(1) - t19) + t87 * t31 + t82) * t29) * g(1), t83 * (mrSges(3,1) * t28 + mrSges(3,2) * t31) + (-t85 * (t29 * t57 + t13) + t86 * t13 + t76 * t29) * g(2) + (-t85 * (t32 * t57 + t14) + t86 * t14 + t76 * t32) * g(1) + (-m(4) * t62 - m(5) * t56 - m(6) * t51 - m(7) * (t51 - t58) + t78 * t31 + t88 * t28 - t84) * g(3) (-t83 * t28 + t70) * t77 ((m(6) * pkin(4) - m(7) * (-pkin(4) - pkin(5)) + t81) * t30 + t90 * t27) * t70 + (-t85 * (t7 * pkin(4) - qJ(5) * t8) - t79 * t8 - t75 * t7) * g(2) + (-t85 * (-t5 * pkin(4) + qJ(5) * t6) + t79 * t6 + t75 * t5) * g(1), t85 * (-g(1) * t5 + g(2) * t7 - t30 * t70) (g(3) * t28 + t83 * t31) * m(7)];
taug  = t1(:);
