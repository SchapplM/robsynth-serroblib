% Calculate Gravitation load on the joints for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2018-11-23 16:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:57:24
% EndTime: 2018-11-23 16:57:25
% DurationCPUTime: 0.87s
% Computational Cost: add. (430->107), mult. (582->119), div. (0->0), fcn. (549->8), ass. (0->56)
t91 = mrSges(5,1) + mrSges(6,1);
t80 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t85 = m(6) + m(7);
t92 = m(5) + t85;
t90 = mrSges(6,2) + mrSges(5,3);
t89 = mrSges(7,3) - t90;
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t88 = -t80 * t27 + t91 * t30;
t29 = sin(qJ(1));
t87 = g(2) * t29;
t25 = qJ(2) + pkin(9);
t22 = sin(t25);
t70 = g(3) * t22;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t47 = t31 * mrSges(3,1) - t28 * mrSges(3,2);
t86 = m(3) * pkin(1) + mrSges(2,1) + t47;
t32 = cos(qJ(1));
t83 = g(1) * t32 + t87;
t81 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t23 = cos(t25);
t49 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t59 = qJ(5) * t27;
t51 = -pkin(3) - t59;
t73 = pkin(4) * t30;
t74 = pkin(2) * t28;
t79 = t92 * t74 + (m(7) * qJ(6) + t89) * t23 + (-m(7) * t51 - t49 * t30 - m(6) * (t51 - t73) + m(5) * pkin(3) + t88) * t22;
t45 = t23 * mrSges(4,1) - t22 * mrSges(4,2);
t78 = t90 * t22 + t45;
t77 = m(7) * pkin(5) + mrSges(7,1);
t76 = pkin(8) * t92;
t75 = t77 + t91;
t17 = t22 * pkin(8);
t18 = t23 * pkin(3);
t24 = t31 * pkin(2);
t67 = t22 * t32;
t66 = t23 * t32;
t26 = -qJ(3) - pkin(7);
t65 = t26 * t32;
t63 = t29 * t27;
t62 = t29 * t30;
t61 = t30 * t32;
t60 = t32 * t27;
t58 = qJ(6) * t22;
t57 = t18 + t17 + t24;
t21 = t24 + pkin(1);
t54 = -t21 - t18;
t50 = t32 * t21 - t29 * t26;
t48 = t57 + (t59 + t73) * t23;
t42 = pkin(3) * t66 + pkin(8) * t67 + t50;
t8 = t23 * t61 + t63;
t7 = t23 * t60 - t62;
t6 = t23 * t62 - t60;
t5 = t23 * t63 + t61;
t1 = [(-m(4) * t50 - m(5) * t42 - t85 * (t8 * pkin(4) + t7 * qJ(5) + t42) - t75 * t8 + t80 * t7 + t89 * t67 + (m(7) * t58 - t45 - t86) * t32 + t81 * t29) * g(2) + (m(5) * t65 - t85 * (-t6 * pkin(4) - t5 * qJ(5) - t65) + t75 * t6 - t80 * t5 + (m(4) * t26 + t81) * t32 + (m(4) * t21 - m(7) * t54 - (m(7) * (-pkin(8) + qJ(6)) + mrSges(7,3)) * t22 + (-m(5) - m(6)) * (t54 - t17) + t78 + t86) * t29) * g(1) (t79 * t32 - t66 * t76) * g(1) + (-t47 - m(4) * t24 - m(5) * t57 - m(6) * t48 - m(7) * (t48 - t58) + t22 * mrSges(7,3) + (-t77 * t30 - t88) * t23 - t78) * g(3) + (m(4) * t74 + mrSges(3,1) * t28 + mrSges(4,1) * t22 + mrSges(3,2) * t31 + mrSges(4,2) * t23) * t83 + (-t23 * t76 + t79) * t87 (-g(1) * t29 + g(2) * t32) * (m(4) + t92) (-t85 * (-t5 * pkin(4) + qJ(5) * t6) + t80 * t6 + t75 * t5) * g(2) + (-t85 * (-t7 * pkin(4) + qJ(5) * t8) + t80 * t8 + t75 * t7) * g(1) + ((-t85 * qJ(5) + t80) * t30 + (m(6) * pkin(4) - t49 + t91) * t27) * t70, t85 * (-g(1) * t7 - g(2) * t5 - t27 * t70) (-g(3) * t23 + t83 * t22) * m(7)];
taug  = t1(:);
