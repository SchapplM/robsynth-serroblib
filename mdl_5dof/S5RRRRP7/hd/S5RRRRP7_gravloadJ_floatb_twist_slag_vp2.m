% Calculate Gravitation load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:10
% DurationCPUTime: 0.71s
% Computational Cost: add. (342->90), mult. (455->107), div. (0->0), fcn. (415->8), ass. (0->51)
t104 = mrSges(5,1) + mrSges(6,1);
t103 = -mrSges(6,2) - mrSges(5,3);
t33 = sin(qJ(4));
t36 = cos(qJ(4));
t102 = t33 * mrSges(6,3) + t104 * t36;
t97 = m(5) + m(6);
t38 = cos(qJ(1));
t101 = t103 * t38;
t32 = qJ(2) + qJ(3);
t29 = sin(t32);
t100 = t103 * t29;
t99 = pkin(4) * t36 + qJ(5) * t33;
t93 = (-m(6) * (-pkin(3) - t99) + t102) * t29;
t30 = cos(t32);
t96 = t30 * mrSges(4,1) - t29 * mrSges(4,2);
t62 = t30 * pkin(3) + t29 * pkin(8);
t83 = pkin(3) * t29;
t34 = sin(qJ(2));
t84 = pkin(2) * t34;
t95 = m(6) * t84 - m(5) * (-t83 - t84) + t93;
t73 = t30 * t38;
t19 = pkin(8) * t73;
t76 = t29 * t33;
t60 = mrSges(5,2) * t76;
t94 = -m(6) * t19 + t101 * t30 - t38 * t60;
t35 = sin(qJ(1));
t92 = g(1) * t38 + g(2) * t35;
t37 = cos(qJ(2));
t51 = t37 * mrSges(3,1) - t34 * mrSges(3,2);
t91 = -m(3) * pkin(1) - mrSges(2,1) - t51 - t96;
t90 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t89 = -t96 + t100 + (mrSges(5,2) * t33 - t102) * t30;
t88 = m(6) * pkin(4) + t104;
t87 = (-t60 + (-t97 * pkin(8) + t103) * t30) * t35;
t86 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t31 = t37 * pkin(2);
t75 = t29 * t38;
t71 = t35 * t33;
t70 = t35 * t36;
t65 = t38 * t33;
t64 = t38 * t36;
t28 = t31 + pkin(1);
t39 = -pkin(7) - pkin(6);
t54 = t38 * t28 - t35 * t39;
t52 = t30 * t99 + t62;
t48 = mrSges(4,1) * t29 + mrSges(4,2) * t30;
t4 = t30 * t64 + t71;
t3 = t30 * t65 - t70;
t2 = t30 * t70 - t65;
t1 = t30 * t71 + t64;
t5 = [(-m(4) * t54 - t97 * (pkin(3) * t73 + pkin(8) * t75 + t54) - t88 * t4 - t86 * t3 + t101 * t29 + t91 * t38 + t90 * t35) * g(2) + (t88 * t2 + t86 * t1 + (m(4) * t28 - t97 * (-t28 - t62) - t91 - t100) * t35 + (t90 + (m(4) + t97) * t39) * t38) * g(1), (t95 * t35 + t87) * g(2) + (-m(5) * t19 + t95 * t38 + t94) * g(1) + (-t51 - m(4) * t31 - m(5) * (t31 + t62) - m(6) * (t31 + t52) + t89) * g(3) + (m(4) * t84 + mrSges(3,1) * t34 + mrSges(3,2) * t37 + t48) * t92, t92 * t48 + ((m(5) * t83 + t93) * t35 + t87) * g(2) + (-m(5) * (-pkin(3) * t75 + t19) + t93 * t38 + t94) * g(1) + (-m(5) * t62 - m(6) * t52 + t89) * g(3), (t88 * t33 - t86 * t36) * g(3) * t29 + (t88 * t1 - t86 * t2) * g(2) + (t88 * t3 - t86 * t4) * g(1), (-g(1) * t3 - g(2) * t1 - g(3) * t76) * m(6)];
taug = t5(:);
