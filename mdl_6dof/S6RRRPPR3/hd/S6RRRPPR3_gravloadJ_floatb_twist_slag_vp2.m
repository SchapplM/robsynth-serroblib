% Calculate Gravitation load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2018-11-23 17:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:34:02
% EndTime: 2018-11-23 17:34:03
% DurationCPUTime: 0.77s
% Computational Cost: add. (414->125), mult. (493->131), div. (0->0), fcn. (415->8), ass. (0->63)
t112 = -mrSges(6,1) - mrSges(5,3);
t107 = m(6) + m(7);
t111 = -m(7) * pkin(5) + t112;
t36 = qJ(2) + qJ(3);
t34 = cos(t36);
t110 = (mrSges(6,2) - mrSges(7,3)) * t34;
t33 = sin(t36);
t109 = (-mrSges(4,1) - mrSges(5,1)) * t34 + (mrSges(4,2) + t112) * t33;
t108 = -m(5) - m(6);
t22 = t33 * qJ(4);
t42 = cos(qJ(1));
t84 = t34 * t42;
t106 = pkin(3) * t84 + t42 * t22;
t103 = t33 * pkin(5) + t34 * pkin(9);
t39 = sin(qJ(1));
t102 = g(1) * t42 + g(2) * t39;
t74 = qJ(4) * t34;
t7 = t39 * t74;
t40 = cos(qJ(6));
t90 = mrSges(7,1) * t40;
t71 = t34 * t90;
t88 = t33 * t39;
t101 = -m(7) * t7 - mrSges(6,2) * t88 + (t111 * t34 - t71) * t39;
t87 = t33 * t42;
t9 = t42 * t74;
t100 = -m(7) * t9 - mrSges(6,2) * t87 + t111 * t84 - t42 * t71;
t95 = -pkin(3) - pkin(4);
t51 = m(7) * (-pkin(9) + t95) - mrSges(7,3);
t37 = sin(qJ(6));
t83 = t37 * mrSges(7,2);
t45 = t33 * t51 - t34 * t83;
t68 = t95 * t33;
t38 = sin(qJ(2));
t94 = pkin(2) * t38;
t99 = m(7) * t94 - t45 - m(5) * (-pkin(3) * t33 - t94) + t33 * mrSges(5,1) - m(6) * (t68 - t94);
t98 = -(-t83 + t90) * t33 + t109 + t110;
t41 = cos(qJ(2));
t56 = t41 * mrSges(3,1) - t38 * mrSges(3,2);
t97 = -m(3) * pkin(1) - mrSges(2,1) + t109 - t56;
t43 = -pkin(8) - pkin(7);
t96 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) - t107 * (-qJ(5) - t43);
t91 = g(3) * t34;
t31 = t34 * pkin(3);
t35 = t41 * pkin(2);
t82 = t37 * t39;
t81 = t37 * t42;
t80 = t39 * t40;
t79 = t40 * t42;
t75 = t31 + t22;
t30 = t34 * pkin(4);
t70 = t30 + t75;
t69 = t35 + t75;
t32 = t35 + pkin(1);
t10 = t42 * t32;
t64 = -t39 * t43 + t10;
t62 = -t32 - t22;
t59 = t70 + t103;
t53 = mrSges(4,1) * t33 + mrSges(4,2) * t34;
t4 = t33 * t79 - t82;
t3 = -t33 * t81 - t80;
t2 = -t33 * t80 - t81;
t1 = t33 * t82 - t79;
t5 = [(-m(4) * t64 - m(5) * (t64 + t106) - t4 * mrSges(7,1) - t3 * mrSges(7,2) - t107 * (pkin(4) * t84 + t10 + t106) + t96 * t39 + (-m(7) * t103 + t110 + t97) * t42) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + ((m(4) + m(5)) * t43 + t96) * t42 + (m(4) * t32 - m(5) * (t62 - t31) - m(6) * t62 - m(7) * (-t32 + (-pkin(5) - qJ(4)) * t33) + (-m(6) * t95 - mrSges(6,2) - t51) * t34 - t97) * t39) * g(1) (t108 * t7 + t99 * t39 + t101) * g(2) + (t108 * t9 + t99 * t42 + t100) * g(1) + (-t56 - m(4) * t35 - m(5) * t69 - m(6) * (t30 + t69) - m(7) * (t35 + t59) + t98) * g(3) + (m(4) * t94 + mrSges(3,1) * t38 + mrSges(3,2) * t41 + t53) * t102, t102 * t53 + (-m(5) * (-pkin(3) * t88 + t7) + mrSges(5,1) * t88 - m(6) * (t39 * t68 + t7) - t39 * t45 + t101) * g(2) + (-m(5) * (-pkin(3) * t87 + t9) + mrSges(5,1) * t87 - m(6) * (t42 * t68 + t9) - t42 * t45 + t100) * g(1) + (-m(5) * t75 - m(6) * t70 - m(7) * t59 + t98) * g(3) (-t102 * t33 + t91) * (m(5) + t107) t107 * (g(1) * t39 - g(2) * t42) -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - (mrSges(7,1) * t37 + mrSges(7,2) * t40) * t91];
taug  = t5(:);
