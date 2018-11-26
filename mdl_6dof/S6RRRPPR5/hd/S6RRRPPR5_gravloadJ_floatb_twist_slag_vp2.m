% Calculate Gravitation load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2018-11-23 17:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:35:29
% EndTime: 2018-11-23 17:35:31
% DurationCPUTime: 1.47s
% Computational Cost: add. (1489->128), mult. (1564->168), div. (0->0), fcn. (1527->18), ass. (0->60)
t54 = pkin(12) + qJ(6);
t49 = sin(t54);
t51 = cos(t54);
t56 = sin(pkin(12));
t58 = cos(pkin(12));
t68 = -mrSges(5,1) - m(6) * pkin(4) - t58 * mrSges(6,1) + t56 * mrSges(6,2) - m(7) * (pkin(5) * t58 + pkin(4));
t141 = -mrSges(7,1) * t51 + mrSges(7,2) * t49 + t68;
t126 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t55 = qJ(3) + pkin(11);
t50 = sin(t55);
t52 = cos(t55);
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t123 = m(4) * pkin(2) + t64 * mrSges(4,1) - t61 * mrSges(4,2) - t126 * t50 - t141 * t52 + mrSges(3,1);
t124 = -m(4) * pkin(9) - t58 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t56 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t139 = mrSges(7,1) * t49 + mrSges(7,2) * t51 - t124;
t133 = m(6) + m(7);
t129 = m(5) + t133;
t106 = cos(pkin(6));
t104 = pkin(6) + qJ(2);
t87 = cos(t104) / 0.2e1;
t105 = pkin(6) - qJ(2);
t93 = cos(t105);
t33 = t87 - t93 / 0.2e1;
t135 = t106 * t64 + t33 * t61;
t57 = sin(pkin(6));
t63 = sin(qJ(1));
t112 = t57 * t63;
t118 = cos(qJ(1));
t86 = sin(t104) / 0.2e1;
t92 = sin(t105);
t32 = t86 - t92 / 0.2e1;
t65 = cos(qJ(2));
t80 = t118 * t65 - t32 * t63;
t9 = t112 * t64 - t61 * t80;
t81 = t118 * t32 + t63 * t65;
t99 = t57 * t118;
t75 = t61 * t81 + t64 * t99;
t107 = pkin(1) * t118 + pkin(8) * t112;
t103 = t61 * t112;
t97 = -pkin(1) * t63 + pkin(8) * t99;
t4 = -t50 * t99 + t52 * t81;
t39 = t61 * t99;
t95 = -t64 * t81 + t39;
t3 = t50 * t81 + t52 * t99;
t69 = t93 / 0.2e1 + t87;
t62 = sin(qJ(2));
t59 = -qJ(4) - pkin(9);
t48 = pkin(3) * t64 + pkin(2);
t31 = t86 + t92 / 0.2e1;
t24 = t118 * t62 + t63 * t69;
t21 = -t118 * t69 + t62 * t63;
t16 = t106 * t50 - t33 * t52;
t15 = -t106 * t52 - t33 * t50;
t10 = t64 * t80 + t103;
t8 = t112 * t50 + t52 * t80;
t7 = -t112 * t52 + t50 * t80;
t2 = t24 * t49 + t51 * t8;
t1 = t24 * t51 - t49 * t8;
t5 = [(-t118 * mrSges(2,1) - m(3) * t107 - t80 * mrSges(3,1) - m(4) * (pkin(2) * t80 + t107) - t10 * mrSges(4,1) - t9 * mrSges(4,2) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t57 + mrSges(2,2)) * t63 + t68 * t8 + t126 * t7 + t124 * t24 - t129 * (pkin(3) * t103 - t24 * t59 + t48 * t80 + t107)) * g(2) + (t63 * mrSges(2,1) + t118 * mrSges(2,2) - m(3) * t97 + t81 * mrSges(3,1) - mrSges(3,3) * t99 - m(4) * (-pkin(2) * t81 + t97) - t95 * mrSges(4,1) - t75 * mrSges(4,2) - t126 * t3 - t141 * t4 + t139 * t21 + t129 * (-pkin(3) * t39 - t21 * t59 + t48 * t81 - t97)) * g(1) (-t129 * (t31 * t48 + t33 * t59) + t139 * t33 - t123 * t31) * g(3) + (-t129 * (-t21 * t48 - t59 * t81) - t139 * t81 + t123 * t21) * g(2) + (-t129 * (-t24 * t48 - t59 * t80) - t139 * t80 + t123 * t24) * g(1) (-t135 * mrSges(4,1) - (-t106 * t61 + t33 * t64) * mrSges(4,2) + t126 * t16 - t141 * t15) * g(3) + (mrSges(4,1) * t75 - mrSges(4,2) * t95 + t126 * t4 - t141 * t3) * g(2) + (-mrSges(4,1) * t9 + mrSges(4,2) * t10 + t126 * t8 - t141 * t7) * g(1) + (-g(1) * t9 + g(2) * t75 - g(3) * t135) * t129 * pkin(3), t129 * (-g(1) * t24 - g(2) * t21 + g(3) * t31) t133 * (-g(1) * t7 - g(2) * t3 - g(3) * t15) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t51 - t4 * t49) * mrSges(7,1) + (-t21 * t49 - t4 * t51) * mrSges(7,2)) - g(3) * ((-t16 * t49 - t31 * t51) * mrSges(7,1) + (-t16 * t51 + t31 * t49) * mrSges(7,2))];
taug  = t5(:);
