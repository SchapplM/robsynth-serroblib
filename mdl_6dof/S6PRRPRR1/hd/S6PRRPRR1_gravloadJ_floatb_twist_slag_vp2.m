% Calculate Gravitation load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:14:23
% EndTime: 2018-11-23 15:14:24
% DurationCPUTime: 0.81s
% Computational Cost: add. (1137->126), mult. (1054->165), div. (0->0), fcn. (997->18), ass. (0->68)
t127 = mrSges(6,2) - mrSges(7,3);
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t126 = -t67 * mrSges(7,1) + t64 * mrSges(7,2) - mrSges(6,1);
t124 = -m(6) - m(7);
t123 = -m(5) * pkin(3) - mrSges(4,1);
t101 = pkin(6) + qJ(2);
t89 = cos(t101) / 0.2e1;
t102 = pkin(6) - qJ(2);
t92 = cos(t102);
t44 = t89 - t92 / 0.2e1;
t59 = qJ(3) + pkin(12);
t56 = qJ(5) + t59;
t51 = sin(t56);
t52 = cos(t56);
t62 = cos(pkin(6));
t24 = t44 * t51 + t52 * t62;
t25 = -t44 * t52 + t51 * t62;
t122 = t126 * t24 + t127 * t25;
t60 = sin(pkin(11));
t61 = sin(pkin(6));
t113 = t60 * t61;
t90 = sin(t101);
t87 = t90 / 0.2e1;
t91 = sin(t102);
t79 = t87 - t91 / 0.2e1;
t103 = cos(pkin(11));
t69 = cos(qJ(2));
t93 = t103 * t69;
t34 = -t60 * t79 + t93;
t13 = t113 * t52 - t34 * t51;
t14 = t113 * t51 + t34 * t52;
t121 = t126 * t13 + t127 * t14;
t112 = t60 * t69;
t31 = t103 * t79 + t112;
t94 = t61 * t103;
t11 = -t31 * t51 - t52 * t94;
t12 = t31 * t52 - t51 * t94;
t120 = t126 * t11 + t127 * t12;
t119 = m(5) - t124;
t54 = sin(t59);
t55 = cos(t59);
t68 = cos(qJ(3));
t57 = t68 * pkin(3);
t65 = sin(qJ(3));
t118 = mrSges(3,1) + m(5) * (t57 + pkin(2)) + t55 * mrSges(5,1) - t54 * mrSges(5,2) + m(4) * pkin(2) + t68 * mrSges(4,1) - t65 * mrSges(4,2) + (m(7) * pkin(5) - t126) * t52 + (m(7) * pkin(10) - t127) * t51;
t63 = -qJ(4) - pkin(8);
t117 = -m(4) * pkin(8) + m(5) * t63 - t64 * mrSges(7,1) - t67 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t45 = -pkin(3) * t65 - pkin(4) * t54;
t46 = pkin(4) * t55 + t57;
t107 = t46 * t113 + t34 * t45;
t104 = -t44 * t45 + t62 * t46;
t98 = t11 * pkin(5) + t12 * pkin(10);
t96 = t13 * pkin(5) + pkin(10) * t14;
t95 = t24 * pkin(5) + pkin(10) * t25;
t88 = t91 / 0.2e1;
t80 = t88 - t90 / 0.2e1;
t77 = t31 * t45 - t46 * t94;
t71 = t92 / 0.2e1 + t89;
t66 = sin(qJ(2));
t58 = -pkin(9) + t63;
t43 = t87 + t88;
t42 = pkin(2) + t46;
t35 = t60 * t80 + t93;
t33 = t103 * t66 + t60 * t71;
t32 = -t103 * t80 + t112;
t30 = -t103 * t71 + t60 * t66;
t1 = [(-m(2) - m(3) - m(4) - t119) * g(3) (t124 * (t43 * t42 + t44 * t58) - t117 * t44 - t118 * t43) * g(3) + (t124 * (-t30 * t42 - t32 * t58) + t117 * t32 + t118 * t30) * g(2) + (t124 * (-t33 * t42 - t35 * t58) + t117 * t35 + t118 * t33) * g(1) (-(t44 * t68 - t62 * t65) * mrSges(4,2) - (t44 * t54 + t55 * t62) * mrSges(5,1) - (t44 * t55 - t54 * t62) * mrSges(5,2) - m(6) * t104 - m(7) * (t95 + t104) + t123 * (t44 * t65 + t62 * t68) + t122) * g(3) + (-(-t31 * t68 + t65 * t94) * mrSges(4,2) - (-t31 * t54 - t55 * t94) * mrSges(5,1) - (-t31 * t55 + t54 * t94) * mrSges(5,2) - m(6) * t77 - m(7) * (t77 + t98) + t123 * (-t31 * t65 - t68 * t94) + t120) * g(2) + (-(-t113 * t65 - t34 * t68) * mrSges(4,2) - (t113 * t55 - t34 * t54) * mrSges(5,1) - (-t113 * t54 - t34 * t55) * mrSges(5,2) - m(6) * t107 - m(7) * (t96 + t107) + t123 * (t113 * t68 - t34 * t65) + t121) * g(1), t119 * (-g(1) * t33 - g(2) * t30 + g(3) * t43) (-m(7) * t95 + t122) * g(3) + (-m(7) * t98 + t120) * g(2) + (-m(7) * t96 + t121) * g(1), -g(1) * ((-t14 * t64 + t33 * t67) * mrSges(7,1) + (-t14 * t67 - t33 * t64) * mrSges(7,2)) - g(2) * ((-t12 * t64 + t30 * t67) * mrSges(7,1) + (-t12 * t67 - t30 * t64) * mrSges(7,2)) - g(3) * ((-t25 * t64 - t43 * t67) * mrSges(7,1) + (-t25 * t67 + t43 * t64) * mrSges(7,2))];
taug  = t1(:);
