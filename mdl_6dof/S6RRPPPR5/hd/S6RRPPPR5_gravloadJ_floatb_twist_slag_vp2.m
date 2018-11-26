% Calculate Gravitation load on the joints for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2018-11-23 16:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:44:31
% EndTime: 2018-11-23 16:44:32
% DurationCPUTime: 1.03s
% Computational Cost: add. (295->119), mult. (688->138), div. (0->0), fcn. (681->8), ass. (0->62)
t98 = mrSges(4,1) - mrSges(5,2);
t87 = mrSges(4,2) - mrSges(6,1) - mrSges(5,3);
t97 = mrSges(7,3) - mrSges(6,2);
t83 = -m(7) * pkin(8) - t97;
t95 = mrSges(5,1) + mrSges(4,3);
t96 = t83 - t95;
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t46 = t27 * t29 - t28 * t32;
t47 = t27 * t32 + t28 * t29;
t94 = t47 * mrSges(7,1) - t46 * mrSges(7,2) - t87 * t27 + t98 * t28;
t93 = -m(4) - m(5);
t92 = m(6) + m(7);
t91 = mrSges(2,2) - mrSges(3,3);
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t90 = g(1) * t34 + g(2) * t31;
t89 = m(5) + t92;
t88 = mrSges(6,3) + t98;
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t51 = t33 * mrSges(3,1) - t30 * mrSges(3,2);
t86 = t30 * t95 + t51;
t62 = qJ(4) * t27;
t58 = -pkin(2) - t62;
t81 = pkin(3) * t28;
t43 = t58 - t81;
t61 = qJ(5) * t28;
t80 = pkin(5) * t27;
t84 = t96 * t33 + (-m(7) * (t43 - t61 - t80) - m(6) * t58 - (m(6) * (-pkin(3) - qJ(5)) - mrSges(6,3)) * t28 - m(5) * t43 + m(4) * pkin(2) + t94) * t30;
t77 = g(3) * t30;
t24 = t33 * pkin(2);
t73 = t30 * t34;
t72 = t31 * t33;
t71 = t33 * t34;
t70 = t34 * t27;
t7 = -t31 * t28 + t33 * t70;
t69 = t7 * qJ(4);
t68 = -pkin(4) - qJ(3);
t21 = t30 * qJ(3);
t65 = t24 + t21;
t64 = t34 * pkin(1) + t31 * pkin(7);
t63 = qJ(3) * t33;
t59 = -pkin(1) - t24;
t57 = t65 + (t62 + t81) * t33;
t56 = pkin(2) * t71 + t34 * t21 + t64;
t25 = t34 * pkin(7);
t5 = t27 * t72 + t28 * t34;
t6 = t28 * t72 - t70;
t55 = -t6 * pkin(3) - qJ(4) * t5 + t25;
t8 = t31 * t27 + t28 * t71;
t54 = t8 * pkin(3) + t56;
t53 = -t29 * t6 - t32 * t5;
t52 = t29 * t5 - t32 * t6;
t39 = pkin(4) * t73 + t8 * qJ(5) + t54;
t13 = t34 * t63;
t11 = t31 * t63;
t2 = t29 * t8 + t32 * t7;
t1 = -t29 * t7 + t32 * t8;
t3 = [(-m(3) * t64 - m(4) * t56 - m(5) * (t54 + t69) - m(6) * (t39 + t69) - m(7) * t39 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(2,1) - t51) * t34 + t91 * t31 - t88 * t8 + (-m(7) * (pkin(5) + qJ(4)) + t87) * t7 + t96 * t73) * g(2) + (-m(5) * t55 - t53 * mrSges(7,1) - t52 * mrSges(7,2) - t92 * (-qJ(5) * t6 + t55) + t91 * t34 + (-m(3) - m(4)) * t25 + t88 * t6 + (m(7) * pkin(5) - t87) * t5 + (m(3) * pkin(1) + mrSges(2,1) - t92 * t59 + t93 * (t59 - t21) + (-m(6) * t68 - m(7) * (-pkin(8) + t68) + t97) * t30 + t86) * t31) * g(1), t90 * (mrSges(3,1) * t30 + mrSges(3,2) * t33) + (-t92 * (pkin(4) * t72 + t11) + t93 * t11 + t84 * t31) * g(2) + (-t92 * (pkin(4) * t71 + t13) + t93 * t13 + t84 * t34) * g(1) + (-m(4) * t65 - m(5) * t57 - t92 * (t30 * pkin(4) + t33 * t61 + t57) + t83 * t30 + (-m(7) * t80 - t28 * mrSges(6,3) - t94) * t33 - t86) * g(3) (t33 * g(3) - t90 * t30) * (m(4) + t89) t89 * (-g(1) * t7 - g(2) * t5 - t27 * t77) t92 * (-g(1) * t8 - g(2) * t6 - t28 * t77) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t52 * mrSges(7,1) + t53 * mrSges(7,2)) - (-t46 * mrSges(7,1) - t47 * mrSges(7,2)) * t77];
taug  = t3(:);
