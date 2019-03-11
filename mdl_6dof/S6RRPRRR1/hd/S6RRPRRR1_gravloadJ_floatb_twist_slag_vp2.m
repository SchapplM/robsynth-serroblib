% Calculate Gravitation load on the joints for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:28
% EndTime: 2019-03-09 13:12:29
% DurationCPUTime: 0.56s
% Computational Cost: add. (568->108), mult. (430->111), div. (0->0), fcn. (354->12), ass. (0->59)
t36 = qJ(2) + pkin(11);
t32 = qJ(4) + t36;
t28 = qJ(5) + t32;
t23 = sin(t28);
t24 = cos(t28);
t38 = sin(qJ(6));
t81 = t38 * mrSges(7,2);
t107 = t23 * t81 + t24 * (m(7) * pkin(10) + mrSges(7,3));
t98 = t24 * pkin(5) + t23 * pkin(10);
t106 = m(7) * t98;
t105 = -t24 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t23;
t26 = sin(t32);
t27 = cos(t32);
t86 = mrSges(6,2) * t24;
t104 = mrSges(5,1) * t26 + mrSges(6,1) * t23 + mrSges(5,2) * t27 + t86;
t30 = sin(t36);
t31 = cos(t36);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t103 = -t42 * mrSges(3,1) - t31 * mrSges(4,1) + t39 * mrSges(3,2) + t30 * mrSges(4,2);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t102 = g(1) * t43 + g(2) * t40;
t99 = m(6) + m(7);
t65 = t27 * mrSges(5,1) - t26 * mrSges(5,2);
t41 = cos(qJ(6));
t78 = t41 * mrSges(7,1);
t97 = -(t78 - t81) * t24 + t105;
t95 = -t65 + t97;
t37 = -qJ(3) - pkin(7);
t35 = -pkin(8) + t37;
t93 = -m(3) * pkin(7) + m(4) * t37 + m(5) * t35 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t34 = t42 * pkin(2);
t74 = pkin(3) * t31 + t34;
t92 = mrSges(2,1) + m(5) * (pkin(1) + t74) + t65 + m(4) * (t34 + pkin(1)) + m(3) * pkin(1) - t103 - t105;
t91 = pkin(4) * t26;
t22 = pkin(4) * t27;
t90 = pkin(5) * t23;
t87 = t39 * pkin(2);
t80 = t40 * t38;
t79 = t40 * t41;
t77 = t43 * t38;
t76 = t43 * t41;
t72 = t23 * t78;
t71 = t22 + t74;
t63 = t107 * t40;
t62 = t107 * t43;
t16 = -pkin(3) * t30 - t87;
t7 = t16 - t91;
t46 = m(7) * (t7 - t90) - t72;
t45 = m(7) * (-t90 - t91) - t72;
t44 = t86 + (m(7) * pkin(5) + mrSges(6,1) + t78) * t23;
t33 = -pkin(9) + t35;
t6 = t24 * t76 + t80;
t5 = -t24 * t77 + t79;
t4 = -t24 * t79 + t77;
t3 = t24 * t80 + t76;
t2 = pkin(1) + t71;
t1 = [(-t6 * mrSges(7,1) - t5 * mrSges(7,2) - t99 * (t43 * t2 - t40 * t33) + t93 * t40 + (-t92 - t106) * t43) * g(2) + (-t4 * mrSges(7,1) - t3 * mrSges(7,2) + (t99 * t33 + t93) * t43 + (m(6) * t2 - m(7) * (-t2 - t98) + t92) * t40) * g(1), -g(1) * (t46 * t43 + t62) - g(2) * (t46 * t40 + t63) + (-m(4) * t34 - m(5) * t74 - m(6) * t71 - m(7) * (t71 + t98) + t95 + t103) * g(3) + t102 * (m(4) * t87 - m(5) * t16 - m(6) * t7 + mrSges(3,1) * t39 + mrSges(4,1) * t30 + mrSges(3,2) * t42 + mrSges(4,2) * t31 + t104) (-g(1) * t40 + g(2) * t43) * (m(4) + m(5) + t99) -g(1) * (t45 * t43 + t62) - g(2) * (t45 * t40 + t63) + (-m(6) * t22 - m(7) * (t22 + t98) + t95) * g(3) + (m(6) * t91 + t104) * t102 (t97 - t106) * g(3) + (t44 * t40 - t63) * g(2) + (t44 * t43 - t62) * g(1), -g(1) * (t5 * mrSges(7,1) - t6 * mrSges(7,2)) - g(2) * (-t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t38 - mrSges(7,2) * t41) * t23];
taug  = t1(:);
