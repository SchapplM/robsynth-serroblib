% Calculate Gravitation load on the joints for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:29
% EndTime: 2018-11-23 15:09:29
% DurationCPUTime: 0.79s
% Computational Cost: add. (968->111), mult. (1171->135), div. (0->0), fcn. (1119->14), ass. (0->69)
t100 = m(6) + m(7);
t43 = sin(qJ(6));
t46 = cos(qJ(6));
t103 = t100 * (pkin(8) - qJ(5)) - t43 * mrSges(7,1) - t46 * mrSges(7,2) - mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t90 = m(7) * pkin(9) + mrSges(4,1) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t102 = -t46 * mrSges(7,1) + t43 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t78 = pkin(6) + qJ(2);
t62 = cos(t78) / 0.2e1;
t79 = pkin(6) - qJ(2);
t66 = cos(t79);
t34 = t66 / 0.2e1 + t62;
t42 = sin(pkin(10));
t45 = sin(qJ(2));
t81 = cos(pkin(10));
t19 = t34 * t81 - t42 * t45;
t44 = sin(qJ(3));
t83 = qJ(4) * t44;
t47 = cos(qJ(3));
t89 = t19 * t47;
t99 = pkin(3) * t89 + t19 * t83;
t22 = -t42 * t34 - t45 * t81;
t88 = t22 * t47;
t98 = pkin(3) * t88 + t22 * t83;
t64 = sin(t78);
t37 = t64 / 0.2e1;
t65 = sin(t79);
t61 = t65 / 0.2e1;
t32 = t37 + t61;
t87 = t32 * t47;
t97 = pkin(3) * t87 + t32 * t83;
t96 = m(5) + t100;
t93 = -m(7) * (qJ(4) + pkin(5)) + t102;
t92 = -mrSges(3,1) - t90 * t47 + (-m(7) * pkin(5) + t102) * t44;
t48 = cos(qJ(2));
t86 = t42 * t48;
t84 = t37 - t65 / 0.2e1;
t82 = cos(pkin(6));
t80 = sin(pkin(6));
t15 = t19 * pkin(2);
t50 = t61 - t64 / 0.2e1;
t21 = -t50 * t81 + t86;
t77 = pkin(8) * t21 + t15;
t16 = t22 * pkin(2);
t71 = t81 * t48;
t24 = t42 * t50 + t71;
t76 = pkin(8) * t24 + t16;
t31 = t32 * pkin(2);
t33 = t62 - t66 / 0.2e1;
t75 = -pkin(8) * t33 + t31;
t20 = t81 * t84 + t86;
t55 = t81 * t80;
t5 = t20 * t44 + t47 * t55;
t2 = t5 * pkin(3);
t6 = t20 * t47 - t44 * t55;
t74 = qJ(4) * t6 - t2;
t23 = -t42 * t84 + t71;
t72 = t42 * t80;
t7 = t23 * t44 - t47 * t72;
t4 = t7 * pkin(3);
t8 = t23 * t47 + t44 * t72;
t73 = qJ(4) * t8 - t4;
t25 = -t33 * t44 - t47 * t82;
t18 = t25 * pkin(3);
t26 = -t33 * t47 + t44 * t82;
t69 = qJ(4) * t26 - t18;
t17 = t25 * pkin(4);
t3 = t7 * pkin(4);
t1 = t5 * pkin(4);
t9 = [(-m(2) - m(3) - m(4) - t96) * g(3) (-m(4) * t75 - m(5) * (t75 + t97) - t100 * (pkin(4) * t87 + t31 + t97) + t103 * t33 + t92 * t32) * g(3) + (-m(4) * t77 - m(5) * (t77 + t99) - t100 * (pkin(4) * t89 + t15 + t99) - t103 * t21 + t92 * t19) * g(2) + (-m(4) * t76 - m(5) * (t76 + t98) - t100 * (pkin(4) * t88 + t16 + t98) - t103 * t24 + t92 * t22) * g(1) (-m(5) * t69 - m(6) * (-t17 + t69) - m(7) * (-t17 - t18) + t93 * t26 + t90 * t25) * g(3) + (-m(5) * t74 - m(6) * (-t1 + t74) - m(7) * (-t1 - t2) + t93 * t6 + t90 * t5) * g(2) + (-m(5) * t73 - m(6) * (-t3 + t73) - m(7) * (-t3 - t4) + t93 * t8 + t90 * t7) * g(1), t96 * (-g(1) * t7 - g(2) * t5 - g(3) * t25) t100 * (-g(1) * t22 - g(2) * t19 - g(3) * t32) -g(1) * ((t22 * t46 - t43 * t7) * mrSges(7,1) + (-t22 * t43 - t46 * t7) * mrSges(7,2)) - g(2) * ((t19 * t46 - t43 * t5) * mrSges(7,1) + (-t19 * t43 - t46 * t5) * mrSges(7,2)) - g(3) * ((-t25 * t43 + t32 * t46) * mrSges(7,1) + (-t25 * t46 - t32 * t43) * mrSges(7,2))];
taug  = t9(:);
