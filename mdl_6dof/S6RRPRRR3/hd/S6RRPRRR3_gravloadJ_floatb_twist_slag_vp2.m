% Calculate Gravitation load on the joints for
% S6RRPRRR3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:21:45
% EndTime: 2018-11-23 17:21:46
% DurationCPUTime: 0.76s
% Computational Cost: add. (550->125), mult. (547->143), div. (0->0), fcn. (499->12), ass. (0->74)
t39 = qJ(4) + qJ(5);
t33 = cos(t39);
t44 = cos(qJ(4));
t35 = t44 * pkin(4);
t23 = pkin(5) * t33 + t35;
t21 = pkin(3) + t23;
t34 = qJ(6) + t39;
t26 = sin(t34);
t27 = cos(t34);
t28 = t35 + pkin(3);
t32 = sin(t39);
t41 = sin(qJ(4));
t108 = -m(5) * pkin(3) - m(6) * t28 - m(7) * t21 - t44 * mrSges(5,1) - t33 * mrSges(6,1) - t27 * mrSges(7,1) + t41 * mrSges(5,2) + t32 * mrSges(6,2) + t26 * mrSges(7,2);
t99 = m(4) + m(5) + m(6) + m(7);
t107 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t37 = qJ(2) + pkin(11);
t30 = sin(t37);
t31 = cos(t37);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t106 = -t45 * mrSges(3,1) - t31 * mrSges(4,1) + t42 * mrSges(3,2) - t107 * t30;
t94 = m(6) * pkin(4);
t89 = pkin(4) * t41;
t105 = m(6) * t89;
t88 = pkin(5) * t32;
t22 = t88 + t89;
t104 = m(7) * t22;
t103 = mrSges(5,1) + t94;
t58 = -mrSges(7,1) * t26 - mrSges(7,2) * t27;
t102 = mrSges(6,1) * t32 + mrSges(6,2) * t33 - t58;
t43 = sin(qJ(1));
t73 = t33 * t43;
t46 = cos(qJ(1));
t74 = t32 * t46;
t15 = -t31 * t74 + t73;
t72 = t33 * t46;
t75 = t32 * t43;
t16 = t31 * t72 + t75;
t76 = t31 * t46;
t7 = -t26 * t76 + t27 * t43;
t8 = t26 * t43 + t27 * t76;
t91 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t101 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t91;
t13 = t31 * t75 + t72;
t14 = -t31 * t73 + t74;
t77 = t31 * t43;
t5 = t26 * t77 + t27 * t46;
t6 = t26 * t46 - t27 * t77;
t92 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t100 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t92;
t96 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t104;
t95 = m(3) * pkin(1) + mrSges(2,1) - t106;
t93 = m(7) * pkin(5);
t47 = -pkin(9) - pkin(8);
t87 = pkin(8) * t30;
t84 = g(3) * t30;
t36 = t45 * pkin(2);
t38 = -pkin(10) + t47;
t79 = t30 * t38;
t78 = t30 * t47;
t71 = t41 * t43;
t70 = t41 * t46;
t69 = t43 * t44;
t68 = t44 * t46;
t63 = pkin(3) * t31 + t87;
t57 = t21 * t31 - t79;
t56 = t31 * t28 - t78;
t19 = -t31 * t70 + t69;
t17 = t31 * t71 + t68;
t40 = -qJ(3) - pkin(7);
t29 = t36 + pkin(1);
t20 = t31 * t68 + t71;
t18 = -t31 * t69 + t70;
t1 = [(-t71 * t94 - t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t19 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) - t99 * (t46 * t29 - t43 * t40) + t96 * t43 + (-m(5) * t63 - m(6) * t56 - m(7) * t57 - t95) * t46) * g(2) + (-t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (t99 * t40 - t105 + t96) * t46 + (m(4) * t29 - m(5) * (-t29 - t63) - m(6) * (-t29 - t56) - m(7) * (-t29 - t57) + t95) * t43) * g(1) (-m(4) * t36 - m(5) * (t36 + t87) - m(6) * (t36 - t78) - m(7) * (t36 - t79) + t108 * t31 + t106) * g(3) + (g(1) * t46 + g(2) * t43) * (mrSges(3,2) * t45 + (-m(5) * pkin(8) + m(6) * t47 + m(7) * t38 - t107) * t31 + (mrSges(4,1) - t108) * t30 + (t99 * pkin(2) + mrSges(3,1)) * t42) (-g(1) * t43 + g(2) * t46) * t99 (mrSges(5,1) * t41 + mrSges(5,2) * t44 + t102 + t104 + t105) * t84 + (-t18 * mrSges(5,2) - m(7) * (-t22 * t77 - t23 * t46) + t103 * t17 + t100) * g(2) + (t20 * mrSges(5,2) - m(7) * (-t22 * t76 + t23 * t43) - t103 * t19 + t101) * g(1) (m(7) * t88 + t102) * t84 + (t13 * t93 + t100) * g(2) + (-t15 * t93 + t101) * g(1), -g(1) * t91 - g(2) * t92 - t58 * t84];
taug  = t1(:);
