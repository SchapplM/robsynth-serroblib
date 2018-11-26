% Calculate Gravitation load on the joints for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:26:32
% EndTime: 2018-11-23 16:26:33
% DurationCPUTime: 0.72s
% Computational Cost: add. (535->101), mult. (501->114), div. (0->0), fcn. (447->10), ass. (0->55)
t111 = mrSges(6,1) + mrSges(7,1);
t109 = -mrSges(6,3) - mrSges(7,2);
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t110 = t42 * mrSges(7,3) + t111 * t44;
t104 = m(6) + m(7);
t108 = m(5) + t104;
t38 = pkin(10) + qJ(3);
t35 = qJ(4) + t38;
t28 = sin(t35);
t107 = t109 * t28;
t106 = pkin(5) * t44 + qJ(6) * t42;
t100 = (-m(7) * (-pkin(4) - t106) + t110) * t28;
t29 = cos(t35);
t103 = t29 * mrSges(5,1) - t28 * mrSges(5,2);
t71 = t29 * pkin(4) + t28 * pkin(9);
t91 = pkin(4) * t28;
t33 = sin(t38);
t92 = pkin(3) * t33;
t102 = m(7) * t92 - m(6) * (-t91 - t92) + t100;
t45 = cos(qJ(1));
t80 = t29 * t45;
t20 = pkin(9) * t80;
t83 = t28 * t42;
t69 = mrSges(6,2) * t83;
t101 = -m(7) * t20 + t109 * t80 - t45 * t69;
t43 = sin(qJ(1));
t99 = g(1) * t45 + g(2) * t43;
t98 = -t103 + t107 + (mrSges(6,2) * t42 - t110) * t29;
t40 = cos(pkin(10));
t30 = t40 * pkin(2) + pkin(1);
t34 = cos(t38);
t58 = t34 * mrSges(4,1) - t33 * mrSges(4,2);
t97 = -m(4) * t30 - mrSges(2,1) - m(3) * pkin(1) - t40 * mrSges(3,1) + sin(pkin(10)) * mrSges(3,2) - t103 - t58;
t41 = -pkin(7) - qJ(2);
t96 = -m(3) * qJ(2) + m(4) * t41 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t95 = m(7) * pkin(5) + t111;
t94 = (-t69 + (-t104 * pkin(9) + t109) * t29) * t43;
t93 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t27 = pkin(3) * t34;
t82 = t28 * t45;
t77 = t43 * t42;
t76 = t43 * t44;
t73 = t44 * t45;
t72 = t45 * t42;
t37 = -pkin(8) + t41;
t8 = t27 + t30;
t63 = -t37 * t43 + t45 * t8;
t60 = t106 * t29 + t71;
t55 = mrSges(5,1) * t28 + mrSges(5,2) * t29;
t4 = t29 * t73 + t77;
t3 = t29 * t72 - t76;
t2 = t29 * t76 - t72;
t1 = t29 * t77 + t73;
t5 = [(-m(5) * t63 + t109 * t82 - t104 * (pkin(4) * t80 + pkin(9) * t82 + t63) - t95 * t4 - t93 * t3 + t97 * t45 + t96 * t43) * g(2) + (t95 * t2 + t93 * t1 + (m(5) * t8 - t104 * (-t8 - t71) - t97 - t107) * t43 + (t108 * t37 + t96) * t45) * g(1) (-g(1) * t43 + g(2) * t45) * (m(3) + m(4) + t108) (t102 * t43 + t94) * g(2) + (-m(6) * t20 + t102 * t45 + t101) * g(1) + (-t58 - m(5) * t27 - m(6) * (t27 + t71) - m(7) * (t27 + t60) + t98) * g(3) + (m(5) * t92 + mrSges(4,1) * t33 + mrSges(4,2) * t34 + t55) * t99, t99 * t55 + ((m(6) * t91 + t100) * t43 + t94) * g(2) + (-m(6) * (-pkin(4) * t82 + t20) + t100 * t45 + t101) * g(1) + (-m(6) * t71 - m(7) * t60 + t98) * g(3) (t95 * t42 - t93 * t44) * g(3) * t28 + (t95 * t1 - t93 * t2) * g(2) + (t95 * t3 - t93 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t83) * m(7)];
taug  = t5(:);
