% Calculate Gravitation load on the joints for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2018-11-23 17:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:10:35
% EndTime: 2018-11-23 17:10:36
% DurationCPUTime: 0.80s
% Computational Cost: add. (504->102), mult. (478->116), div. (0->0), fcn. (410->10), ass. (0->58)
t116 = mrSges(6,1) + mrSges(7,1);
t44 = cos(qJ(5));
t115 = t116 * t44;
t114 = -mrSges(6,3) - mrSges(7,3);
t107 = mrSges(6,2) + mrSges(7,2);
t38 = qJ(2) + pkin(10);
t35 = qJ(4) + t38;
t29 = sin(t35);
t30 = cos(t35);
t31 = pkin(5) * t44 + pkin(4);
t39 = -qJ(6) - pkin(9);
t113 = m(7) * t30 * t39 + (m(7) * t31 + t115) * t29;
t112 = t29 * t107;
t33 = sin(t38);
t34 = cos(t38);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t111 = -t45 * mrSges(3,1) - t34 * mrSges(4,1) + t42 * mrSges(3,2) + t33 * mrSges(4,2);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t101 = g(1) * t46 + g(2) * t43;
t110 = t30 * mrSges(5,1) + (-mrSges(5,2) - t114) * t29;
t104 = t30 * pkin(4) + t29 * pkin(9);
t106 = -t29 * t39 + t30 * t31;
t109 = -m(6) * t104 - m(7) * t106;
t94 = m(7) * pkin(5);
t41 = sin(qJ(5));
t82 = t41 * t46;
t83 = t30 * t46;
t103 = -t82 * t112 + t113 * t46 + t114 * t83;
t81 = t43 * t41;
t84 = t30 * t43;
t102 = -t81 * t112 + t113 * t43 + t114 * t84;
t100 = m(5) + m(6) + m(7);
t99 = t94 + t116;
t98 = -t110 + (t107 * t41 - t115) * t30;
t40 = -qJ(3) - pkin(7);
t96 = -m(3) * pkin(7) + m(4) * t40 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t36 = t45 * pkin(2);
t95 = mrSges(2,1) + m(4) * (t36 + pkin(1)) + m(3) * pkin(1) + t110 - t111;
t93 = pkin(2) * t42;
t92 = pkin(4) * t29;
t80 = t43 * t44;
t78 = t44 * t46;
t76 = pkin(3) * t34 + t36;
t66 = pkin(9) * t84 - t43 * t92;
t65 = pkin(9) * t83 - t46 * t92;
t59 = mrSges(5,1) * t29 + mrSges(5,2) * t30;
t3 = -t30 * t82 + t80;
t1 = t30 * t81 + t78;
t37 = -pkin(8) + t40;
t14 = -pkin(3) * t33 - t93;
t13 = pkin(1) + t76;
t7 = t46 * t14;
t6 = t43 * t14;
t4 = t30 * t78 + t81;
t2 = -t30 * t80 + t82;
t5 = [(-t81 * t94 - t100 * (t46 * t13 - t43 * t37) - t116 * t4 - t107 * t3 + t96 * t43 + (t109 - t95) * t46) * g(2) + (-t116 * t2 - t107 * t1 + (t100 * t37 - t41 * t94 + t96) * t46 + (m(5) * t13 - m(6) * (-t13 - t104) - m(7) * (-t13 - t106) + t95) * t43) * g(1) (-m(6) * (t6 + t66) - m(7) * t6 + t102) * g(2) + (-m(6) * (t65 + t7) - m(7) * t7 + t103) * g(1) + (-m(4) * t36 - m(5) * t76 - m(6) * (t76 + t104) - m(7) * (t106 + t76) + t98 + t111) * g(3) + t101 * (m(4) * t93 - m(5) * t14 + mrSges(3,1) * t42 + mrSges(4,1) * t33 + mrSges(3,2) * t45 + mrSges(4,2) * t34 + t59) (-g(1) * t43 + g(2) * t46) * (m(4) + t100) t101 * t59 + (-m(6) * t66 + t102) * g(2) + (-m(6) * t65 + t103) * g(1) + (t109 + t98) * g(3) (t107 * t44 + t99 * t41) * g(3) * t29 + (t99 * t1 - t107 * t2) * g(2) + (t107 * t4 - t99 * t3) * g(1) (g(3) * t30 - t101 * t29) * m(7)];
taug  = t5(:);
