% Calculate Gravitation load on the joints for
% S6RPRRRP7
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
% Datum: 2018-11-23 16:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:27:35
% EndTime: 2018-11-23 16:27:36
% DurationCPUTime: 0.90s
% Computational Cost: add. (525->117), mult. (553->131), div. (0->0), fcn. (519->10), ass. (0->63)
t113 = mrSges(6,3) + mrSges(7,2);
t36 = pkin(10) + qJ(3);
t32 = sin(t36);
t107 = t113 * t32;
t112 = m(6) + m(7);
t43 = cos(qJ(4));
t31 = pkin(4) * t43 + pkin(3);
t33 = cos(t36);
t21 = t33 * t31;
t45 = -pkin(9) - pkin(8);
t80 = t32 * t45;
t114 = -t112 * (t21 - t80) - t107;
t111 = mrSges(6,1) + mrSges(7,1);
t110 = mrSges(6,2) - mrSges(7,3);
t37 = qJ(4) + qJ(5);
t34 = sin(t37);
t35 = cos(t37);
t41 = sin(qJ(4));
t108 = m(5) * pkin(3) + t43 * mrSges(5,1) - t41 * mrSges(5,2) - t110 * t34 + t111 * t35;
t42 = sin(qJ(1));
t73 = t42 * t43;
t44 = cos(qJ(1));
t75 = t41 * t44;
t17 = -t33 * t75 + t73;
t106 = -m(7) * qJ(6) - mrSges(7,3);
t104 = -m(4) - m(5);
t71 = t44 * t34;
t74 = t42 * t35;
t13 = t33 * t71 - t74;
t78 = t35 * t44;
t79 = t34 * t42;
t14 = t33 * t78 + t79;
t101 = t110 * t14 + t111 * t13;
t11 = t33 * t79 + t78;
t12 = t33 * t74 - t71;
t100 = t11 * t111 + t110 * t12;
t99 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t39 = cos(pkin(10));
t58 = t33 * mrSges(4,1) - t32 * mrSges(4,2);
t98 = mrSges(2,1) + m(3) * pkin(1) + t39 * mrSges(3,1) - sin(pkin(10)) * mrSges(3,2) + t58 + t32 * mrSges(5,3);
t96 = m(7) * pkin(5) + t111;
t95 = -mrSges(6,2) - t106;
t90 = pkin(4) * t41;
t87 = g(3) * t32;
t86 = mrSges(6,2) * t35;
t40 = -pkin(7) - qJ(2);
t77 = t40 * t44;
t76 = t41 * t42;
t72 = t43 * t44;
t69 = m(5) * pkin(8) + mrSges(5,3);
t68 = t106 * t32 * t35;
t66 = -pkin(5) * t11 + qJ(6) * t12;
t30 = pkin(2) * t39 + pkin(1);
t64 = t30 * t44 - t40 * t42;
t62 = -pkin(5) * t13 + qJ(6) * t14;
t59 = pkin(3) * t33 + pkin(8) * t32;
t54 = pkin(5) * t35 + qJ(6) * t34;
t53 = t17 * pkin(4);
t15 = t33 * t76 + t72;
t50 = t15 * pkin(4);
t18 = t33 * t72 + t76;
t16 = -t33 * t73 + t75;
t1 = [(-t18 * mrSges(5,1) - t17 * mrSges(5,2) + t104 * t64 - t112 * (pkin(4) * t76 + t64) - t96 * t14 - t95 * t13 + t99 * t42) * g(2) + (m(5) * t77 - t16 * mrSges(5,1) - t15 * mrSges(5,2) - t112 * (pkin(4) * t75 + t42 * t80 - t77) + t96 * t12 + t95 * t11 + (m(4) * t30 - m(5) * (-t30 - t59) - t112 * (-t30 - t21) + t98 + t107) * t42) * g(1) + ((-m(5) * t59 + t114 - t98) * g(2) + (m(4) * t40 + t99) * g(1)) * t44 (-g(1) * t42 + g(2) * t44) * (m(3) + t112 - t104) (-t69 * t32 - t58 + (-m(7) * t54 - t108) * t33 + t114) * g(3) + (g(1) * t44 + g(2) * t42) * ((t112 * t45 + mrSges(4,2) - t113 - t69) * t33 + (mrSges(4,1) + m(6) * t31 - m(7) * (-t31 - t54) + t108) * t32) -g(3) * ((m(7) * (-pkin(5) * t34 - t90) - t34 * mrSges(7,1)) * t32 - t68) + (m(6) * t90 + mrSges(5,1) * t41 + mrSges(6,1) * t34 + mrSges(5,2) * t43 + t86) * t87 + (t15 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * t50 - m(7) * (-t50 + t66) + t100) * g(2) + (-t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t53 - m(7) * (t53 + t62) + t101) * g(1) ((t34 * t96 + t86) * t32 + t68) * g(3) + (-m(7) * t66 + t100) * g(2) + (-m(7) * t62 + t101) * g(1) (-g(1) * t13 - g(2) * t11 - t34 * t87) * m(7)];
taug  = t1(:);
