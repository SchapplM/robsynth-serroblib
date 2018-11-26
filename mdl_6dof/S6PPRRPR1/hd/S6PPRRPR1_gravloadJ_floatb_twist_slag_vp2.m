% Calculate Gravitation load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:44
% EndTime: 2018-11-23 14:49:45
% DurationCPUTime: 0.74s
% Computational Cost: add. (2485->93), mult. (2574->131), div. (0->0), fcn. (2518->24), ass. (0->75)
t110 = m(6) + m(7);
t108 = m(5) + t110;
t41 = pkin(13) + qJ(6);
t39 = sin(t41);
t40 = cos(t41);
t42 = sin(pkin(13));
t43 = cos(pkin(13));
t107 = mrSges(5,1) + m(7) * (pkin(5) * t43 + pkin(4)) + t40 * mrSges(7,1) - t39 * mrSges(7,2) + m(6) * pkin(4) + t43 * mrSges(6,1) - t42 * mrSges(6,2);
t106 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t114 = pkin(3) * t108 - t106 * t46 + t107 * t48 + mrSges(4,1);
t90 = pkin(7) + qJ(3);
t78 = sin(t90) / 0.2e1;
t91 = pkin(7) - qJ(3);
t83 = sin(t91);
t109 = t78 - t83 / 0.2e1;
t89 = pkin(6) - pkin(12);
t74 = cos(t89) / 0.2e1;
t88 = pkin(6) + pkin(12);
t82 = cos(t88);
t35 = t74 - t82 / 0.2e1;
t49 = cos(qJ(3));
t73 = sin(t88) / 0.2e1;
t81 = sin(t89);
t59 = t73 + t81 / 0.2e1;
t113 = t109 * t59 + t35 * t49;
t34 = t73 - t81 / 0.2e1;
t44 = cos(pkin(12));
t93 = sin(pkin(11));
t96 = cos(pkin(11));
t29 = -t34 * t93 + t44 * t96;
t60 = t74 + t82 / 0.2e1;
t92 = sin(pkin(12));
t53 = t60 * t93 + t92 * t96;
t112 = -t109 * t53 + t29 * t49;
t28 = t34 * t96 + t44 * t93;
t52 = -t60 * t96 + t92 * t93;
t111 = -t109 * t52 + t28 * t49;
t103 = m(3) + m(4) + t108;
t102 = -t39 * mrSges(7,1) - t43 * mrSges(6,2) - t40 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t42 - t108 * pkin(9);
t98 = cos(pkin(6));
t97 = cos(pkin(7));
t95 = sin(pkin(6));
t94 = sin(pkin(7));
t85 = cos(t91);
t84 = cos(t90);
t80 = t85 / 0.2e1;
t79 = t84 / 0.2e1;
t76 = t96 * t95;
t75 = t95 * t93;
t72 = t80 - t84 / 0.2e1;
t66 = t72 * t95;
t65 = t80 + t79;
t64 = t79 - t85 / 0.2e1;
t62 = t78 + t83 / 0.2e1;
t58 = t64 * t95;
t57 = t62 * t95;
t54 = -t59 * t94 + t97 * t98;
t51 = t53 * t94 + t75 * t97;
t50 = t52 * t94 - t76 * t97;
t47 = sin(qJ(3));
t17 = t72 * t98 + t113;
t16 = t35 * t47 - t59 * t65 - t62 * t98;
t13 = t66 * t93 + t112;
t12 = t29 * t47 + t53 * t65 - t57 * t93;
t10 = -t66 * t96 + t111;
t9 = t28 * t47 + t52 * t65 + t57 * t96;
t6 = t17 * t48 + t46 * t54;
t5 = t17 * t46 - t48 * t54;
t4 = t13 * t48 + t46 * t51;
t3 = t13 * t46 - t48 * t51;
t2 = t10 * t48 + t46 * t50;
t1 = t10 * t46 - t48 * t50;
t7 = [(-m(2) - t103) * g(3) (-g(1) * t75 + g(2) * t76 - g(3) * t98) * t103 (t102 * (-t64 * t98 + t113) + t114 * t16) * g(3) + (t102 * (t58 * t96 + t111) + t114 * t9) * g(2) + (t102 * (-t58 * t93 + t112) + t114 * t12) * g(1) (t106 * t6 + t107 * t5) * g(3) + (t1 * t107 + t106 * t2) * g(2) + (t106 * t4 + t107 * t3) * g(1), t110 * (-g(1) * t3 - g(2) * t1 - g(3) * t5) -g(1) * ((t12 * t40 - t39 * t4) * mrSges(7,1) + (-t12 * t39 - t4 * t40) * mrSges(7,2)) - g(2) * ((-t2 * t39 + t40 * t9) * mrSges(7,1) + (-t2 * t40 - t39 * t9) * mrSges(7,2)) - g(3) * ((t16 * t40 - t39 * t6) * mrSges(7,1) + (-t16 * t39 - t40 * t6) * mrSges(7,2))];
taug  = t7(:);
