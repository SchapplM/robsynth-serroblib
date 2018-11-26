% Calculate Gravitation load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:56:30
% EndTime: 2018-11-23 14:56:30
% DurationCPUTime: 0.82s
% Computational Cost: add. (1090->112), mult. (928->147), div. (0->0), fcn. (850->20), ass. (0->76)
t110 = m(6) + m(7);
t116 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t115 = t56 * mrSges(7,1) + t59 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t51 = pkin(6) - qJ(2);
t40 = sin(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t44 = sin(t50);
t112 = t44 / 0.2e1 + t40;
t111 = m(7) * pkin(9) + pkin(4) * t110 + t116;
t49 = qJ(2) + pkin(11);
t43 = sin(t49);
t52 = sin(pkin(10));
t54 = cos(pkin(10));
t80 = pkin(6) - t49;
t70 = cos(t80) / 0.2e1;
t79 = pkin(6) + t49;
t73 = cos(t79);
t62 = t73 / 0.2e1 + t70;
t13 = t43 * t52 - t54 * t62;
t57 = sin(qJ(4));
t84 = qJ(5) * t57;
t60 = cos(qJ(4));
t97 = t13 * t60;
t109 = -pkin(4) * t97 - t13 * t84;
t16 = t54 * t43 + t52 * t62;
t96 = t16 * t60;
t108 = -pkin(4) * t96 - t16 * t84;
t39 = sin(t79);
t101 = t39 / 0.2e1;
t72 = sin(t80);
t69 = t72 / 0.2e1;
t28 = t69 + t101;
t95 = t28 * t60;
t107 = pkin(4) * t95 + t28 * t84;
t41 = cos(t50) / 0.2e1;
t48 = cos(t51);
t32 = t48 / 0.2e1 + t41;
t104 = -t110 * qJ(5) - t115;
t103 = t115 * t57 + t116 * t60 + mrSges(4,1);
t102 = -m(7) * (pkin(5) + pkin(8)) - t59 * mrSges(7,1) + t56 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t46 = cos(t49);
t94 = t52 * t46;
t58 = sin(qJ(2));
t93 = t52 * t58;
t53 = sin(pkin(6));
t92 = t53 * t57;
t91 = t53 * t60;
t90 = t54 * t46;
t89 = t54 * t58;
t85 = t112 * pkin(2);
t83 = t101 - t72 / 0.2e1;
t82 = t28 * pkin(3) + t85;
t81 = m(4) + m(5) + t110;
t30 = t32 * pkin(2);
t78 = -pkin(2) * t93 + t54 * t30;
t75 = -t13 * pkin(3) + t78;
t29 = t70 - t73 / 0.2e1;
t74 = t29 * pkin(8) + t82;
t71 = -pkin(2) * t89 - t30 * t52;
t68 = -t16 * pkin(3) + t71;
t64 = -t39 / 0.2e1 + t69;
t15 = -t54 * t64 + t94;
t65 = t15 * pkin(8) + t75;
t18 = t52 * t64 + t90;
t63 = pkin(8) * t18 + t68;
t61 = cos(qJ(2));
t55 = cos(pkin(6));
t31 = t40 - t44 / 0.2e1;
t20 = t29 * t57 - t55 * t60;
t17 = -t52 * t83 + t90;
t14 = t54 * t83 + t94;
t5 = t17 * t57 - t52 * t91;
t3 = t14 * t57 + t54 * t91;
t1 = [(-m(2) - m(3) - t81) * g(3) (-t112 * mrSges(3,1) - (t41 - t48 / 0.2e1) * mrSges(3,2) - m(4) * t85 - m(5) * t74 - m(6) * (t74 + t107) - m(7) * (pkin(9) * t95 + t107 + t82) + t102 * t29 - t103 * t28) * g(3) + (-(t54 * t32 - t93) * mrSges(3,1) - (t54 * t31 - t52 * t61) * mrSges(3,2) - m(4) * t78 - m(5) * t65 - m(6) * (t65 + t109) - m(7) * (-pkin(9) * t97 + t109 + t75) + t102 * t15 + t103 * t13) * g(2) + (-(-t32 * t52 - t89) * mrSges(3,1) - (-t52 * t31 - t54 * t61) * mrSges(3,2) - m(4) * t71 - m(5) * t63 - m(6) * (t63 + t108) - m(7) * (-pkin(9) * t96 + t108 + t68) + t102 * t18 + t103 * t16) * g(1) (-g(3) * t55 + (-g(1) * t52 + g(2) * t54) * t53) * t81 (t104 * (t29 * t60 + t55 * t57) + t111 * t20) * g(3) + (t104 * (t14 * t60 - t54 * t92) + t111 * t3) * g(2) + (t104 * (t17 * t60 + t52 * t92) + t111 * t5) * g(1), t110 * (-g(1) * t5 - g(2) * t3 - g(3) * t20) -g(1) * ((-t16 * t56 + t5 * t59) * mrSges(7,1) + (-t16 * t59 - t5 * t56) * mrSges(7,2)) - g(2) * ((-t13 * t56 + t3 * t59) * mrSges(7,1) + (-t13 * t59 - t3 * t56) * mrSges(7,2)) - g(3) * ((t20 * t59 + t28 * t56) * mrSges(7,1) + (-t20 * t56 + t28 * t59) * mrSges(7,2))];
taug  = t1(:);
