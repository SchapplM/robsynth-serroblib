% Calculate Gravitation load on the joints for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2018-11-23 14:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:55:04
% EndTime: 2018-11-23 14:55:05
% DurationCPUTime: 0.90s
% Computational Cost: add. (1058->110), mult. (842->153), div. (0->0), fcn. (754->22), ass. (0->65)
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t111 = m(7) * pkin(5) + t59 * mrSges(7,1) - t56 * mrSges(7,2) + mrSges(6,1);
t103 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t106 = m(6) + m(7);
t51 = sin(pkin(10));
t57 = sin(qJ(4));
t48 = qJ(2) + pkin(11);
t82 = pkin(6) + t48;
t68 = sin(t82) / 0.2e1;
t83 = pkin(6) - t48;
t73 = sin(t83);
t22 = t68 - t73 / 0.2e1;
t44 = cos(t48);
t53 = cos(pkin(10));
t75 = -t22 * t51 + t44 * t53;
t52 = sin(pkin(6));
t60 = cos(qJ(4));
t92 = t52 * t60;
t110 = t51 * t92 - t75 * t57;
t69 = cos(t83) / 0.2e1;
t74 = cos(t82);
t24 = t69 - t74 / 0.2e1;
t54 = cos(pkin(6));
t109 = -t24 * t57 + t54 * t60;
t50 = pkin(6) - qJ(2);
t36 = sin(t50) / 0.2e1;
t49 = pkin(6) + qJ(2);
t41 = sin(t49);
t108 = t41 / 0.2e1 + t36;
t107 = -m(4) - m(5);
t76 = t53 * t22 + t44 * t51;
t70 = -t53 * t92 - t76 * t57;
t37 = cos(t49) / 0.2e1;
t46 = cos(t50);
t27 = t46 / 0.2e1 + t37;
t47 = qJ(4) + pkin(12);
t39 = sin(t47);
t43 = cos(t47);
t104 = m(5) * pkin(3) + t60 * mrSges(5,1) - t57 * mrSges(5,2) - t103 * t39 + t111 * t43 + mrSges(4,1);
t102 = -m(5) * pkin(8) - t56 * mrSges(7,1) - t59 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t96 = t51 * t52;
t58 = sin(qJ(2));
t95 = t51 * t58;
t94 = t52 * t53;
t93 = t52 * t57;
t91 = t53 * t58;
t89 = t108 * pkin(2);
t86 = t106 - t107;
t25 = t27 * pkin(2);
t81 = -pkin(2) * t95 + t53 * t25;
t72 = -pkin(2) * t91 - t25 * t51;
t63 = t74 / 0.2e1 + t69;
t61 = cos(qJ(2));
t55 = -qJ(5) - pkin(8);
t40 = sin(t48);
t38 = pkin(4) * t60 + pkin(3);
t26 = t36 - t41 / 0.2e1;
t23 = t73 / 0.2e1 + t68;
t14 = t53 * t40 + t51 * t63;
t11 = t40 * t51 - t53 * t63;
t10 = t24 * t43 + t39 * t54;
t4 = t39 * t96 + t43 * t75;
t2 = -t39 * t94 + t43 * t76;
t1 = [(-m(2) - m(3) - t86) * g(3) (-t108 * mrSges(3,1) - (t37 - t46 / 0.2e1) * mrSges(3,2) + t107 * t89 - t106 * (t23 * t38 - t24 * t55 + t89) + t102 * t24 - t104 * t23) * g(3) + (-(t27 * t53 - t95) * mrSges(3,1) - (t53 * t26 - t51 * t61) * mrSges(3,2) + t107 * t81 - t106 * (-t11 * t38 - t55 * t76 + t81) + t102 * t76 + t104 * t11) * g(2) + (-(-t27 * t51 - t91) * mrSges(3,1) - (-t51 * t26 - t53 * t61) * mrSges(3,2) + t107 * t72 - t106 * (-t14 * t38 - t55 * t75 + t72) + t102 * t75 + t104 * t14) * g(1) (-g(3) * t54 + (-g(1) * t51 + g(2) * t53) * t52) * t86 (-t109 * mrSges(5,1) - (-t24 * t60 - t54 * t57) * mrSges(5,2) - t111 * (-t24 * t39 + t43 * t54) + t103 * t10) * g(3) + (-t70 * mrSges(5,1) - (t53 * t93 - t60 * t76) * mrSges(5,2) - t111 * (-t39 * t76 - t43 * t94) + t103 * t2) * g(2) + (-t110 * mrSges(5,1) - (-t51 * t93 - t60 * t75) * mrSges(5,2) + t103 * t4 - t111 * (-t39 * t75 + t43 * t96)) * g(1) + (-g(1) * t110 - g(2) * t70 - g(3) * t109) * t106 * pkin(4), t106 * (-g(1) * t14 - g(2) * t11 + g(3) * t23) -g(1) * ((t14 * t59 - t4 * t56) * mrSges(7,1) + (-t14 * t56 - t4 * t59) * mrSges(7,2)) - g(2) * ((t11 * t59 - t2 * t56) * mrSges(7,1) + (-t11 * t56 - t2 * t59) * mrSges(7,2)) - g(3) * ((-t10 * t56 - t23 * t59) * mrSges(7,1) + (-t10 * t59 + t23 * t56) * mrSges(7,2))];
taug  = t1(:);
