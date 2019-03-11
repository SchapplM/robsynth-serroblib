% Calculate Gravitation load on the joints for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:05
% EndTime: 2019-03-08 21:51:08
% DurationCPUTime: 0.92s
% Computational Cost: add. (651->117), mult. (892->162), div. (0->0), fcn. (997->14), ass. (0->56)
t115 = mrSges(6,2) - mrSges(7,3);
t55 = sin(qJ(6));
t58 = cos(qJ(6));
t114 = -mrSges(7,1) * t58 + mrSges(7,2) * t55 - mrSges(6,1);
t54 = -qJ(4) - pkin(8);
t101 = -m(4) * pkin(8) + m(5) * t54 - t55 * mrSges(7,1) - t58 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t51 = qJ(3) + pkin(12);
t48 = qJ(5) + t51;
t43 = sin(t48);
t44 = cos(t48);
t46 = sin(t51);
t47 = cos(t51);
t59 = cos(qJ(3));
t49 = t59 * pkin(3);
t56 = sin(qJ(3));
t112 = -mrSges(3,1) - m(5) * (t49 + pkin(2)) - t47 * mrSges(5,1) + t46 * mrSges(5,2) - m(4) * pkin(2) - t59 * mrSges(4,1) + t56 * mrSges(4,2) + (-m(7) * pkin(5) + t114) * t44 + (-m(7) * pkin(10) + t115) * t43;
t111 = m(6) + m(7);
t85 = cos(pkin(6));
t53 = sin(pkin(6));
t57 = sin(qJ(2));
t95 = t53 * t57;
t24 = -t43 * t95 + t44 * t85;
t25 = t43 * t85 + t44 * t95;
t109 = t114 * t24 + t115 * t25;
t60 = cos(qJ(2));
t52 = sin(pkin(11));
t77 = t52 * t85;
t84 = cos(pkin(11));
t34 = -t57 * t77 + t60 * t84;
t96 = t52 * t53;
t13 = -t34 * t43 + t44 * t96;
t14 = t34 * t44 + t43 * t96;
t108 = t114 * t13 + t115 * t14;
t71 = t85 * t84;
t32 = t52 * t60 + t57 * t71;
t76 = t53 * t84;
t11 = -t32 * t43 - t44 * t76;
t12 = t32 * t44 - t43 * t76;
t107 = t11 * t114 + t115 * t12;
t106 = -m(5) * pkin(3) - mrSges(4,1);
t105 = m(5) + t111;
t94 = t53 * t59;
t93 = t53 * t60;
t38 = -pkin(3) * t56 - pkin(4) * t46;
t39 = pkin(4) * t47 + t49;
t88 = t34 * t38 + t39 * t96;
t86 = t38 * t95 + t39 * t85;
t81 = pkin(5) * t11 + pkin(10) * t12;
t79 = pkin(5) * t13 + pkin(10) * t14;
t78 = pkin(5) * t24 + pkin(10) * t25;
t68 = t32 * t38 - t39 * t76;
t50 = -pkin(9) + t54;
t37 = pkin(2) + t39;
t33 = t57 * t84 + t60 * t77;
t31 = t52 * t57 - t60 * t71;
t1 = [(-m(2) - m(3) - m(4) - t105) * g(3) (-t111 * (-t31 * t37 - t32 * t50) + t101 * t32 - t112 * t31) * g(2) + (-t111 * (-t33 * t37 - t34 * t50) + t101 * t34 - t112 * t33) * g(1) + (-t111 * t37 * t93 + (t112 * t60 + (t111 * t50 + t101) * t57) * t53) * g(3) (-(-t56 * t85 - t57 * t94) * mrSges(4,2) - (-t46 * t95 + t47 * t85) * mrSges(5,1) - (-t46 * t85 - t47 * t95) * mrSges(5,2) - m(6) * t86 - m(7) * (t78 + t86) + t106 * (-t56 * t95 + t59 * t85) + t109) * g(3) + (-(-t32 * t59 + t56 * t76) * mrSges(4,2) - (-t32 * t46 - t47 * t76) * mrSges(5,1) - (-t32 * t47 + t46 * t76) * mrSges(5,2) - m(6) * t68 - m(7) * (t68 + t81) + t106 * (-t32 * t56 - t59 * t76) + t107) * g(2) + (-(-t34 * t59 - t56 * t96) * mrSges(4,2) - (-t34 * t46 + t47 * t96) * mrSges(5,1) - (-t34 * t47 - t46 * t96) * mrSges(5,2) - m(6) * t88 - m(7) * (t79 + t88) + t106 * (-t34 * t56 + t52 * t94) + t108) * g(1), t105 * (-g(1) * t33 - g(2) * t31 + g(3) * t93) (-m(7) * t78 + t109) * g(3) + (-m(7) * t81 + t107) * g(2) + (-m(7) * t79 + t108) * g(1), -g(1) * ((-t14 * t55 + t33 * t58) * mrSges(7,1) + (-t14 * t58 - t33 * t55) * mrSges(7,2)) - g(2) * ((-t12 * t55 + t31 * t58) * mrSges(7,1) + (-t12 * t58 - t31 * t55) * mrSges(7,2)) - g(3) * ((-t25 * t55 - t58 * t93) * mrSges(7,1) + (-t25 * t58 + t55 * t93) * mrSges(7,2))];
taug  = t1(:);
