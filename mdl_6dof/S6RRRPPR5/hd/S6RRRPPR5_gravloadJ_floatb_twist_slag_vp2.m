% Calculate Gravitation load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:21
% EndTime: 2019-03-09 15:38:26
% DurationCPUTime: 1.57s
% Computational Cost: add. (790->122), mult. (1331->168), div. (0->0), fcn. (1527->14), ass. (0->55)
t45 = pkin(12) + qJ(6);
t40 = sin(t45);
t42 = cos(t45);
t47 = sin(pkin(12));
t49 = cos(pkin(12));
t59 = -mrSges(5,1) - m(6) * pkin(4) - t49 * mrSges(6,1) + t47 * mrSges(6,2) - m(7) * (pkin(5) * t49 + pkin(4));
t132 = -t42 * mrSges(7,1) + t40 * mrSges(7,2) + t59;
t112 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t111 = -m(4) * pkin(9) - t49 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t47 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t131 = t40 * mrSges(7,1) + t42 * mrSges(7,2) - t111;
t46 = qJ(3) + pkin(11);
t41 = sin(t46);
t43 = cos(t46);
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t110 = m(4) * pkin(2) + t55 * mrSges(4,1) - t52 * mrSges(4,2) - t112 * t41 - t132 * t43 + mrSges(3,1);
t121 = m(6) + m(7);
t115 = m(5) + t121;
t91 = cos(pkin(6));
t48 = sin(pkin(6));
t53 = sin(qJ(2));
t99 = t48 * t53;
t123 = -t52 * t99 + t91 * t55;
t104 = cos(qJ(1));
t56 = cos(qJ(2));
t54 = sin(qJ(1));
t80 = t54 * t91;
t24 = t104 * t56 - t53 * t80;
t97 = t48 * t55;
t9 = -t24 * t52 + t54 * t97;
t75 = t91 * t104;
t22 = t53 * t75 + t54 * t56;
t85 = t48 * t104;
t66 = t22 * t52 + t55 * t85;
t98 = t48 * t54;
t96 = t48 * t56;
t92 = t104 * pkin(1) + pkin(8) * t98;
t90 = t52 * t98;
t83 = -pkin(1) * t54 + pkin(8) * t85;
t4 = t22 * t43 - t41 * t85;
t32 = t52 * t85;
t81 = -t22 * t55 + t32;
t3 = t22 * t41 + t43 * t85;
t50 = -qJ(4) - pkin(9);
t39 = pkin(3) * t55 + pkin(2);
t23 = t104 * t53 + t56 * t80;
t21 = t53 * t54 - t56 * t75;
t16 = t91 * t41 + t43 * t99;
t15 = t41 * t99 - t91 * t43;
t10 = t24 * t55 + t90;
t8 = t24 * t43 + t41 * t98;
t7 = t24 * t41 - t43 * t98;
t2 = t23 * t40 + t42 * t8;
t1 = t23 * t42 - t40 * t8;
t5 = [(-t104 * mrSges(2,1) - m(3) * t92 - t24 * mrSges(3,1) - m(4) * (pkin(2) * t24 + t92) - t10 * mrSges(4,1) - t9 * mrSges(4,2) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t48 + mrSges(2,2)) * t54 + t59 * t8 + t112 * t7 + t111 * t23 - t115 * (pkin(3) * t90 - t23 * t50 + t24 * t39 + t92)) * g(2) + (t54 * mrSges(2,1) + t104 * mrSges(2,2) - m(3) * t83 + t22 * mrSges(3,1) - mrSges(3,3) * t85 - m(4) * (-pkin(2) * t22 + t83) - t81 * mrSges(4,1) - t66 * mrSges(4,2) - t112 * t3 - t132 * t4 + t131 * t21 + t115 * (-pkin(3) * t32 - t21 * t50 + t22 * t39 - t83)) * g(1) (-t115 * (-t21 * t39 - t22 * t50) - t131 * t22 + t110 * t21) * g(2) + (-t115 * (-t23 * t39 - t24 * t50) - t131 * t24 + t110 * t23) * g(1) + (-t115 * t39 * t96 + (-t110 * t56 + (t115 * t50 - t131) * t53) * t48) * g(3) (-t123 * mrSges(4,1) - (-t91 * t52 - t53 * t97) * mrSges(4,2) + t112 * t16 - t132 * t15) * g(3) + (t66 * mrSges(4,1) - t81 * mrSges(4,2) + t112 * t4 - t132 * t3) * g(2) + (-mrSges(4,1) * t9 + mrSges(4,2) * t10 + t112 * t8 - t132 * t7) * g(1) + (-g(1) * t9 + g(2) * t66 - g(3) * t123) * t115 * pkin(3), t115 * (-g(1) * t23 - g(2) * t21 + g(3) * t96) t121 * (-g(1) * t7 - g(2) * t3 - g(3) * t15) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t42 - t4 * t40) * mrSges(7,1) + (-t21 * t40 - t4 * t42) * mrSges(7,2)) - g(3) * ((-t16 * t40 - t42 * t96) * mrSges(7,1) + (-t16 * t42 + t40 * t96) * mrSges(7,2))];
taug  = t5(:);
