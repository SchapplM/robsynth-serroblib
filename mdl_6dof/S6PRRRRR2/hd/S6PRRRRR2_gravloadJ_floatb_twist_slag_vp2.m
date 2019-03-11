% Calculate Gravitation load on the joints for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:27
% EndTime: 2019-03-09 00:42:30
% DurationCPUTime: 1.19s
% Computational Cost: add. (784->124), mult. (1241->170), div. (0->0), fcn. (1437->14), ass. (0->62)
t60 = qJ(5) + qJ(6);
t56 = sin(t60);
t58 = cos(t60);
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t158 = -t67 * mrSges(6,1) - t58 * mrSges(7,1) + t64 * mrSges(6,2) + t56 * mrSges(7,2) - mrSges(5,1);
t153 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t138 = -m(7) * pkin(5) - mrSges(6,1);
t54 = pkin(5) * t67 + pkin(4);
t61 = qJ(3) + qJ(4);
t57 = sin(t61);
t59 = cos(t61);
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t70 = -pkin(11) - pkin(10);
t156 = -m(4) * pkin(2) - t68 * mrSges(4,1) + t65 * mrSges(4,2) - mrSges(3,1) + (-m(6) * pkin(4) - m(7) * t54 + t158) * t59 + (-m(6) * pkin(10) + m(7) * t70 + t153) * t57;
t102 = cos(pkin(6));
t63 = sin(pkin(6));
t66 = sin(qJ(2));
t117 = t63 * t66;
t150 = t102 * t68 - t65 * t117;
t116 = t63 * t68;
t101 = cos(pkin(12));
t69 = cos(qJ(2));
t62 = sin(pkin(12));
t92 = t62 * t102;
t47 = t101 * t69 - t66 * t92;
t148 = t62 * t116 - t47 * t65;
t81 = t102 * t101;
t45 = t62 * t69 + t66 * t81;
t91 = t63 * t101;
t25 = -t45 * t57 - t59 * t91;
t26 = t45 * t59 - t57 * t91;
t147 = t153 * t26 + t158 * t25;
t118 = t62 * t63;
t27 = t59 * t118 - t47 * t57;
t28 = t57 * t118 + t47 * t59;
t146 = t153 * t28 + t158 * t27;
t40 = t102 * t59 - t57 * t117;
t41 = t102 * t57 + t59 * t117;
t145 = t153 * t41 + t158 * t40;
t132 = -m(4) * pkin(8) - t56 * mrSges(7,1) - t67 * mrSges(6,2) - t58 * mrSges(7,2) + t138 * t64 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t137 = m(5) + m(6) + m(7);
t44 = t62 * t66 - t69 * t81;
t130 = (-t26 * t56 + t44 * t58) * mrSges(7,1) + (-t26 * t58 - t44 * t56) * mrSges(7,2);
t46 = t101 * t66 + t69 * t92;
t129 = (-t28 * t56 + t46 * t58) * mrSges(7,1) + (-t28 * t58 - t46 * t56) * mrSges(7,2);
t124 = t25 * t54 - t26 * t70;
t123 = t27 * t54 - t28 * t70;
t115 = t63 * t69;
t122 = (-t58 * t115 - t41 * t56) * mrSges(7,1) + (t56 * t115 - t41 * t58) * mrSges(7,2);
t108 = t40 * t54 - t41 * t70;
t95 = t25 * pkin(4) + t26 * pkin(10);
t94 = t27 * pkin(4) + pkin(10) * t28;
t93 = t40 * pkin(4) + pkin(10) * t41;
t89 = t148 * pkin(3);
t82 = t150 * pkin(3);
t76 = -t45 * t65 - t68 * t91;
t74 = t76 * pkin(3);
t71 = -pkin(9) - pkin(8);
t55 = pkin(3) * t68 + pkin(2);
t1 = [(-m(2) - m(3) - m(4) - t137) * g(3) (-t137 * (-t44 * t55 - t45 * t71) + t132 * t45 - t156 * t44) * g(2) + (-t137 * (-t46 * t55 - t47 * t71) + t132 * t47 - t156 * t46) * g(1) + (-t137 * t55 * t115 + (t156 * t69 + (t137 * t71 + t132) * t66) * t63) * g(3) (-t150 * mrSges(4,1) - (-t102 * t65 - t66 * t116) * mrSges(4,2) - m(5) * t82 - m(6) * (t82 + t93) - m(7) * (t82 + t108) + t145) * g(3) + (-t76 * mrSges(4,1) - (-t45 * t68 + t65 * t91) * mrSges(4,2) - m(5) * t74 - m(6) * (t74 + t95) - m(7) * (t74 + t124) + t147) * g(2) + (-t148 * mrSges(4,1) - (-t65 * t118 - t47 * t68) * mrSges(4,2) - m(5) * t89 - m(6) * (t89 + t94) - m(7) * (t89 + t123) + t146) * g(1) (-m(6) * t93 - m(7) * t108 + t145) * g(3) + (-m(6) * t95 - m(7) * t124 + t147) * g(2) + (-m(6) * t94 - m(7) * t123 + t146) * g(1) (-(t64 * t115 - t41 * t67) * mrSges(6,2) - t122 + t138 * (-t67 * t115 - t41 * t64)) * g(3) + (-(-t26 * t67 - t44 * t64) * mrSges(6,2) - t130 + t138 * (-t26 * t64 + t44 * t67)) * g(2) + (-(-t28 * t67 - t46 * t64) * mrSges(6,2) - t129 + t138 * (-t28 * t64 + t46 * t67)) * g(1), -g(1) * t129 - g(2) * t130 - g(3) * t122];
taug  = t1(:);
