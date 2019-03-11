% Calculate Gravitation load on the joints for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:23
% EndTime: 2019-03-09 23:12:28
% DurationCPUTime: 1.68s
% Computational Cost: add. (876->159), mult. (1684->204), div. (0->0), fcn. (1987->14), ass. (0->74)
t51 = qJ(4) + pkin(12);
t46 = cos(t51);
t57 = cos(qJ(4));
t48 = t57 * pkin(4);
t31 = pkin(5) * t46 + t48;
t29 = pkin(3) + t31;
t47 = qJ(6) + t51;
t42 = sin(t47);
t43 = cos(t47);
t44 = t48 + pkin(3);
t45 = sin(t51);
t54 = sin(qJ(4));
t114 = m(5) * pkin(3) + m(6) * t44 + m(7) * t29 + t57 * mrSges(5,1) + t46 * mrSges(6,1) + t43 * mrSges(7,1) - t54 * mrSges(5,2) - t45 * mrSges(6,2) - t42 * mrSges(7,2) + mrSges(4,1);
t53 = -qJ(5) - pkin(10);
t132 = mrSges(4,2) + m(7) * (-pkin(11) + t53) - mrSges(7,3) + m(6) * t53 - mrSges(6,3) - m(5) * pkin(10) - mrSges(5,3);
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t138 = -t114 * t58 + t132 * t55 - mrSges(3,1);
t136 = -m(4) - m(5);
t104 = pkin(4) * t54;
t30 = pkin(5) * t45 + t104;
t134 = m(7) * (pkin(9) + t30) + m(6) * (pkin(9) + t104);
t123 = m(6) + m(7);
t115 = t123 - t136;
t131 = pkin(2) * t115 - t138;
t101 = sin(qJ(1));
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t102 = cos(qJ(1));
t89 = cos(pkin(6));
t80 = t89 * t102;
t26 = t101 * t59 + t56 * t80;
t52 = sin(pkin(6));
t85 = t52 * t102;
t14 = t26 * t58 - t55 * t85;
t25 = t101 * t56 - t59 * t80;
t130 = t14 * t42 - t25 * t43;
t129 = t14 * t45 - t25 * t46;
t128 = t14 * t46 + t25 * t45;
t127 = t14 * t54 - t25 * t57;
t126 = t14 * t57 + t25 * t54;
t125 = -t14 * t43 - t25 * t42;
t119 = mrSges(4,3) - mrSges(3,2);
t124 = -t54 * mrSges(5,1) - t45 * mrSges(6,1) - t42 * mrSges(7,1) - t57 * mrSges(5,2) - t46 * mrSges(6,2) - t43 * mrSges(7,2) - t119;
t118 = -m(6) * pkin(4) - mrSges(5,1);
t110 = t119 + t134;
t109 = t136 * pkin(9) + t124 - t134;
t106 = -t130 * mrSges(7,1) + t125 * mrSges(7,2);
t79 = t89 * t101;
t28 = t102 * t59 - t56 * t79;
t84 = t52 * t101;
t18 = t28 * t58 + t55 * t84;
t27 = t102 * t56 + t59 * t79;
t5 = -t18 * t42 + t27 * t43;
t6 = t18 * t43 + t27 * t42;
t105 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t94 = t52 * t56;
t93 = t52 * t59;
t24 = t55 * t89 + t58 * t94;
t92 = (-t24 * t42 - t43 * t93) * mrSges(7,1) + (-t24 * t43 + t42 * t93) * mrSges(7,2);
t90 = t102 * pkin(1) + pkin(8) * t84;
t88 = t28 * pkin(2) + t90;
t81 = -pkin(1) * t101 + pkin(8) * t85;
t9 = -t18 * t54 + t27 * t57;
t73 = pkin(9) * t27 + t88;
t72 = -t26 * pkin(2) + t81;
t13 = t26 * t55 + t58 * t85;
t65 = -t25 * pkin(9) + t72;
t23 = t55 * t94 - t58 * t89;
t17 = t28 * t55 - t58 * t84;
t10 = t18 * t57 + t27 * t54;
t8 = t18 * t46 + t27 * t45;
t7 = -t18 * t45 + t27 * t46;
t1 = [(-t102 * mrSges(2,1) + t101 * mrSges(2,2) - m(3) * t90 - t28 * mrSges(3,1) - mrSges(3,3) * t84 - m(4) * t73 - t18 * mrSges(4,1) - m(5) * (pkin(3) * t18 + t73) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * (t18 * t44 + t88) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t18 * t29 + t88) - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t110 * t27 + t132 * t17) * g(2) + (t101 * mrSges(2,1) + t102 * mrSges(2,2) - m(3) * t81 + t26 * mrSges(3,1) - mrSges(3,3) * t85 - m(4) * t65 + t14 * mrSges(4,1) - m(5) * (-pkin(3) * t14 + t65) + t126 * mrSges(5,1) - t127 * mrSges(5,2) - m(6) * (-t14 * t44 + t72) + t128 * mrSges(6,1) - t129 * mrSges(6,2) - m(7) * (-t14 * t29 + t72) - t125 * mrSges(7,1) - t130 * mrSges(7,2) + t110 * t25 - t132 * t13) * g(1) (-t115 * (pkin(2) * t93 + pkin(9) * t94) + (t138 * t59 + (-m(6) * t104 - m(7) * t30 + t124) * t56) * t52) * g(3) + (t109 * t26 + t131 * t25) * g(2) + (t109 * t28 + t131 * t27) * g(1) (t114 * t23 + t132 * t24) * g(3) + (t114 * t13 + t132 * t14) * g(2) + (t114 * t17 + t132 * t18) * g(1) (-(-t24 * t57 + t54 * t93) * mrSges(5,2) - (-t24 * t45 - t46 * t93) * mrSges(6,1) - (-t24 * t46 + t45 * t93) * mrSges(6,2) - m(7) * (-t24 * t30 - t31 * t93) - t92 + t118 * (-t24 * t54 - t57 * t93)) * g(3) + (t126 * mrSges(5,2) + t129 * mrSges(6,1) + t128 * mrSges(6,2) - m(7) * (-t14 * t30 + t25 * t31) - t106 - t118 * t127) * g(2) + (t10 * mrSges(5,2) - t7 * mrSges(6,1) + t8 * mrSges(6,2) - m(7) * (-t18 * t30 + t27 * t31) - t105 + t118 * t9) * g(1), t123 * (-g(1) * t17 - g(2) * t13 - g(3) * t23) -g(1) * t105 - g(2) * t106 - g(3) * t92];
taug  = t1(:);
