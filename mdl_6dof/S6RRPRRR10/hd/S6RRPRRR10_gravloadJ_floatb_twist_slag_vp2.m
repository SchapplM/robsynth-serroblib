% Calculate Gravitation load on the joints for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:21
% EndTime: 2019-03-09 14:18:24
% DurationCPUTime: 1.30s
% Computational Cost: add. (834->129), mult. (1364->172), div. (0->0), fcn. (1583->14), ass. (0->57)
t57 = cos(qJ(5));
t42 = pkin(5) * t57 + pkin(4);
t49 = qJ(5) + qJ(6);
t45 = sin(t49);
t46 = cos(t49);
t54 = sin(qJ(5));
t102 = m(6) * pkin(4) + m(7) * t42 + t57 * mrSges(6,1) + t46 * mrSges(7,1) - t54 * mrSges(6,2) - t45 * mrSges(7,2) + mrSges(5,1);
t63 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t50 = sin(pkin(12));
t52 = cos(pkin(12));
t108 = -m(4) * pkin(2) - t52 * mrSges(4,1) + t50 * mrSges(4,2) - mrSges(3,1);
t48 = pkin(12) + qJ(4);
t43 = sin(t48);
t44 = cos(t48);
t101 = t102 * t44 - t63 * t43 - t108;
t69 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t115 = t45 * mrSges(7,1) + t57 * mrSges(6,2) + t46 * mrSges(7,2) + t69;
t111 = m(7) * pkin(5);
t84 = t54 * t111;
t100 = -t54 * mrSges(6,1) - t115 - t84;
t107 = mrSges(6,1) + t111;
t103 = m(5) + m(6) + m(7);
t55 = sin(qJ(2));
t56 = sin(qJ(1));
t58 = cos(qJ(2));
t85 = cos(pkin(6));
t93 = cos(qJ(1));
t74 = t85 * t93;
t28 = t55 * t74 + t56 * t58;
t51 = sin(pkin(6));
t81 = t51 * t93;
t12 = t28 * t44 - t43 * t81;
t27 = t55 * t56 - t58 * t74;
t98 = (-t12 * t45 + t27 * t46) * mrSges(7,1) + (-t12 * t46 - t27 * t45) * mrSges(7,2);
t77 = t56 * t85;
t30 = -t55 * t77 + t58 * t93;
t91 = t51 * t56;
t16 = t30 * t44 + t43 * t91;
t29 = t55 * t93 + t58 * t77;
t5 = -t16 * t45 + t29 * t46;
t6 = t16 * t46 + t29 * t45;
t97 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t92 = t51 * t55;
t22 = t43 * t85 + t44 * t92;
t90 = t51 * t58;
t94 = (-t22 * t45 - t46 * t90) * mrSges(7,1) + (-t22 * t46 + t45 * t90) * mrSges(7,2);
t86 = t93 * pkin(1) + pkin(8) * t91;
t79 = -t56 * pkin(1) + pkin(8) * t81;
t11 = -t28 * t43 - t44 * t81;
t76 = t50 * t81;
t41 = pkin(3) * t52 + pkin(2);
t53 = -pkin(9) - qJ(3);
t75 = t50 * pkin(3) * t91 - t29 * t53 + t30 * t41 + t86;
t7 = -t16 * t54 + t29 * t57;
t15 = t30 * t43 - t44 * t91;
t8 = t16 * t57 + t29 * t54;
t1 = [(-t93 * mrSges(2,1) - m(5) * t75 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t75) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t42 + t75) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t108 * t30 + (mrSges(2,2) + (-mrSges(4,1) * t50 - mrSges(4,2) * t52 - mrSges(3,3)) * t51) * t56 + (-t69 - t84) * t29 + t63 * t15 + (-m(3) - m(4)) * t86) * g(2) + (t56 * mrSges(2,1) + t93 * mrSges(2,2) - m(3) * t79 + t28 * mrSges(3,1) - mrSges(3,3) * t81 - m(4) * (-pkin(2) * t28 + t79) - (-t28 * t52 + t76) * mrSges(4,1) - (t28 * t50 + t52 * t81) * mrSges(4,2) + t102 * t12 + (t107 * t54 + t115) * t27 + t63 * t11 + t103 * (-pkin(3) * t76 - t27 * t53 + t28 * t41 - t79)) * g(1) (-t103 * (-t27 * t41 - t28 * t53) + t100 * t28 + t101 * t27) * g(2) + (-t103 * (-t29 * t41 - t30 * t53) + t100 * t30 + t101 * t29) * g(1) + (-t103 * t41 * t90 + (-t101 * t58 + (t103 * t53 + t100) * t55) * t51) * g(3) (-g(1) * t29 - g(2) * t27 + g(3) * t90) * (m(4) + t103) (t63 * t22 - t102 * (-t43 * t92 + t44 * t85)) * g(3) + (-t102 * t11 + t63 * t12) * g(2) + (t102 * t15 + t63 * t16) * g(1) (-(-t22 * t57 + t54 * t90) * mrSges(6,2) - t94 - t107 * (-t22 * t54 - t57 * t90)) * g(3) + (-(-t12 * t57 - t27 * t54) * mrSges(6,2) - t98 - t107 * (-t12 * t54 + t27 * t57)) * g(2) + (mrSges(6,2) * t8 - t107 * t7 - t97) * g(1), -g(1) * t97 - g(2) * t98 - g(3) * t94];
taug  = t1(:);
