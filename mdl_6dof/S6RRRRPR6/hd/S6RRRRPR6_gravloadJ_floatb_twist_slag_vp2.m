% Calculate Gravitation load on the joints for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:15:21
% EndTime: 2018-11-23 18:15:21
% DurationCPUTime: 0.96s
% Computational Cost: add. (651->149), mult. (672->177), div. (0->0), fcn. (626->12), ass. (0->77)
t119 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t54 = qJ(3) + qJ(4);
t47 = cos(t54);
t42 = pkin(4) * t47;
t58 = cos(qJ(3));
t50 = t58 * pkin(3);
t36 = t42 + t50;
t45 = pkin(11) + t54;
t41 = cos(t45);
t37 = pkin(5) * t41;
t27 = t37 + t36;
t21 = pkin(2) + t27;
t34 = pkin(2) + t36;
t43 = qJ(6) + t45;
t38 = sin(t43);
t39 = cos(t43);
t40 = sin(t45);
t44 = t50 + pkin(2);
t46 = sin(t54);
t55 = sin(qJ(3));
t118 = -m(4) * pkin(2) - m(5) * t44 - m(6) * t34 - m(7) * t21 - t58 * mrSges(4,1) - t47 * mrSges(5,1) - t41 * mrSges(6,1) - t39 * mrSges(7,1) + t55 * mrSges(4,2) + t46 * mrSges(5,2) + t40 * mrSges(6,2) + t38 * mrSges(7,2);
t61 = -pkin(9) - pkin(8);
t53 = -qJ(5) + t61;
t48 = -pkin(10) + t53;
t117 = -m(4) * pkin(8) + m(5) * t61 + m(6) * t53 + m(7) * t48 - t119;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t116 = g(1) * t60 + g(2) * t57;
t105 = m(5) * pkin(3);
t115 = m(6) + m(7);
t114 = mrSges(4,1) + t105;
t73 = -mrSges(7,1) * t38 - mrSges(7,2) * t39;
t113 = mrSges(5,1) * t46 + mrSges(6,1) * t40 + mrSges(5,2) * t47 + mrSges(6,2) * t41 - t73;
t59 = cos(qJ(2));
t88 = t59 * t60;
t7 = -t38 * t88 + t39 * t57;
t8 = t38 * t57 + t39 * t88;
t102 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t15 = -t40 * t88 + t41 * t57;
t16 = t40 * t57 + t41 * t88;
t24 = -t46 * t88 + t47 * t57;
t25 = t46 * t57 + t47 * t88;
t112 = -t24 * mrSges(5,1) - t15 * mrSges(6,1) + t25 * mrSges(5,2) + t16 * mrSges(6,2) - t102;
t89 = t57 * t59;
t5 = t38 * t89 + t39 * t60;
t6 = t38 * t60 - t39 * t89;
t103 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t13 = t40 * t89 + t41 * t60;
t14 = t40 * t60 - t41 * t89;
t22 = t46 * t89 + t47 * t60;
t23 = t46 * t60 - t47 * t89;
t111 = t22 * mrSges(5,1) + t13 * mrSges(6,1) - t23 * mrSges(5,2) - t14 * mrSges(6,2) - t103;
t110 = -m(3) - m(4) - m(5) - t115;
t101 = pkin(3) * t55;
t100 = pkin(4) * t46;
t32 = -pkin(5) * t40 - t100;
t26 = -t32 + t101;
t35 = t100 + t101;
t108 = -m(6) * t35 - m(7) * t26;
t56 = sin(qJ(2));
t77 = t59 * mrSges(3,1) - t56 * mrSges(3,2);
t107 = t119 * t56 + mrSges(2,1) + t77;
t106 = mrSges(2,2) - mrSges(3,3) + t108;
t104 = m(6) * pkin(4);
t97 = g(3) * t56;
t95 = t55 * t57;
t94 = t55 * t60;
t79 = pkin(2) * t59 + pkin(8) * t56;
t72 = t21 * t59 - t48 * t56;
t71 = t34 * t59 - t53 * t56;
t70 = t59 * t44 - t56 * t61;
t30 = -t55 * t88 + t57 * t58;
t28 = t55 * t89 + t58 * t60;
t33 = t37 + t42;
t31 = t58 * t88 + t95;
t29 = -t58 * t89 + t94;
t1 = [(-t95 * t105 - t31 * mrSges(4,1) - t25 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t30 * mrSges(4,2) - t24 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) + t110 * (t60 * pkin(1) + t57 * pkin(7)) + t106 * t57 + (-m(4) * t79 - m(5) * t70 - m(6) * t71 - m(7) * t72 - t107) * t60) * g(2) + (-t94 * t105 - t29 * mrSges(4,1) - t23 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t28 * mrSges(4,2) - t22 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t79) - m(5) * (-pkin(1) - t70) - m(6) * (-pkin(1) - t71) - m(7) * (-pkin(1) - t72) + t107) * t57 + (t110 * pkin(7) + t106) * t60) * g(1), -g(3) * t77 + (t118 * g(3) + t116 * (mrSges(3,2) + t117)) * t59 + (t117 * g(3) + t116 * (mrSges(3,1) - t118)) * t56 (m(5) * t101 + mrSges(4,1) * t55 + mrSges(4,2) * t58 - t108 + t113) * t97 + (-t29 * mrSges(4,2) - m(6) * (-t35 * t89 - t36 * t60) - m(7) * (-t26 * t89 - t27 * t60) + t114 * t28 + t111) * g(2) + (t31 * mrSges(4,2) - m(6) * (-t35 * t88 + t36 * t57) - m(7) * (-t26 * t88 + t27 * t57) - t114 * t30 + t112) * g(1) (m(6) * t100 - m(7) * t32 + t113) * t97 + (t22 * t104 - m(7) * (t32 * t89 - t33 * t60) + t111) * g(2) + (-t24 * t104 - m(7) * (t32 * t88 + t33 * t57) + t112) * g(1) (t59 * g(3) - t56 * t116) * t115, -g(1) * t102 - g(2) * t103 - t73 * t97];
taug  = t1(:);
