% Calculate Gravitation load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:01
% EndTime: 2019-03-09 06:57:03
% DurationCPUTime: 0.71s
% Computational Cost: add. (530->112), mult. (444->131), div. (0->0), fcn. (389->12), ass. (0->65)
t121 = mrSges(6,3) + mrSges(7,3);
t43 = qJ(5) + qJ(6);
t36 = sin(t43);
t44 = qJ(3) + qJ(4);
t37 = sin(t44);
t39 = cos(t44);
t45 = sin(qJ(5));
t120 = t39 * (-m(6) * pkin(9) - t121) + (-mrSges(6,2) * t45 - mrSges(7,2) * t36) * t37;
t38 = cos(t43);
t48 = cos(qJ(5));
t116 = (mrSges(6,1) * t48 + mrSges(7,1) * t38) * t37;
t115 = t39 * mrSges(5,1) + (-mrSges(5,2) + t121) * t37;
t109 = pkin(4) * t39 + pkin(9) * t37;
t32 = pkin(5) * t48 + pkin(4);
t51 = -pkin(10) - pkin(9);
t111 = t32 * t39 - t37 * t51;
t114 = -m(6) * t109 - m(7) * t111;
t99 = m(7) * pkin(5);
t113 = m(3) + m(4);
t112 = mrSges(6,1) + t99;
t60 = -t32 * t37 - t39 * t51;
t95 = pkin(4) * t37;
t46 = sin(qJ(3));
t96 = pkin(3) * t46;
t108 = -m(7) * (t60 - t96) - m(6) * (-t95 - t96) + t116;
t42 = qJ(1) + pkin(11);
t34 = sin(t42);
t35 = cos(t42);
t107 = g(1) * t35 + g(2) * t34;
t106 = -m(5) - m(6) - m(7);
t81 = t39 * t48;
t82 = t39 * t45;
t83 = t38 * t39;
t86 = t36 * t39;
t105 = -mrSges(6,1) * t81 - mrSges(7,1) * t83 + mrSges(6,2) * t82 + mrSges(7,2) * t86 - t115;
t104 = t120 * t34;
t103 = t120 * t35;
t102 = m(6) * t95 - m(7) * t60 + t116;
t49 = cos(qJ(3));
t66 = t49 * mrSges(4,1) - t46 * mrSges(4,2);
t101 = m(4) * pkin(2) + mrSges(3,1) + t115 + t66;
t100 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t5 = t34 * t86 + t35 * t38;
t6 = -t34 * t83 + t35 * t36;
t98 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t7 = t34 * t38 - t35 * t86;
t8 = t34 * t36 + t35 * t83;
t97 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t92 = g(3) * t37;
t47 = sin(qJ(1));
t91 = t47 * pkin(1);
t40 = t49 * pkin(3);
t50 = cos(qJ(1));
t41 = t50 * pkin(1);
t89 = t34 * t45;
t87 = t35 * t45;
t63 = mrSges(5,1) * t37 + mrSges(5,2) * t39;
t62 = -mrSges(7,1) * t36 - mrSges(7,2) * t38;
t11 = t34 * t48 - t35 * t82;
t9 = t34 * t82 + t35 * t48;
t52 = -pkin(8) - pkin(7);
t33 = t40 + pkin(2);
t12 = t35 * t81 + t89;
t10 = -t34 * t81 + t87;
t1 = [(-t89 * t99 - t50 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t47 * mrSges(2,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + t106 * (t33 * t35 - t34 * t52 + t41) - t113 * t41 + t100 * t34 + (-t101 + t114) * t35) * g(2) + (-t87 * t99 + t47 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t50 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t113 * t91 + t106 * (-t35 * t52 - t91) + t100 * t35 + (m(5) * t33 - m(6) * (-t33 - t109) - m(7) * (-t33 - t111) + t101) * t34) * g(1) (t106 - t113) * g(3) (t108 * t34 + t104) * g(2) + (t108 * t35 + t103) * g(1) + (-t66 - m(5) * t40 - m(6) * (t40 + t109) - m(7) * (t40 + t111) + t105) * g(3) + (m(5) * t96 + mrSges(4,1) * t46 + mrSges(4,2) * t49 + t63) * t107, t107 * t63 + (t102 * t34 + t104) * g(2) + (t102 * t35 + t103) * g(1) + (t105 + t114) * g(3) (mrSges(6,2) * t48 + t112 * t45 - t62) * t92 + (-t10 * mrSges(6,2) + t112 * t9 - t98) * g(2) + (t12 * mrSges(6,2) - t11 * t112 - t97) * g(1), -g(1) * t97 - g(2) * t98 - t62 * t92];
taug  = t1(:);
