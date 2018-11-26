% Calculate Gravitation load on the joints for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:55:19
% EndTime: 2018-11-23 16:55:19
% DurationCPUTime: 0.79s
% Computational Cost: add. (356->114), mult. (520->123), div. (0->0), fcn. (462->10), ass. (0->59)
t98 = mrSges(3,1) - mrSges(4,2);
t97 = -mrSges(3,2) + mrSges(4,3);
t96 = -mrSges(6,3) - mrSges(7,3);
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t86 = g(1) * t40 + g(2) * t38;
t33 = pkin(10) + qJ(5);
t23 = sin(t33);
t34 = sin(pkin(10));
t76 = pkin(4) * t34;
t14 = pkin(5) * t23 + t76;
t25 = qJ(6) + t33;
t20 = sin(t25);
t21 = cos(t25);
t24 = cos(t33);
t35 = cos(pkin(10));
t48 = t34 * mrSges(5,1) + t35 * mrSges(5,2);
t95 = -m(6) * t76 - m(7) * t14 - t23 * mrSges(6,1) - t20 * mrSges(7,1) - t24 * mrSges(6,2) - t21 * mrSges(7,2) - t48;
t36 = -pkin(8) - qJ(4);
t32 = -pkin(9) + t36;
t94 = -m(5) * (-pkin(2) - qJ(4)) + mrSges(5,3) - m(6) * (-pkin(2) + t36) - m(7) * (-pkin(2) + t32) - t96;
t37 = sin(qJ(2));
t26 = t37 * qJ(3);
t39 = cos(qJ(2));
t60 = t39 * pkin(2) + t26;
t90 = m(4) * t60;
t56 = -m(5) - m(6) - m(7);
t89 = m(4) - t56;
t88 = -m(7) * pkin(5) - mrSges(6,1);
t87 = t97 * t37 + t98 * t39;
t83 = t96 * t39 - t87;
t82 = m(3) + t89;
t22 = t35 * pkin(4) + pkin(3);
t75 = pkin(5) * t24;
t80 = -m(5) * pkin(3) - m(6) * t22 - m(7) * (t22 + t75) - mrSges(5,1) * t35 + mrSges(5,2) * t34 - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t67 = t38 * t20;
t68 = t37 * t40;
t5 = t21 * t68 - t67;
t66 = t38 * t21;
t6 = t20 * t68 + t66;
t78 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t20 * t40 + t37 * t66;
t8 = t21 * t40 - t37 * t67;
t77 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t72 = g(3) * t39;
t71 = mrSges(7,1) * t21;
t70 = t32 * t39;
t69 = t36 * t39;
t65 = t38 * t23;
t64 = t38 * t24;
t61 = t39 * t40;
t59 = t40 * pkin(1) + t38 * pkin(7);
t57 = qJ(4) * t39;
t9 = t24 * t68 - t65;
t11 = t23 * t40 + t37 * t64;
t15 = t39 * t20 * mrSges(7,2);
t12 = t24 * t40 - t37 * t65;
t10 = t23 * t68 + t64;
t1 = [(-m(3) * t59 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) - mrSges(5,3) * t61 - t48 * t68 - t89 * (pkin(2) * t61 + t40 * t26 + t59) + t80 * t38 + (-mrSges(2,1) - m(5) * t57 - m(6) * (t37 * t76 - t69) - m(7) * (t14 * t37 - t70) + t83) * t40) * g(2) + (-t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + (mrSges(2,1) + t90 + (m(5) * qJ(3) + t48 - m(6) * (-qJ(3) - t76) - m(7) * (-qJ(3) - t14)) * t37 + t82 * pkin(1) + t94 * t39 + t87) * t38 + (-t82 * pkin(7) + t80) * t40) * g(1) (-t90 - m(5) * (t57 + t60) - t39 * mrSges(5,3) - m(6) * (t60 - t69) - m(7) * (t60 - t70) + t95 * t37 + t83) * g(3) + ((m(4) * pkin(2) + t94 + t98) * t37 + (-qJ(3) * t89 + t95 - t97) * t39) * t86 (-t86 * t37 + t72) * t89 (g(3) * t37 + t86 * t39) * t56 -(-mrSges(6,1) * t24 + mrSges(6,2) * t23) * t72 - g(3) * (t15 + (-m(7) * t75 - t71) * t39) + (-t12 * mrSges(6,2) + t88 * t11 - t77) * g(2) + (t10 * mrSges(6,2) + t88 * t9 - t78) * g(1), -g(1) * t78 - g(2) * t77 - g(3) * (-t39 * t71 + t15)];
taug  = t1(:);
