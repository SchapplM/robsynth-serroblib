% Calculate Gravitation load on the joints for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:00:57
% EndTime: 2018-11-23 17:00:58
% DurationCPUTime: 0.77s
% Computational Cost: add. (487->122), mult. (502->140), div. (0->0), fcn. (449->12), ass. (0->71)
t86 = m(6) + m(7);
t92 = m(4) + m(5) + t86;
t35 = qJ(2) + pkin(10);
t27 = sin(t35);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t96 = g(1) * t43 + g(2) * t40;
t100 = t27 * t96;
t34 = qJ(4) + pkin(11);
t28 = cos(t34);
t41 = cos(qJ(4));
t31 = t41 * pkin(4);
t19 = pkin(5) * t28 + t31;
t17 = pkin(3) + t19;
t30 = qJ(6) + t34;
t22 = sin(t30);
t23 = cos(t30);
t24 = t31 + pkin(3);
t26 = sin(t34);
t38 = sin(qJ(4));
t99 = -m(5) * pkin(3) - m(6) * t24 - m(7) * t17 - t41 * mrSges(5,1) - t28 * mrSges(6,1) - t23 * mrSges(7,1) + t38 * mrSges(5,2) + t26 * mrSges(6,2) + t22 * mrSges(7,2);
t98 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t29 = cos(t35);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t97 = -t42 * mrSges(3,1) - t29 * mrSges(4,1) + t39 * mrSges(3,2) - t98 * t27;
t87 = m(6) * pkin(4);
t94 = mrSges(5,1) + t87;
t93 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t82 = pkin(4) * t38;
t18 = pkin(5) * t26 + t82;
t89 = -m(6) * t82 - m(7) * t18;
t88 = m(3) * pkin(1) + mrSges(2,1) - t97;
t69 = t40 * t22;
t5 = t23 * t43 + t29 * t69;
t68 = t40 * t23;
t6 = t22 * t43 - t29 * t68;
t85 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t72 = t29 * t43;
t7 = -t22 * t72 + t68;
t8 = t23 * t72 + t69;
t84 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t81 = pkin(8) * t27;
t78 = g(3) * t27;
t32 = t42 * pkin(2);
t36 = -qJ(5) - pkin(8);
t33 = -pkin(9) + t36;
t74 = t27 * t33;
t73 = t27 * t36;
t71 = t38 * t43;
t70 = t40 * t18;
t67 = t40 * t26;
t66 = t40 * t28;
t65 = t40 * t38;
t64 = t40 * t41;
t63 = t41 * t43;
t59 = pkin(3) * t29 + t81;
t54 = -mrSges(7,1) * t22 - mrSges(7,2) * t23;
t53 = t17 * t29 - t74;
t52 = t24 * t29 - t73;
t15 = -t29 * t71 + t64;
t13 = t29 * t65 + t63;
t37 = -qJ(3) - pkin(7);
t25 = t32 + pkin(1);
t16 = t29 * t63 + t65;
t14 = -t29 * t64 + t71;
t12 = t28 * t72 + t67;
t11 = -t26 * t72 + t66;
t10 = t26 * t43 - t29 * t66;
t9 = t28 * t43 + t29 * t67;
t1 = [(-t65 * t87 - m(7) * t70 - t16 * mrSges(5,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t15 * mrSges(5,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t92 * (t43 * t25 - t40 * t37) + t93 * t40 + (-m(5) * t59 - m(6) * t52 - m(7) * t53 - t88) * t43) * g(2) + (-t14 * mrSges(5,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t13 * mrSges(5,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t92 * t37 + t89 + t93) * t43 + (m(4) * t25 - m(5) * (-t25 - t59) - m(6) * (-t25 - t52) - m(7) * (-t25 - t53) + t88) * t40) * g(1) (mrSges(4,1) - t99) * t100 + (-m(4) * t32 - m(5) * (t32 + t81) - m(6) * (t32 - t73) - m(7) * (t32 - t74) + t97 + t99 * t29) * g(3) + (mrSges(3,2) * t42 + (-m(5) * pkin(8) + m(6) * t36 + m(7) * t33 - t98) * t29 + (t92 * pkin(2) + mrSges(3,1)) * t39) * t96 (-g(1) * t40 + g(2) * t43) * t92 (mrSges(5,1) * t38 + mrSges(6,1) * t26 + mrSges(5,2) * t41 + mrSges(6,2) * t28 - t54 - t89) * t78 + (-t14 * mrSges(5,2) + t9 * mrSges(6,1) - t10 * mrSges(6,2) - m(7) * (-t19 * t43 - t29 * t70) - t85 + t94 * t13) * g(2) + (t16 * mrSges(5,2) - t11 * mrSges(6,1) + t12 * mrSges(6,2) - m(7) * (-t18 * t72 + t40 * t19) - t84 - t94 * t15) * g(1) (t29 * g(3) - t100) * t86, -g(1) * t84 - g(2) * t85 - t54 * t78];
taug  = t1(:);
