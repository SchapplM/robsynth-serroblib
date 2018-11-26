% Calculate Gravitation load on the joints for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:27:01
% EndTime: 2018-11-23 16:27:01
% DurationCPUTime: 0.73s
% Computational Cost: add. (470->105), mult. (500->121), div. (0->0), fcn. (452->10), ass. (0->60)
t97 = mrSges(6,1) + mrSges(7,1);
t96 = mrSges(6,2) + mrSges(7,2);
t98 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t31 = qJ(4) + qJ(5);
t27 = cos(t31);
t37 = cos(qJ(4));
t28 = t37 * pkin(4);
t19 = pkin(5) * t27 + t28;
t17 = pkin(3) + t19;
t23 = t28 + pkin(3);
t26 = sin(t31);
t35 = sin(qJ(4));
t95 = -m(5) * pkin(3) - m(6) * t23 - m(7) * t17 - t37 * mrSges(5,1) + t35 * mrSges(5,2) + t96 * t26 - t97 * t27;
t39 = -pkin(9) - pkin(8);
t29 = -qJ(6) + t39;
t94 = -m(5) * pkin(8) + m(6) * t39 + m(7) * t29 - t98;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t93 = g(1) * t38 + g(2) * t36;
t81 = m(6) * pkin(4);
t92 = mrSges(5,1) + t81;
t89 = mrSges(6,1) * t26 + t96 * t27;
t30 = pkin(10) + qJ(3);
t25 = cos(t30);
t64 = t27 * t36;
t65 = t26 * t38;
t11 = -t25 * t65 + t64;
t63 = t27 * t38;
t66 = t26 * t36;
t12 = t25 * t63 + t66;
t88 = -t97 * t11 + t96 * t12;
t10 = -t25 * t64 + t65;
t9 = t25 * t66 + t63;
t87 = -t96 * t10 + t97 * t9;
t86 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t85 = m(4) + m(5) + m(6) + m(7);
t75 = pkin(4) * t35;
t18 = pkin(5) * t26 + t75;
t83 = -m(6) * t75 - m(7) * t18;
t24 = sin(t30);
t33 = cos(pkin(10));
t51 = t25 * mrSges(4,1) - t24 * mrSges(4,2);
t82 = m(3) * pkin(1) + t33 * mrSges(3,1) - sin(pkin(10)) * mrSges(3,2) + mrSges(2,1) + t51 + t98 * t24;
t80 = m(7) * pkin(5);
t72 = g(3) * t24;
t70 = t18 * t36;
t62 = t35 * t36;
t61 = t35 * t38;
t60 = t36 * t37;
t59 = t37 * t38;
t52 = pkin(3) * t25 + pkin(8) * t24;
t48 = t17 * t25 - t24 * t29;
t47 = t25 * t23 - t24 * t39;
t15 = -t25 * t61 + t60;
t13 = t25 * t62 + t59;
t34 = -pkin(7) - qJ(2);
t21 = pkin(2) * t33 + pkin(1);
t16 = t25 * t59 + t62;
t14 = -t25 * t60 + t61;
t1 = [(-t62 * t81 - m(7) * t70 - t16 * mrSges(5,1) - t15 * mrSges(5,2) - t85 * (t38 * t21 - t36 * t34) - t97 * t12 - t96 * t11 + t86 * t36 + (-m(5) * t52 - m(6) * t47 - m(7) * t48 - t82) * t38) * g(2) + (-t14 * mrSges(5,1) - t13 * mrSges(5,2) - t96 * t9 - t97 * t10 + (t85 * t34 + t83 + t86) * t38 + (m(4) * t21 - m(5) * (-t21 - t52) - m(6) * (-t21 - t47) - m(7) * (-t21 - t48) + t82) * t36) * g(1) (-g(1) * t36 + g(2) * t38) * (m(3) + t85) -g(3) * t51 + (t95 * g(3) + t93 * (mrSges(4,2) + t94)) * t25 + (t94 * g(3) + t93 * (mrSges(4,1) - t95)) * t24 (mrSges(5,1) * t35 + mrSges(7,1) * t26 + mrSges(5,2) * t37 - t83 + t89) * t72 + (-t14 * mrSges(5,2) - m(7) * (-t19 * t38 - t25 * t70) + t92 * t13 + t87) * g(2) + (t16 * mrSges(5,2) - m(7) * (-t18 * t25 * t38 + t19 * t36) - t92 * t15 + t88) * g(1) (-(-mrSges(7,1) - t80) * t26 + t89) * t72 + (t9 * t80 + t87) * g(2) + (-t11 * t80 + t88) * g(1) (g(3) * t25 - t93 * t24) * m(7)];
taug  = t1(:);
