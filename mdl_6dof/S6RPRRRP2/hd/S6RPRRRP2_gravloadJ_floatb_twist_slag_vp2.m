% Calculate Gravitation load on the joints for
% S6RPRRRP2
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:32
% EndTime: 2018-11-23 16:24:33
% DurationCPUTime: 0.73s
% Computational Cost: add. (486->104), mult. (477->126), div. (0->0), fcn. (434->10), ass. (0->59)
t95 = mrSges(6,1) + mrSges(7,1);
t94 = mrSges(6,2) + mrSges(7,2);
t96 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t33 = qJ(4) + qJ(5);
t28 = cos(t33);
t37 = cos(qJ(4));
t29 = t37 * pkin(4);
t19 = pkin(5) * t28 + t29;
t17 = pkin(3) + t19;
t24 = t29 + pkin(3);
t27 = sin(t33);
t34 = sin(qJ(4));
t93 = -m(5) * pkin(3) - m(6) * t24 - m(7) * t17 - t37 * mrSges(5,1) + t34 * mrSges(5,2) + t94 * t27 - t95 * t28;
t40 = -pkin(9) - pkin(8);
t31 = -qJ(6) + t40;
t92 = -m(5) * pkin(8) + m(6) * t40 + m(7) * t31 - t96;
t32 = qJ(1) + pkin(10);
t25 = sin(t32);
t26 = cos(t32);
t91 = g(1) * t26 + g(2) * t25;
t79 = m(6) * pkin(4);
t73 = pkin(4) * t34;
t18 = pkin(5) * t27 + t73;
t90 = m(7) * t18;
t89 = mrSges(5,1) + t79;
t86 = mrSges(6,1) * t27 + t94 * t28;
t38 = cos(qJ(3));
t64 = t27 * t38;
t11 = t25 * t28 - t26 * t64;
t63 = t28 * t38;
t12 = t25 * t27 + t26 * t63;
t85 = -t95 * t11 + t94 * t12;
t10 = -t25 * t63 + t26 * t27;
t9 = t25 * t64 + t26 * t28;
t84 = -t94 * t10 + t95 * t9;
t83 = -m(4) - m(5) - m(6) - m(7);
t81 = mrSges(3,2) - mrSges(4,3) - t90;
t35 = sin(qJ(3));
t51 = t38 * mrSges(4,1) - t35 * mrSges(4,2);
t80 = t96 * t35 + mrSges(3,1) + t51;
t78 = m(7) * pkin(5);
t70 = g(3) * t35;
t36 = sin(qJ(1));
t69 = t36 * pkin(1);
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t67 = t18 * t38;
t66 = t25 * t34;
t65 = t26 * t34;
t62 = t34 * t38;
t58 = t37 * t38;
t52 = t38 * pkin(3) + t35 * pkin(8);
t48 = t38 * t17 - t35 * t31;
t47 = t38 * t24 - t35 * t40;
t15 = t25 * t37 - t26 * t62;
t13 = t25 * t62 + t26 * t37;
t16 = t26 * t58 + t66;
t14 = -t25 * t58 + t65;
t1 = [(-t66 * t79 - m(3) * t30 - t39 * mrSges(2,1) - t16 * mrSges(5,1) + t36 * mrSges(2,2) - t15 * mrSges(5,2) + t83 * (t26 * pkin(2) + t25 * pkin(7) + t30) + t81 * t25 - t95 * t12 - t94 * t11 + (-m(5) * t52 - m(6) * t47 - m(7) * t48 - t80) * t26) * g(2) + (-t65 * t79 + m(3) * t69 + t36 * mrSges(2,1) - t14 * mrSges(5,1) + t39 * mrSges(2,2) - t13 * mrSges(5,2) - t94 * t9 + t83 * (t26 * pkin(7) - t69) + t81 * t26 - t95 * t10 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t52) - m(6) * (-pkin(2) - t47) - m(7) * (-pkin(2) - t48) + t80) * t25) * g(1) (-m(3) + t83) * g(3), -g(3) * t51 + (t93 * g(3) + t91 * (mrSges(4,2) + t92)) * t38 + (t92 * g(3) + t91 * (mrSges(4,1) - t93)) * t35 (m(6) * t73 + mrSges(5,1) * t34 + mrSges(7,1) * t27 + mrSges(5,2) * t37 + t86 + t90) * t70 + (-t14 * mrSges(5,2) - m(7) * (-t19 * t26 - t25 * t67) + t89 * t13 + t84) * g(2) + (t16 * mrSges(5,2) - m(7) * (t19 * t25 - t26 * t67) - t89 * t15 + t85) * g(1) (-(-mrSges(7,1) - t78) * t27 + t86) * t70 + (t9 * t78 + t84) * g(2) + (-t11 * t78 + t85) * g(1) (g(3) * t38 - t91 * t35) * m(7)];
taug  = t1(:);
