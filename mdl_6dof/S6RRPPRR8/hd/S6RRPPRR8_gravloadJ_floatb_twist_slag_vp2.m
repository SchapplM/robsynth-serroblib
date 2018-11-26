% Calculate Gravitation load on the joints for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2018-11-23 16:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:53:59
% EndTime: 2018-11-23 16:54:00
% DurationCPUTime: 1.15s
% Computational Cost: add. (398->128), mult. (832->147), div. (0->0), fcn. (864->10), ass. (0->69)
t117 = mrSges(4,1) + mrSges(5,1);
t103 = -mrSges(4,2) + mrSges(5,3);
t114 = -mrSges(5,2) - mrSges(4,3);
t116 = -mrSges(6,3) - mrSges(7,3);
t96 = m(6) * pkin(8) - t116;
t115 = t96 + t114;
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t102 = g(1) * t41 + g(2) * t38;
t34 = sin(pkin(10));
t35 = cos(pkin(10));
t39 = cos(qJ(5));
t36 = sin(qJ(5));
t81 = t34 * t36;
t52 = t35 * t39 + t81;
t53 = t34 * t39 - t35 * t36;
t33 = qJ(5) + qJ(6);
t26 = sin(t33);
t27 = cos(t33);
t54 = t26 * t34 + t27 * t35;
t55 = t26 * t35 - t27 * t34;
t113 = t52 * mrSges(6,1) + t54 * mrSges(7,1) + t53 * mrSges(6,2) - t55 * mrSges(7,2) + t103 * t34 + t117 * t35;
t42 = -pkin(9) - pkin(8);
t112 = m(7) * t42;
t106 = -m(6) - m(7);
t101 = m(5) - t106;
t100 = m(4) + t101;
t109 = pkin(3) * t35 + qJ(4) * t34;
t37 = sin(qJ(2));
t108 = t114 * t37;
t91 = m(7) * pkin(5);
t105 = -mrSges(6,1) - t91;
t104 = mrSges(2,2) - mrSges(3,3);
t40 = cos(qJ(2));
t63 = t40 * mrSges(3,1) - t37 * mrSges(3,2);
t98 = -t112 * t37 - t63;
t25 = pkin(5) * t39 + pkin(4);
t92 = m(6) * pkin(4) + m(7) * t25 + t117;
t76 = t38 * t40;
t13 = t34 * t76 + t35 * t41;
t74 = t41 * t34;
t14 = t35 * t76 - t74;
t58 = -t13 * t26 - t14 * t27;
t59 = t13 * t27 - t14 * t26;
t90 = t59 * mrSges(7,1) + t58 * mrSges(7,2);
t15 = -t38 * t35 + t40 * t74;
t75 = t40 * t41;
t16 = t34 * t38 + t35 * t75;
t5 = t15 * t27 - t16 * t26;
t6 = t15 * t26 + t16 * t27;
t89 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t87 = pkin(4) * t35;
t84 = g(3) * t37;
t30 = t40 * pkin(2);
t83 = (-mrSges(7,1) * t55 - mrSges(7,2) * t54) * t37;
t82 = t13 * t36;
t72 = t41 * pkin(1) + t38 * pkin(7);
t28 = t37 * qJ(3);
t68 = -pkin(1) - t30;
t66 = pkin(2) * t75 + t41 * t28 + t72;
t64 = t16 * pkin(3) + t66;
t57 = t13 * t39 - t14 * t36;
t56 = -t14 * t39 - t82;
t7 = t15 * t39 - t16 * t36;
t51 = pkin(5) * t81 + t25 * t35;
t49 = -pkin(2) - t109;
t31 = t41 * pkin(7);
t8 = t15 * t36 + t16 * t39;
t1 = [(-m(3) * t72 - m(4) * t66 - m(7) * t64 - t8 * mrSges(6,1) - t6 * mrSges(7,1) - t7 * mrSges(6,2) - t5 * mrSges(7,2) + (-m(5) - m(6)) * (t15 * qJ(4) + t64) + t104 * t38 - t92 * t16 + (-m(7) * (pkin(5) * t36 + qJ(4)) - t103) * t15 + (t115 * t37 - mrSges(2,1) + t98) * t41) * g(2) + (t82 * t91 - t56 * mrSges(6,1) - t58 * mrSges(7,1) + t57 * mrSges(6,2) + t59 * mrSges(7,2) - t101 * (-t14 * pkin(3) - t13 * qJ(4) + t31) + t104 * t41 + (-m(3) - m(4)) * t31 + t92 * t14 + t103 * t13 + (m(3) * pkin(1) + mrSges(2,1) + t63 + t106 * t68 + (-m(4) - m(5)) * (t68 - t28) + (-m(6) * (pkin(8) - qJ(3)) - m(7) * (-qJ(3) - t42) + t116) * t37 - t108) * t38) * g(1) (t96 * t37 + t98 - t100 * (t30 + t28) + (-m(6) * t87 - m(7) * t51 - t101 * t109 - t113) * t40 + t108) * g(3) + ((mrSges(3,1) - m(7) * (t49 - t51) - m(6) * (t49 - t87) - m(5) * t49 + m(4) * pkin(2) + t113) * t37 + (-qJ(3) * t100 + mrSges(3,2) - t112 + t115) * t40) * t102 (g(3) * t40 - t102 * t37) * t100, t101 * (-g(1) * t15 - g(2) * t13 - t34 * t84) -(mrSges(6,1) * t53 - mrSges(6,2) * t52) * t84 - g(3) * (t37 * t53 * t91 + t83) + (-t56 * mrSges(6,2) + t105 * t57 - t90) * g(2) + (t8 * mrSges(6,2) + t105 * t7 - t89) * g(1), -g(1) * t89 - g(2) * t90 - g(3) * t83];
taug  = t1(:);
