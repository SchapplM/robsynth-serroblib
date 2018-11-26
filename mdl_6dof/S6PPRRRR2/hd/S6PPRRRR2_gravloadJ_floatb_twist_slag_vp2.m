% Calculate Gravitation load on the joints for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:52:14
% EndTime: 2018-11-23 14:52:15
% DurationCPUTime: 0.83s
% Computational Cost: add. (2940->104), mult. (3035->148), div. (0->0), fcn. (2987->24), ass. (0->74)
t116 = m(5) + m(6) + m(7);
t47 = qJ(5) + qJ(6);
t45 = sin(t47);
t46 = cos(t47);
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t123 = mrSges(5,1) + m(7) * (pkin(5) * t55 + pkin(4)) + t46 * mrSges(7,1) - t45 * mrSges(7,2) + m(6) * pkin(4) + t55 * mrSges(6,1) - t52 * mrSges(6,2);
t114 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t117 = -m(7) * pkin(5) - mrSges(6,1);
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t122 = pkin(3) * t116 - t114 * t53 + t123 * t56 + mrSges(4,1);
t96 = pkin(7) + qJ(3);
t83 = sin(t96) / 0.2e1;
t97 = pkin(7) - qJ(3);
t89 = sin(t97);
t118 = t83 - t89 / 0.2e1;
t95 = pkin(6) - pkin(13);
t79 = cos(t95) / 0.2e1;
t94 = pkin(6) + pkin(13);
t88 = cos(t94);
t41 = t79 - t88 / 0.2e1;
t57 = cos(qJ(3));
t78 = sin(t94) / 0.2e1;
t87 = sin(t95);
t65 = t78 + t87 / 0.2e1;
t121 = t118 * t65 + t41 * t57;
t100 = cos(pkin(12));
t40 = t78 - t87 / 0.2e1;
t50 = cos(pkin(13));
t99 = sin(pkin(12));
t34 = t100 * t50 - t40 * t99;
t66 = t79 + t88 / 0.2e1;
t98 = sin(pkin(13));
t60 = t100 * t98 + t66 * t99;
t120 = -t118 * t60 + t34 * t57;
t33 = t100 * t40 + t50 * t99;
t59 = -t100 * t66 + t98 * t99;
t119 = -t118 * t59 + t33 * t57;
t111 = m(3) + m(4) + t116;
t110 = -t45 * mrSges(7,1) - t55 * mrSges(6,2) - t46 * mrSges(7,2) - t116 * pkin(9) + t117 * t52 + mrSges(4,2) - mrSges(5,3);
t54 = sin(qJ(3));
t49 = sin(pkin(6));
t68 = t83 + t89 / 0.2e1;
t63 = t49 * t68;
t90 = cos(t96);
t84 = t90 / 0.2e1;
t91 = cos(t97);
t85 = t91 / 0.2e1;
t71 = t85 + t84;
t15 = t100 * t63 + t33 * t54 + t59 * t71;
t77 = t85 - t90 / 0.2e1;
t72 = t49 * t77;
t16 = -t100 * t72 + t119;
t51 = cos(pkin(7));
t102 = t49 * t51;
t48 = sin(pkin(7));
t25 = -t100 * t102 + t48 * t59;
t8 = t16 * t56 + t25 * t53;
t108 = (t15 * t46 - t45 * t8) * mrSges(7,1) + (-t15 * t45 - t46 * t8) * mrSges(7,2);
t19 = t72 * t99 + t120;
t26 = t102 * t99 + t48 * t60;
t10 = t19 * t56 + t26 * t53;
t18 = t34 * t54 + t60 * t71 - t63 * t99;
t107 = (-t10 * t45 + t18 * t46) * mrSges(7,1) + (-t10 * t46 - t18 * t45) * mrSges(7,2);
t101 = cos(pkin(6));
t23 = t101 * t77 + t121;
t35 = t101 * t51 - t48 * t65;
t12 = t23 * t56 + t35 * t53;
t22 = -t101 * t68 + t41 * t54 - t65 * t71;
t106 = (-t12 * t45 + t22 * t46) * mrSges(7,1) + (-t12 * t46 - t22 * t45) * mrSges(7,2);
t70 = t84 - t91 / 0.2e1;
t64 = t49 * t70;
t1 = [(-m(2) - t111) * g(3) (-t101 * g(3) + (-g(1) * t99 + g(2) * t100) * t49) * t111 (t110 * (-t101 * t70 + t121) + t122 * t22) * g(3) + (t110 * (t100 * t64 + t119) + t122 * t15) * g(2) + (t110 * (-t64 * t99 + t120) + t122 * t18) * g(1) (t114 * t12 - t123 * (-t23 * t53 + t35 * t56)) * g(3) + (t114 * t8 - t123 * (-t16 * t53 + t25 * t56)) * g(2) + (-t123 * (-t19 * t53 + t26 * t56) + t114 * t10) * g(1) (-(-t12 * t55 - t22 * t52) * mrSges(6,2) - t106 + t117 * (-t12 * t52 + t22 * t55)) * g(3) + (-(-t15 * t52 - t55 * t8) * mrSges(6,2) - t108 + t117 * (t15 * t55 - t52 * t8)) * g(2) + (-(-t10 * t55 - t18 * t52) * mrSges(6,2) - t107 + t117 * (-t10 * t52 + t18 * t55)) * g(1), -g(1) * t107 - g(2) * t108 - g(3) * t106];
taug  = t1(:);
