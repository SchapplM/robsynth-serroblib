% Calculate Gravitation load on the joints for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:11:11
% EndTime: 2018-11-23 15:11:12
% DurationCPUTime: 0.98s
% Computational Cost: add. (1111->105), mult. (1176->146), div. (0->0), fcn. (1123->16), ass. (0->56)
t56 = cos(qJ(5));
t118 = -m(6) * pkin(4) - m(7) * (pkin(5) * t56 + pkin(4)) - mrSges(5,1);
t103 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t117 = mrSges(6,1) + mrSges(7,1);
t111 = -mrSges(6,2) - mrSges(7,2);
t109 = m(5) + m(6) + m(7);
t84 = pkin(6) + qJ(2);
t73 = cos(t84) / 0.2e1;
t85 = pkin(6) - qJ(2);
t79 = cos(t85);
t35 = t73 - t79 / 0.2e1;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t87 = cos(pkin(6));
t116 = t35 * t54 + t87 * t57;
t72 = sin(t84) / 0.2e1;
t78 = sin(t85);
t34 = t72 - t78 / 0.2e1;
t49 = sin(pkin(10));
t58 = cos(qJ(2));
t86 = cos(pkin(10));
t68 = -t49 * t34 + t86 * t58;
t50 = sin(pkin(6));
t91 = t49 * t50;
t115 = -t68 * t54 + t57 * t91;
t48 = qJ(3) + pkin(11);
t46 = sin(t48);
t47 = cos(t48);
t114 = -m(4) * pkin(2) - t57 * mrSges(4,1) + t54 * mrSges(4,2) + t103 * t46 + t118 * t47 - mrSges(3,1);
t102 = m(7) * pkin(5);
t69 = t86 * t34 + t49 * t58;
t81 = t50 * t86;
t110 = -t69 * t54 - t57 * t81;
t53 = sin(qJ(5));
t107 = t111 * t53 + t117 * t56 - t118;
t106 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t105 = -t102 - t117;
t100 = t69 * t53;
t98 = t68 * t53;
t97 = t35 * t53;
t93 = t47 * t53;
t92 = t47 * t56;
t59 = t79 / 0.2e1 + t73;
t55 = sin(qJ(2));
t52 = -qJ(4) - pkin(8);
t45 = pkin(3) * t57 + pkin(2);
t33 = t72 + t78 / 0.2e1;
t26 = t49 * t59 + t86 * t55;
t23 = t49 * t55 - t86 * t59;
t20 = -t35 * t47 + t87 * t46;
t19 = -t35 * t46 - t87 * t47;
t14 = t46 * t91 + t47 * t68;
t13 = t46 * t68 - t47 * t91;
t12 = -t46 * t81 + t47 * t69;
t11 = t46 * t69 + t47 * t81;
t1 = [(-m(2) - m(3) - m(4) - t109) * g(3) (t97 * t102 - t109 * (t33 * t45 + t35 * t52) - t117 * (t33 * t92 - t97) + t111 * (-t33 * t93 - t35 * t56) - t106 * t35 + t114 * t33) * g(3) + (-t100 * t102 - t109 * (-t23 * t45 - t52 * t69) - t117 * (-t23 * t92 + t100) + t111 * (t23 * t93 + t56 * t69) + t106 * t69 - t114 * t23) * g(2) + (-t98 * t102 - t109 * (-t26 * t45 - t52 * t68) - t117 * (-t26 * t92 + t98) + t111 * (t26 * t93 + t56 * t68) + t106 * t68 - t114 * t26) * g(1) (-t116 * mrSges(4,1) - (t35 * t57 - t87 * t54) * mrSges(4,2) + t103 * t20 + t107 * t19) * g(3) + (-(t54 * t81 - t57 * t69) * mrSges(4,2) - mrSges(4,1) * t110 + t103 * t12 + t107 * t11) * g(2) + (-t115 * mrSges(4,1) - (-t54 * t91 - t57 * t68) * mrSges(4,2) + t103 * t14 + t107 * t13) * g(1) + (-g(1) * t115 - t110 * g(2) - g(3) * t116) * t109 * pkin(3), t109 * (-g(1) * t26 - g(2) * t23 + g(3) * t33) (t111 * (-t20 * t56 + t33 * t53) + t105 * (-t20 * t53 - t33 * t56)) * g(3) + (t111 * (-t12 * t56 - t23 * t53) + t105 * (-t12 * t53 + t23 * t56)) * g(2) + (t111 * (-t14 * t56 - t26 * t53) + t105 * (-t14 * t53 + t26 * t56)) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t19) * m(7)];
taug  = t1(:);
