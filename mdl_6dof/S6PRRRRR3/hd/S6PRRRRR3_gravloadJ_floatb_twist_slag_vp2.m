% Calculate Gravitation load on the joints for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:33:33
% EndTime: 2018-11-23 15:33:34
% DurationCPUTime: 1.19s
% Computational Cost: add. (1474->120), mult. (1590->155), div. (0->0), fcn. (1576->18), ass. (0->63)
t48 = qJ(4) + qJ(5);
t44 = cos(t48);
t54 = cos(qJ(4));
t46 = t54 * pkin(4);
t34 = pkin(5) * t44 + t46;
t45 = qJ(6) + t48;
t40 = sin(t45);
t41 = cos(t45);
t43 = sin(t48);
t51 = sin(qJ(4));
t127 = mrSges(4,1) + m(7) * (pkin(3) + t34) + t41 * mrSges(7,1) - t40 * mrSges(7,2) + m(6) * (t46 + pkin(3)) + t44 * mrSges(6,1) - t43 * mrSges(6,2) + m(5) * pkin(3) + t54 * mrSges(5,1) - t51 * mrSges(5,2);
t57 = -pkin(10) - pkin(9);
t111 = mrSges(4,2) + m(7) * (-pkin(11) + t57) - mrSges(7,3) + m(6) * t57 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t125 = -m(4) - m(5);
t113 = m(6) + m(7) - t125;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t121 = pkin(2) * t113 - t111 * t52 + t127 * t55 + mrSges(3,1);
t100 = pkin(4) * t51;
t33 = pkin(5) * t43 + t100;
t107 = -t51 * mrSges(5,1) - t54 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + t125 * pkin(8) - m(6) * (pkin(8) + t100) - t43 * mrSges(6,1) - t44 * mrSges(6,2) - m(7) * (pkin(8) + t33) - t40 * mrSges(7,1) - t41 * mrSges(7,2);
t92 = pkin(6) + qJ(2);
t82 = cos(t92) / 0.2e1;
t93 = pkin(6) - qJ(2);
t85 = cos(t93);
t32 = t82 - t85 / 0.2e1;
t50 = cos(pkin(6));
t26 = -t32 * t55 + t50 * t52;
t83 = sin(t92);
t80 = t83 / 0.2e1;
t84 = sin(t93);
t81 = t84 / 0.2e1;
t31 = t80 + t81;
t72 = -t26 * t43 - t31 * t44;
t98 = (-t26 * t40 - t31 * t41) * mrSges(7,1) + (-t26 * t41 + t31 * t40) * mrSges(7,2);
t119 = -t72 * mrSges(6,1) - (-t26 * t44 + t31 * t43) * mrSges(6,2) - t98;
t118 = -m(6) * pkin(4) - mrSges(5,1);
t49 = sin(pkin(12));
t68 = t80 - t84 / 0.2e1;
t56 = cos(qJ(2));
t95 = cos(pkin(12));
t86 = t95 * t56;
t23 = -t49 * t68 + t86;
t94 = sin(pkin(6));
t87 = t49 * t94;
t16 = t23 * t55 + t52 * t87;
t53 = sin(qJ(2));
t62 = t85 / 0.2e1 + t82;
t22 = t49 * t62 + t95 * t53;
t103 = (-t16 * t40 + t22 * t41) * mrSges(7,1) + (-t16 * t41 - t22 * t40) * mrSges(7,2);
t74 = -t16 * t43 + t22 * t44;
t117 = -t74 * mrSges(6,1) - (-t16 * t44 - t22 * t43) * mrSges(6,2) - t103;
t97 = t49 * t56;
t20 = t95 * t68 + t97;
t70 = t95 * t94;
t14 = t20 * t55 - t52 * t70;
t19 = t49 * t53 - t95 * t62;
t104 = (-t14 * t40 + t19 * t41) * mrSges(7,1) + (-t14 * t41 - t19 * t40) * mrSges(7,2);
t76 = -t14 * t43 + t19 * t44;
t116 = -t76 * mrSges(6,1) - (-t14 * t44 - t19 * t43) * mrSges(6,2) - t104;
t105 = m(7) * pkin(5);
t69 = t81 - t83 / 0.2e1;
t1 = [(-m(2) - m(3) - t113) * g(3) (-t107 * t32 - t121 * t31) * g(3) + (t107 * (-t95 * t69 + t97) + t121 * t19) * g(2) + (t107 * (t49 * t69 + t86) + t121 * t22) * g(1) (t111 * t26 - t127 * (t32 * t52 + t50 * t55)) * g(3) + (t111 * t14 - t127 * (-t20 * t52 - t55 * t70)) * g(2) + (t111 * t16 - t127 * (-t23 * t52 + t55 * t87)) * g(1) (-(-t26 * t54 + t31 * t51) * mrSges(5,2) - m(7) * (-t26 * t33 - t31 * t34) + t118 * (-t26 * t51 - t31 * t54) + t119) * g(3) + (-(-t14 * t54 - t19 * t51) * mrSges(5,2) - m(7) * (-t14 * t33 + t19 * t34) + t118 * (-t14 * t51 + t19 * t54) + t116) * g(2) + (-(-t16 * t54 - t22 * t51) * mrSges(5,2) - m(7) * (-t16 * t33 + t22 * t34) + t118 * (-t16 * t51 + t22 * t54) + t117) * g(1) (-t72 * t105 + t119) * g(3) + (-t76 * t105 + t116) * g(2) + (-t74 * t105 + t117) * g(1), -g(1) * t103 - g(2) * t104 - g(3) * t98];
taug  = t1(:);
