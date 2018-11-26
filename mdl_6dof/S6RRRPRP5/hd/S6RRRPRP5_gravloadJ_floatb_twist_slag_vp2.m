% Calculate Gravitation load on the joints for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:43:43
% EndTime: 2018-11-23 17:43:44
% DurationCPUTime: 1.05s
% Computational Cost: add. (601->143), mult. (675->165), div. (0->0), fcn. (638->10), ass. (0->74)
t124 = mrSges(6,1) + mrSges(7,1);
t123 = mrSges(6,2) - mrSges(7,3);
t126 = m(6) + m(7);
t125 = mrSges(4,3) + mrSges(5,3);
t122 = -mrSges(6,3) - mrSges(7,2);
t51 = -qJ(4) - pkin(8);
t121 = -m(4) * pkin(8) + m(5) * t51 - t125;
t50 = qJ(3) + pkin(10);
t44 = qJ(5) + t50;
t39 = sin(t44);
t40 = cos(t44);
t55 = cos(qJ(3));
t46 = t55 * pkin(3);
t41 = t46 + pkin(2);
t42 = sin(t50);
t43 = cos(t50);
t52 = sin(qJ(3));
t120 = m(4) * pkin(2) + m(5) * t41 + t55 * mrSges(4,1) + t43 * mrSges(5,1) - t52 * mrSges(4,2) - t42 * mrSges(5,2) - t123 * t39 + t124 * t40;
t53 = sin(qJ(2));
t119 = t122 * t53;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t118 = g(1) * t57 + g(2) * t54;
t117 = -m(7) * qJ(6) - mrSges(7,3);
t105 = m(5) * pkin(3);
t114 = mrSges(2,2) - mrSges(3,3);
t113 = mrSges(4,1) + t105;
t56 = cos(qJ(2));
t84 = t57 * t39;
t13 = -t54 * t40 + t56 * t84;
t85 = t56 * t57;
t14 = t54 * t39 + t40 * t85;
t112 = t123 * t14 + t124 * t13;
t86 = t54 * t56;
t11 = t39 * t86 + t40 * t57;
t12 = t40 * t86 - t84;
t111 = t124 * t11 + t123 * t12;
t110 = -m(3) - m(4) - m(5);
t70 = t56 * mrSges(3,1) - t53 * mrSges(3,2);
t109 = t125 * t53 + mrSges(2,1) + t70;
t108 = m(7) * pkin(5) + t124;
t107 = -mrSges(6,2) - t117;
t100 = pkin(3) * t52;
t97 = g(3) * t53;
t96 = mrSges(6,2) * t40;
t49 = -pkin(9) + t51;
t94 = t49 * t53;
t93 = t52 * t57;
t88 = t53 * t57;
t87 = t54 * t52;
t32 = pkin(4) * t43 + t46;
t30 = pkin(2) + t32;
t22 = t56 * t30;
t31 = pkin(4) * t42 + t100;
t29 = t57 * t31;
t83 = -t56 * t29 + t54 * t32;
t82 = t57 * pkin(1) + t54 * pkin(7);
t77 = t117 * t40 * t53;
t75 = -t11 * pkin(5) + t12 * qJ(6);
t74 = -t31 * t86 - t32 * t57;
t72 = -t13 * pkin(5) + qJ(6) * t14;
t71 = pkin(2) * t56 + pkin(8) * t53;
t66 = pkin(5) * t40 + qJ(6) * t39;
t65 = t41 * t56 - t51 * t53;
t25 = -t52 * t85 + t54 * t55;
t23 = t52 * t86 + t55 * t57;
t47 = t57 * pkin(7);
t26 = t55 * t85 + t87;
t24 = -t55 * t86 + t93;
t18 = t54 * t42 + t43 * t85;
t17 = -t42 * t85 + t54 * t43;
t16 = t42 * t57 - t43 * t86;
t15 = t42 * t86 + t43 * t57;
t1 = [(-t87 * t105 - t26 * mrSges(4,1) - t18 * mrSges(5,1) - t25 * mrSges(4,2) - t17 * mrSges(5,2) + t122 * t88 + t110 * t82 - t126 * (t30 * t85 + t54 * t31 - t49 * t88 + t82) + t114 * t54 - t108 * t14 - t107 * t13 + (-m(4) * t71 - m(5) * t65 - t109) * t57) * g(2) + (-t93 * t105 - t24 * mrSges(4,1) - t16 * mrSges(5,1) - t23 * mrSges(4,2) - t15 * mrSges(5,2) - t126 * (t54 * t94 + t29 + t47) + t114 * t57 + t110 * t47 + t108 * t12 + t107 * t11 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t71) - m(5) * (-pkin(1) - t65) - t126 * (-pkin(1) - t22) + t109 - t119) * t54) * g(1) (-t70 - t126 * (t22 - t94) + t119) * g(3) + ((-m(7) * t66 - t120) * g(3) + t118 * (t126 * t49 + mrSges(3,2) + t121 + t122)) * t56 + (t121 * g(3) + t118 * (mrSges(3,1) + m(6) * t30 - m(7) * (-t30 - t66) + t120)) * t53, -g(3) * ((m(7) * (-pkin(5) * t39 - t31) - t39 * mrSges(7,1)) * t53 - t77) + (m(5) * t100 + m(6) * t31 + mrSges(4,1) * t52 + mrSges(5,1) * t42 + mrSges(6,1) * t39 + mrSges(4,2) * t55 + mrSges(5,2) * t43 + t96) * t97 + (-t24 * mrSges(4,2) + t15 * mrSges(5,1) - t16 * mrSges(5,2) - m(6) * t74 - m(7) * (t74 + t75) + t113 * t23 + t111) * g(2) + (t26 * mrSges(4,2) - t17 * mrSges(5,1) + t18 * mrSges(5,2) - m(6) * t83 - m(7) * (t72 + t83) - t113 * t25 + t112) * g(1) (t56 * g(3) - t118 * t53) * (m(5) + t126) ((t108 * t39 + t96) * t53 + t77) * g(3) + (-m(7) * t75 + t111) * g(2) + (-m(7) * t72 + t112) * g(1) (-g(1) * t13 - g(2) * t11 - t39 * t97) * m(7)];
taug  = t1(:);
