% Calculate Gravitation load on the joints for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:04:42
% EndTime: 2018-11-23 15:04:43
% DurationCPUTime: 0.68s
% Computational Cost: add. (1088->107), mult. (986->141), div. (0->0), fcn. (928->18), ass. (0->67)
t122 = mrSges(6,2) - mrSges(7,3);
t60 = sin(qJ(6));
t62 = cos(qJ(6));
t121 = -t62 * mrSges(7,1) + t60 * mrSges(7,2) - mrSges(6,1);
t55 = sin(pkin(11));
t56 = sin(pkin(6));
t106 = t55 * t56;
t96 = pkin(6) + qJ(2);
t83 = sin(t96);
t78 = t83 / 0.2e1;
t97 = pkin(6) - qJ(2);
t84 = sin(t97);
t71 = t78 - t84 / 0.2e1;
t63 = cos(qJ(2));
t98 = cos(pkin(11));
t86 = t98 * t63;
t31 = -t55 * t71 + t86;
t53 = pkin(12) + qJ(4);
t48 = sin(t53);
t49 = cos(t53);
t120 = t49 * t106 - t31 * t48;
t80 = cos(t96) / 0.2e1;
t85 = cos(t97);
t40 = t80 - t85 / 0.2e1;
t58 = cos(pkin(6));
t118 = t40 * t48 + t49 * t58;
t117 = -m(6) - m(7);
t50 = qJ(5) + t53;
t45 = sin(t50);
t46 = cos(t50);
t22 = t40 * t45 + t46 * t58;
t23 = -t40 * t46 + t45 * t58;
t116 = t121 * t22 + t122 * t23;
t13 = t46 * t106 - t31 * t45;
t14 = t45 * t106 + t31 * t46;
t115 = t121 * t13 + t122 * t14;
t105 = t55 * t63;
t28 = t98 * t71 + t105;
t87 = t56 * t98;
t11 = -t28 * t45 - t46 * t87;
t12 = t28 * t46 - t45 * t87;
t114 = t121 * t11 + t122 * t12;
t57 = cos(pkin(12));
t47 = t57 * pkin(3) + pkin(2);
t113 = m(4) * pkin(2) + m(5) * t47 + t57 * mrSges(4,1) + t49 * mrSges(5,1) - sin(pkin(12)) * mrSges(4,2) - t48 * mrSges(5,2) + mrSges(3,1) + (m(7) * pkin(5) - t121) * t46 + (m(7) * pkin(10) - t122) * t45;
t59 = -pkin(8) - qJ(3);
t112 = -m(4) * qJ(3) + m(5) * t59 - mrSges(7,1) * t60 - mrSges(7,2) * t62 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t94 = m(4) + m(5) - t117;
t92 = t11 * pkin(5) + t12 * pkin(10);
t90 = t13 * pkin(5) + pkin(10) * t14;
t89 = t22 * pkin(5) + pkin(10) * t23;
t82 = t120 * pkin(4);
t81 = t118 * pkin(4);
t79 = t84 / 0.2e1;
t72 = t79 - t83 / 0.2e1;
t67 = -t28 * t48 - t49 * t87;
t66 = t67 * pkin(4);
t65 = t85 / 0.2e1 + t80;
t61 = sin(qJ(2));
t52 = -pkin(9) + t59;
t39 = t78 + t79;
t37 = pkin(4) * t49 + t47;
t32 = t55 * t72 + t86;
t30 = t55 * t65 + t98 * t61;
t29 = -t98 * t72 + t105;
t27 = t55 * t61 - t98 * t65;
t1 = [(-m(2) - m(3) - t94) * g(3) (t117 * (t39 * t37 + t40 * t52) - t112 * t40 - t113 * t39) * g(3) + (t117 * (-t27 * t37 - t29 * t52) + t112 * t29 + t113 * t27) * g(2) + (t117 * (-t30 * t37 - t32 * t52) + t112 * t32 + t113 * t30) * g(1) (-g(1) * t30 - g(2) * t27 + g(3) * t39) * t94 (-t118 * mrSges(5,1) - (t40 * t49 - t48 * t58) * mrSges(5,2) - m(6) * t81 - m(7) * (t81 + t89) + t116) * g(3) + (-t67 * mrSges(5,1) - (-t28 * t49 + t48 * t87) * mrSges(5,2) - m(6) * t66 - m(7) * (t66 + t92) + t114) * g(2) + (-t120 * mrSges(5,1) - (-t48 * t106 - t31 * t49) * mrSges(5,2) - m(6) * t82 - m(7) * (t82 + t90) + t115) * g(1) (-m(7) * t89 + t116) * g(3) + (-m(7) * t92 + t114) * g(2) + (-m(7) * t90 + t115) * g(1), -g(1) * ((-t14 * t60 + t30 * t62) * mrSges(7,1) + (-t14 * t62 - t30 * t60) * mrSges(7,2)) - g(2) * ((-t12 * t60 + t27 * t62) * mrSges(7,1) + (-t12 * t62 - t27 * t60) * mrSges(7,2)) - g(3) * ((-t23 * t60 - t39 * t62) * mrSges(7,1) + (-t23 * t62 + t39 * t60) * mrSges(7,2))];
taug  = t1(:);
