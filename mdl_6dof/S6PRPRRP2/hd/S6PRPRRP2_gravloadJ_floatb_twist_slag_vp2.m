% Calculate Gravitation load on the joints for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2018-11-23 15:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:18
% EndTime: 2018-11-23 15:00:19
% DurationCPUTime: 0.90s
% Computational Cost: add. (1413->100), mult. (1198->138), div. (0->0), fcn. (1145->20), ass. (0->70)
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t95 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t97 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t122 = t95 * t72 - t97 * t75 - mrSges(5,1);
t121 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t124 = -m(6) - m(7);
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t131 = t124 * (pkin(4) * t76 + pkin(9) * t73) + t121 * t73 - mrSges(4,1) + t122 * t76;
t130 = -m(5) + t124;
t67 = pkin(6) - qJ(2);
t57 = sin(t67) / 0.2e1;
t66 = pkin(6) + qJ(2);
t60 = sin(t66);
t129 = t60 / 0.2e1 + t57;
t127 = -t97 * t72 - t95 * t75 + mrSges(4,2) - mrSges(5,3);
t58 = cos(t66) / 0.2e1;
t64 = cos(t67);
t50 = t64 / 0.2e1 + t58;
t65 = qJ(2) + pkin(11);
t62 = cos(t65);
t68 = sin(pkin(10));
t115 = t68 * t62;
t74 = sin(qJ(2));
t114 = t68 * t74;
t69 = sin(pkin(6));
t113 = t69 * t73;
t112 = t69 * t76;
t70 = cos(pkin(10));
t111 = t70 * t62;
t110 = t70 * t74;
t98 = pkin(6) + t65;
t89 = sin(t98);
t52 = t89 / 0.2e1;
t99 = pkin(6) - t65;
t90 = sin(t99);
t105 = t52 - t90 / 0.2e1;
t104 = t129 * pkin(2);
t103 = m(4) - t130;
t48 = t50 * pkin(2);
t96 = -pkin(2) * t114 + t70 * t48;
t91 = cos(t98);
t88 = -pkin(2) * t110 - t48 * t68;
t87 = cos(t99) / 0.2e1;
t86 = t90 / 0.2e1;
t79 = t91 / 0.2e1 + t87;
t78 = -t89 / 0.2e1 + t86;
t77 = cos(qJ(2));
t71 = cos(pkin(6));
t59 = sin(t65);
t49 = t57 - t60 / 0.2e1;
t47 = t87 - t91 / 0.2e1;
t46 = t86 + t52;
t37 = t47 * t76 + t71 * t73;
t36 = -t47 * t73 + t71 * t76;
t34 = t68 * t78 + t111;
t33 = -t105 * t68 + t111;
t32 = t70 * t59 + t68 * t79;
t31 = -t70 * t78 + t115;
t30 = t105 * t70 + t115;
t29 = t59 * t68 - t70 * t79;
t16 = t113 * t68 + t33 * t76;
t15 = t112 * t68 - t33 * t73;
t14 = -t113 * t70 + t30 * t76;
t13 = -t112 * t70 - t30 * t73;
t9 = t37 * t72 + t46 * t75;
t3 = t16 * t72 - t32 * t75;
t1 = t14 * t72 - t29 * t75;
t2 = [(-m(2) - m(3) - t103) * g(3) (-t129 * mrSges(3,1) - (t58 - t64 / 0.2e1) * mrSges(3,2) - m(4) * t104 + t130 * (t46 * pkin(3) + pkin(8) * t47 + t104) + t127 * t47 + t131 * t46) * g(3) + (-(t50 * t70 - t114) * mrSges(3,1) - (t70 * t49 - t68 * t77) * mrSges(3,2) - m(4) * t96 + t130 * (-t29 * pkin(3) + pkin(8) * t31 + t96) + t127 * t31 - t131 * t29) * g(2) + (-(-t50 * t68 - t110) * mrSges(3,1) - (-t68 * t49 - t70 * t77) * mrSges(3,2) - m(4) * t88 + t130 * (-t32 * pkin(3) + pkin(8) * t34 + t88) + t127 * t34 - t131 * t32) * g(1) (-g(3) * t71 + (-g(1) * t68 + g(2) * t70) * t69) * t103 (t124 * (t36 * pkin(4) + pkin(9) * t37) + t121 * t37 + t122 * t36) * g(3) + (t124 * (t13 * pkin(4) + pkin(9) * t14) + t121 * t14 + t122 * t13) * g(2) + (t124 * (t15 * pkin(4) + pkin(9) * t16) + t121 * t16 + t122 * t15) * g(1) (t97 * t9 + t95 * (t37 * t75 - t46 * t72)) * g(3) + (t95 * (t14 * t75 + t29 * t72) + t97 * t1) * g(2) + (t95 * (t16 * t75 + t32 * t72) + t97 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
