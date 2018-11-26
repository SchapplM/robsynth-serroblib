% Calculate Gravitation load on the joints for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:26:01
% EndTime: 2018-11-23 15:26:02
% DurationCPUTime: 1.34s
% Computational Cost: add. (3772->132), mult. (3875->188), div. (0->0), fcn. (3780->24), ass. (0->83)
t146 = m(6) + m(7);
t143 = m(5) + t146;
t66 = pkin(13) + qJ(6);
t64 = sin(t66);
t65 = cos(t66);
t67 = sin(pkin(13));
t69 = cos(pkin(13));
t154 = -t64 * mrSges(7,1) - t69 * mrSges(6,2) - t65 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t67 - t143 * pkin(10);
t142 = mrSges(5,1) + m(7) * (pkin(5) * t69 + pkin(4)) + mrSges(7,1) * t65 - mrSges(7,2) * t64 + m(6) * pkin(4) + mrSges(6,1) * t69 - mrSges(6,2) * t67;
t96 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(11) - qJ(5)) - mrSges(7,3);
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t151 = pkin(3) * t143 + t142 * t73 - t71 * t96 + mrSges(4,1);
t119 = pkin(7) + qJ(3);
t106 = sin(t119) / 0.2e1;
t120 = pkin(7) - qJ(3);
t111 = sin(t120);
t144 = t106 - t111 / 0.2e1;
t121 = pkin(6) + qJ(2);
t110 = cos(t121) / 0.2e1;
t122 = pkin(6) - qJ(2);
t115 = cos(t122);
t58 = t110 - t115 / 0.2e1;
t74 = cos(qJ(3));
t107 = sin(t121) / 0.2e1;
t112 = sin(t122);
t90 = t107 + t112 / 0.2e1;
t149 = t144 * t90 - t58 * t74;
t124 = sin(pkin(12));
t126 = cos(pkin(12));
t57 = t107 - t112 / 0.2e1;
t75 = cos(qJ(2));
t50 = -t124 * t57 + t126 * t75;
t138 = sin(qJ(2));
t93 = t115 / 0.2e1 + t110;
t79 = t124 * t93 + t126 * t138;
t148 = -t144 * t79 + t50 * t74;
t48 = t124 * t75 + t126 * t57;
t78 = t124 * t138 - t126 * t93;
t147 = -t144 * t78 + t48 * t74;
t68 = sin(pkin(7));
t134 = t68 * t71;
t133 = t68 * t73;
t128 = cos(pkin(6));
t127 = cos(pkin(7));
t125 = sin(pkin(6));
t118 = m(4) + t143;
t114 = cos(t120);
t113 = cos(t119);
t109 = t114 / 0.2e1;
t108 = t113 / 0.2e1;
t103 = t127 * t125;
t102 = t109 - t113 / 0.2e1;
t95 = t102 * t125;
t94 = -mrSges(3,2) + (t118 * pkin(9) + mrSges(4,3)) * t68;
t92 = t109 + t108;
t91 = t108 - t114 / 0.2e1;
t88 = t106 + t111 / 0.2e1;
t86 = t91 * t125;
t85 = t88 * t125;
t81 = t128 * t127 - t90 * t68;
t77 = t124 * t103 + t79 * t68;
t76 = -t126 * t103 + t78 * t68;
t72 = sin(qJ(3));
t56 = t90 * pkin(2);
t47 = t79 * pkin(2);
t46 = t78 * pkin(2);
t33 = t144 * t58 + t90 * t74;
t29 = t128 * t102 + t149;
t28 = -t128 * t88 - t58 * t72 - t90 * t92;
t26 = -t144 * t50 - t79 * t74;
t24 = -t144 * t48 - t78 * t74;
t17 = t124 * t95 + t148;
t16 = -t124 * t85 + t50 * t72 + t79 * t92;
t14 = -t126 * t95 + t147;
t13 = t126 * t85 + t48 * t72 + t78 * t92;
t10 = t29 * t73 + t81 * t71;
t9 = t29 * t71 - t81 * t73;
t4 = t17 * t73 + t77 * t71;
t3 = t17 * t71 - t77 * t73;
t2 = t14 * t73 + t76 * t71;
t1 = t14 * t71 - t76 * t73;
t5 = [(-m(2) - m(3) - t118) * g(3) (-t90 * mrSges(3,1) - m(4) * t56 - t33 * mrSges(4,1) + t94 * t58 - t142 * (-t58 * t134 + t33 * t73) + t154 * (-t58 * t92 + t90 * t72) + t96 * (t58 * t133 + t33 * t71)) * g(3) + (t78 * mrSges(3,1) + m(4) * t46 - t24 * mrSges(4,1) - t94 * t48 + t96 * (-t133 * t48 + t24 * t71) - t142 * (t134 * t48 + t24 * t73) + t154 * (t48 * t92 - t78 * t72)) * g(2) + (t79 * mrSges(3,1) + m(4) * t47 - t26 * mrSges(4,1) - t94 * t50 + t96 * (-t133 * t50 + t26 * t71) - t142 * (t134 * t50 + t26 * t73) + t154 * (t50 * t92 - t79 * t72)) * g(1) + (-g(2) * (t24 * pkin(3) - t46) - g(1) * (t26 * pkin(3) - t47) - g(3) * (t33 * pkin(3) + t56)) * t143 (t154 * (-t128 * t91 + t149) + t151 * t28) * g(3) + (t154 * (t126 * t86 + t147) + t151 * t13) * g(2) + (t154 * (-t124 * t86 + t148) + t151 * t16) * g(1) (t10 * t96 + t142 * t9) * g(3) + (t1 * t142 + t2 * t96) * g(2) + (t142 * t3 + t4 * t96) * g(1), t146 * (-g(1) * t3 - g(2) * t1 - g(3) * t9) -g(1) * ((t16 * t65 - t4 * t64) * mrSges(7,1) + (-t16 * t64 - t4 * t65) * mrSges(7,2)) - g(2) * ((t13 * t65 - t2 * t64) * mrSges(7,1) + (-t13 * t64 - t2 * t65) * mrSges(7,2)) - g(3) * ((-t10 * t64 + t28 * t65) * mrSges(7,1) + (-t10 * t65 - t28 * t64) * mrSges(7,2))];
taug  = t5(:);
