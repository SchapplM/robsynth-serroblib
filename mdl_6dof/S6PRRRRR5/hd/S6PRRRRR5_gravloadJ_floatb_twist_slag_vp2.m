% Calculate Gravitation load on the joints for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:35:22
% EndTime: 2018-11-23 15:35:23
% DurationCPUTime: 1.46s
% Computational Cost: add. (4227->143), mult. (4336->205), div. (0->0), fcn. (4249->24), ass. (0->83)
t151 = m(5) + m(6) + m(7);
t153 = -m(7) * pkin(5) - mrSges(6,1);
t72 = qJ(5) + qJ(6);
t70 = sin(t72);
t71 = cos(t72);
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t163 = -t70 * mrSges(7,1) - t78 * mrSges(6,2) - t71 * mrSges(7,2) - pkin(10) * t151 + t153 * t75 + mrSges(4,2) - mrSges(5,3);
t160 = mrSges(5,1) + m(7) * (pkin(5) * t78 + pkin(4)) + t71 * mrSges(7,1) - t70 * mrSges(7,2) + m(6) * pkin(4) + t78 * mrSges(6,1) - t75 * mrSges(6,2);
t100 = mrSges(5,2) + m(7) * (-pkin(12) - pkin(11)) - mrSges(7,3) - m(6) * pkin(11) - mrSges(6,3);
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t159 = pkin(3) * t151 - t100 * t76 + t160 * t79 + mrSges(4,1);
t125 = pkin(7) + qJ(3);
t110 = sin(t125) / 0.2e1;
t126 = pkin(7) - qJ(3);
t116 = sin(t126);
t152 = t110 - t116 / 0.2e1;
t127 = pkin(6) + qJ(2);
t114 = cos(t127) / 0.2e1;
t128 = pkin(6) - qJ(2);
t120 = cos(t128);
t64 = t114 - t120 / 0.2e1;
t80 = cos(qJ(3));
t111 = sin(t127) / 0.2e1;
t117 = sin(t128);
t94 = t111 + t117 / 0.2e1;
t157 = t152 * t94 - t64 * t80;
t130 = sin(pkin(13));
t132 = cos(pkin(13));
t63 = t111 - t117 / 0.2e1;
t81 = cos(qJ(2));
t56 = -t130 * t63 + t132 * t81;
t143 = sin(qJ(2));
t97 = t120 / 0.2e1 + t114;
t84 = t130 * t97 + t132 * t143;
t156 = -t152 * t84 + t56 * t80;
t54 = t130 * t81 + t132 * t63;
t83 = t130 * t143 - t132 * t97;
t155 = -t152 * t83 + t54 * t80;
t77 = sin(qJ(3));
t131 = sin(pkin(6));
t92 = t110 + t116 / 0.2e1;
t89 = t92 * t131;
t118 = cos(t125);
t112 = t118 / 0.2e1;
t119 = cos(t126);
t113 = t119 / 0.2e1;
t96 = t113 + t112;
t19 = t132 * t89 + t54 * t77 + t83 * t96;
t105 = t113 - t118 / 0.2e1;
t99 = t105 * t131;
t20 = -t132 * t99 + t155;
t74 = cos(pkin(7));
t121 = t74 * t131;
t73 = sin(pkin(7));
t40 = -t121 * t132 + t73 * t83;
t8 = t20 * t79 + t40 * t76;
t146 = (t19 * t71 - t70 * t8) * mrSges(7,1) + (-t19 * t70 - t71 * t8) * mrSges(7,2);
t23 = t130 * t99 + t156;
t41 = t121 * t130 + t73 * t84;
t10 = t23 * t79 + t41 * t76;
t22 = -t130 * t89 + t56 * t77 + t84 * t96;
t145 = (-t10 * t70 + t22 * t71) * mrSges(7,1) + (-t10 * t71 - t22 * t70) * mrSges(7,2);
t133 = cos(pkin(6));
t35 = t105 * t133 + t157;
t53 = t133 * t74 - t73 * t94;
t16 = t35 * t79 + t53 * t76;
t34 = -t133 * t92 - t64 * t77 - t94 * t96;
t144 = (-t16 * t70 + t34 * t71) * mrSges(7,1) + (-t16 * t71 - t34 * t70) * mrSges(7,2);
t139 = t73 * t76;
t138 = t73 * t79;
t124 = m(4) + t151;
t98 = -mrSges(3,2) + (pkin(9) * t124 + mrSges(4,3)) * t73;
t95 = t112 - t119 / 0.2e1;
t90 = t95 * t131;
t62 = t94 * pkin(2);
t52 = t84 * pkin(2);
t51 = t83 * pkin(2);
t39 = t152 * t64 + t80 * t94;
t32 = -t152 * t56 - t80 * t84;
t30 = -t152 * t54 - t80 * t83;
t1 = [(-m(2) - m(3) - t124) * g(3) (-t94 * mrSges(3,1) - m(4) * t62 - t39 * mrSges(4,1) + t98 * t64 - t160 * (-t139 * t64 + t39 * t79) + t163 * (-t64 * t96 + t77 * t94) + t100 * (t138 * t64 + t39 * t76)) * g(3) + (t83 * mrSges(3,1) + m(4) * t51 - t30 * mrSges(4,1) - t98 * t54 - t160 * (t139 * t54 + t30 * t79) + t163 * (t54 * t96 - t77 * t83) + t100 * (-t138 * t54 + t30 * t76)) * g(2) + (t84 * mrSges(3,1) + m(4) * t52 - t32 * mrSges(4,1) - t98 * t56 - t160 * (t139 * t56 + t32 * t79) + t163 * (t56 * t96 - t77 * t84) + t100 * (-t138 * t56 + t32 * t76)) * g(1) + (-g(2) * (t30 * pkin(3) - t51) - g(1) * (t32 * pkin(3) - t52) - g(3) * (t39 * pkin(3) + t62)) * t151 (t163 * (-t133 * t95 + t157) + t159 * t34) * g(3) + (t163 * (t132 * t90 + t155) + t159 * t19) * g(2) + (t163 * (-t130 * t90 + t156) + t159 * t22) * g(1) (t100 * t16 - t160 * (-t35 * t76 + t53 * t79)) * g(3) + (t100 * t8 - t160 * (-t20 * t76 + t40 * t79)) * g(2) + (-t160 * (-t23 * t76 + t41 * t79) + t100 * t10) * g(1) (-(-t16 * t78 - t34 * t75) * mrSges(6,2) - t144 + t153 * (-t16 * t75 + t34 * t78)) * g(3) + (-(-t19 * t75 - t78 * t8) * mrSges(6,2) - t146 + t153 * (t19 * t78 - t75 * t8)) * g(2) + (-(-t10 * t78 - t22 * t75) * mrSges(6,2) - t145 + t153 * (-t10 * t75 + t22 * t78)) * g(1), -g(1) * t145 - g(2) * t146 - g(3) * t144];
taug  = t1(:);
