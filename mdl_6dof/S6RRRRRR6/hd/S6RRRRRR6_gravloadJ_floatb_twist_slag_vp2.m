% Calculate Gravitation load on the joints for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:47
% EndTime: 2019-03-10 04:11:51
% DurationCPUTime: 1.82s
% Computational Cost: add. (1034->165), mult. (1663->217), div. (0->0), fcn. (1941->14), ass. (0->79)
t80 = qJ(5) + qJ(6);
t75 = sin(t80);
t77 = cos(t80);
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t189 = -t87 * mrSges(6,1) - t77 * mrSges(7,1) + t83 * mrSges(6,2) + t75 * mrSges(7,2) - mrSges(5,1);
t73 = pkin(5) * t87 + pkin(4);
t193 = m(6) * pkin(4) + m(7) * t73 - t189;
t183 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t90 = -pkin(12) - pkin(11);
t96 = -m(6) * pkin(11) + m(7) * t90 + t183;
t81 = qJ(3) + qJ(4);
t76 = sin(t81);
t78 = cos(t81);
t84 = sin(qJ(3));
t88 = cos(qJ(3));
t162 = m(4) * pkin(2) + mrSges(4,1) * t88 - mrSges(4,2) * t84 + t193 * t78 - t96 * t76 + mrSges(3,1);
t104 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t190 = t75 * mrSges(7,1) + t87 * mrSges(6,2) + t77 * mrSges(7,2) + t104;
t181 = m(7) * pkin(5);
t127 = t83 * t181;
t161 = -t83 * mrSges(6,1) - t127 - t190;
t128 = cos(pkin(6));
t82 = sin(pkin(6));
t85 = sin(qJ(2));
t150 = t82 * t85;
t180 = t128 * t88 - t84 * t150;
t148 = t82 * t88;
t86 = sin(qJ(1));
t114 = t86 * t128;
t154 = cos(qJ(1));
t89 = cos(qJ(2));
t59 = -t85 * t114 + t154 * t89;
t35 = t86 * t148 - t59 * t84;
t121 = t82 * t154;
t110 = t128 * t154;
t57 = t110 * t85 + t86 * t89;
t29 = -t78 * t121 - t57 * t76;
t30 = -t76 * t121 + t57 * t78;
t179 = t183 * t30 + t189 * t29;
t149 = t82 * t86;
t33 = -t78 * t149 + t59 * t76;
t34 = t76 * t149 + t59 * t78;
t178 = t183 * t34 - t189 * t33;
t50 = t128 * t78 - t76 * t150;
t51 = t128 * t76 + t78 * t150;
t177 = t183 * t51 + t189 * t50;
t169 = mrSges(6,1) + t181;
t166 = m(5) + m(6) + m(7);
t56 = -t110 * t89 + t85 * t86;
t159 = (-t30 * t75 + t56 * t77) * mrSges(7,1) + (-t30 * t77 - t56 * t75) * mrSges(7,2);
t58 = t89 * t114 + t154 * t85;
t5 = -t34 * t75 + t58 * t77;
t6 = t34 * t77 + t58 * t75;
t158 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t155 = t29 * t73 - t30 * t90;
t147 = t82 * t89;
t143 = -t33 * t73 - t34 * t90;
t141 = (-t77 * t147 - t51 * t75) * mrSges(7,1) + (t75 * t147 - t51 * t77) * mrSges(7,2);
t135 = t50 * t73 - t51 * t90;
t129 = t154 * pkin(1) + pkin(8) * t149;
t126 = t84 * t149;
t119 = -pkin(1) * t86 + pkin(8) * t121;
t118 = t29 * pkin(4) + t30 * pkin(11);
t117 = -t33 * pkin(4) + pkin(11) * t34;
t116 = t50 * pkin(4) + pkin(11) * t51;
t68 = t84 * t121;
t115 = -t57 * t88 + t68;
t112 = t35 * pkin(3);
t74 = pkin(3) * t88 + pkin(2);
t91 = -pkin(10) - pkin(9);
t111 = pkin(3) * t126 - t58 * t91 + t59 * t74 + t129;
t7 = -t34 * t83 + t58 * t87;
t103 = t180 * pkin(3);
t98 = t121 * t88 + t57 * t84;
t94 = t98 * pkin(3);
t36 = t59 * t88 + t126;
t8 = t34 * t87 + t58 * t83;
t1 = [(-t154 * mrSges(2,1) - m(3) * t129 - t59 * mrSges(3,1) - m(4) * (pkin(2) * t59 + t129) - t36 * mrSges(4,1) - t35 * mrSges(4,2) - m(5) * t111 - t34 * mrSges(5,1) - m(6) * (pkin(4) * t34 + t111) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t34 * t73 + t111) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t82 + mrSges(2,2)) * t86 + (-t104 - t127) * t58 + t96 * t33) * g(2) + (t86 * mrSges(2,1) + t154 * mrSges(2,2) - m(3) * t119 + t57 * mrSges(3,1) - mrSges(3,3) * t121 - m(4) * (-pkin(2) * t57 + t119) - t115 * mrSges(4,1) - t98 * mrSges(4,2) + t193 * t30 + (t169 * t83 + t190) * t56 + t96 * t29 + t166 * (-pkin(3) * t68 - t56 * t91 + t57 * t74 - t119)) * g(1) (-t166 * (-t56 * t74 - t57 * t91) + t161 * t57 + t162 * t56) * g(2) + (-t166 * (-t58 * t74 - t59 * t91) + t161 * t59 + t162 * t58) * g(1) + (-t166 * t74 * t147 + (-t162 * t89 + (t166 * t91 + t161) * t85) * t82) * g(3) (-t180 * mrSges(4,1) - (-t128 * t84 - t85 * t148) * mrSges(4,2) - m(5) * t103 - m(6) * (t103 + t116) - m(7) * (t103 + t135) + t177) * g(3) + (mrSges(4,1) * t98 - mrSges(4,2) * t115 + m(5) * t94 - m(6) * (t118 - t94) - m(7) * (-t94 + t155) + t179) * g(2) + (-mrSges(4,1) * t35 + mrSges(4,2) * t36 - m(5) * t112 - m(6) * (t112 + t117) - m(7) * (t112 + t143) + t178) * g(1) (-m(6) * t116 - m(7) * t135 + t177) * g(3) + (-m(6) * t118 - m(7) * t155 + t179) * g(2) + (-m(6) * t117 - m(7) * t143 + t178) * g(1) (-(t83 * t147 - t51 * t87) * mrSges(6,2) - t141 - t169 * (-t87 * t147 - t51 * t83)) * g(3) + (-(-t30 * t87 - t56 * t83) * mrSges(6,2) - t159 - t169 * (-t30 * t83 + t56 * t87)) * g(2) + (t8 * mrSges(6,2) - t169 * t7 - t158) * g(1), -g(1) * t158 - g(2) * t159 - g(3) * t141];
taug  = t1(:);
