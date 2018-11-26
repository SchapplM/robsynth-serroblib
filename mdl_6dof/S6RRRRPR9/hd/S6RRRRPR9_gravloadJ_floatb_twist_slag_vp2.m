% Calculate Gravitation load on the joints for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:18:07
% EndTime: 2018-11-23 18:18:08
% DurationCPUTime: 1.58s
% Computational Cost: add. (1735->154), mult. (1798->192), div. (0->0), fcn. (1782->18), ass. (0->79)
t184 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t86 = -pkin(11) - qJ(5);
t195 = m(6) * qJ(5) - m(7) * t86 - t184;
t81 = pkin(12) + qJ(6);
t76 = sin(t81);
t77 = cos(t81);
t172 = -t77 * mrSges(7,1) + t76 * mrSges(7,2);
t194 = t172 - mrSges(5,1);
t83 = sin(pkin(12));
t85 = cos(pkin(12));
t193 = t85 * mrSges(6,1) - t83 * mrSges(6,2);
t74 = pkin(5) * t85 + pkin(4);
t191 = -m(6) * pkin(4) - m(7) * t74 - t193;
t82 = qJ(3) + qJ(4);
t78 = sin(t82);
t79 = cos(t82);
t87 = sin(qJ(3));
t90 = cos(qJ(3));
t164 = m(4) * pkin(2) + t90 * mrSges(4,1) - t87 * mrSges(4,2) + mrSges(3,1) - (t191 + t194) * t79 + t195 * t78;
t190 = t193 - t194;
t165 = -m(4) * pkin(9) - t85 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t83 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t189 = t76 * mrSges(7,1) + t77 * mrSges(7,2) - t165;
t135 = cos(pkin(6));
t133 = pkin(6) + qJ(2);
t114 = cos(t133) / 0.2e1;
t134 = pkin(6) - qJ(2);
t119 = cos(t134);
t60 = t114 - t119 / 0.2e1;
t182 = t135 * t90 + t60 * t87;
t154 = cos(qJ(1));
t113 = sin(t133) / 0.2e1;
t118 = sin(t134);
t59 = t113 - t118 / 0.2e1;
t89 = sin(qJ(1));
t91 = cos(qJ(2));
t106 = t154 * t91 - t89 * t59;
t84 = sin(pkin(6));
t148 = t84 * t89;
t27 = -t106 * t87 + t90 * t148;
t107 = t154 * t59 + t89 * t91;
t128 = t84 * t154;
t21 = t107 * t78 + t128 * t79;
t22 = t107 * t79 - t78 * t128;
t180 = t184 * t22 + t190 * t21;
t25 = t106 * t78 - t148 * t79;
t26 = t106 * t79 + t148 * t78;
t179 = t184 * t26 + t190 * t25;
t42 = -t135 * t79 - t60 * t78;
t43 = t135 * t78 - t60 * t79;
t178 = t184 * t43 + t190 * t42;
t177 = m(6) + m(7);
t171 = m(5) + t177;
t161 = -t21 * t74 - t22 * t86;
t155 = -t25 * t74 - t26 * t86;
t143 = -t42 * t74 - t43 * t86;
t136 = t154 * pkin(1) + pkin(8) * t148;
t132 = t87 * t148;
t126 = -pkin(1) * t89 + pkin(8) * t128;
t66 = t87 * t128;
t124 = -t107 * t90 + t66;
t122 = -t21 * pkin(4) + t22 * qJ(5);
t121 = -t25 * pkin(4) + qJ(5) * t26;
t120 = -t42 * pkin(4) + qJ(5) * t43;
t117 = t27 * pkin(3);
t116 = t182 * pkin(3);
t102 = t107 * t87 + t128 * t90;
t97 = t102 * pkin(3);
t96 = t119 / 0.2e1 + t114;
t95 = -mrSges(5,1) + t191;
t92 = -pkin(10) - pkin(9);
t88 = sin(qJ(2));
t75 = pkin(3) * t90 + pkin(2);
t58 = t113 + t118 / 0.2e1;
t51 = t154 * t88 + t89 * t96;
t48 = -t154 * t96 + t88 * t89;
t28 = t106 * t90 + t132;
t2 = t26 * t77 + t51 * t76;
t1 = -t26 * t76 + t51 * t77;
t3 = [(-t154 * mrSges(2,1) - m(3) * t136 - t106 * mrSges(3,1) - m(4) * (pkin(2) * t106 + t136) - t28 * mrSges(4,1) - t27 * mrSges(4,2) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t84 + mrSges(2,2)) * t89 + t95 * t26 + t165 * t51 - t195 * t25 - t171 * (pkin(3) * t132 + t106 * t75 - t51 * t92 + t136)) * g(2) + (t89 * mrSges(2,1) + t154 * mrSges(2,2) - m(3) * t126 + t107 * mrSges(3,1) - mrSges(3,3) * t128 - m(4) * (-pkin(2) * t107 + t126) - t124 * mrSges(4,1) - t102 * mrSges(4,2) - (t95 + t172) * t22 + t189 * t48 + t195 * t21 + t171 * (-pkin(3) * t66 + t107 * t75 - t48 * t92 - t126)) * g(1) (-t171 * (t58 * t75 + t60 * t92) + t189 * t60 - t164 * t58) * g(3) + (-t171 * (-t107 * t92 - t48 * t75) - t189 * t107 + t164 * t48) * g(2) + (-t171 * (-t106 * t92 - t51 * t75) - t189 * t106 + t164 * t51) * g(1) (-t182 * mrSges(4,1) - (-t135 * t87 + t60 * t90) * mrSges(4,2) - m(5) * t116 - m(6) * (t116 + t120) - m(7) * (t116 + t143) + t178) * g(3) + (mrSges(4,1) * t102 - mrSges(4,2) * t124 + m(5) * t97 - m(6) * (t122 - t97) - m(7) * (-t97 + t161) + t180) * g(2) + (-mrSges(4,1) * t27 + mrSges(4,2) * t28 - m(5) * t117 - m(6) * (t117 + t121) - m(7) * (t117 + t155) + t179) * g(1) (-m(6) * t120 - m(7) * t143 + t178) * g(3) + (-m(6) * t122 - m(7) * t161 + t180) * g(2) + (-m(6) * t121 - m(7) * t155 + t179) * g(1), t177 * (-g(1) * t25 - g(2) * t21 - g(3) * t42) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t22 * t76 + t48 * t77) * mrSges(7,1) + (-t22 * t77 - t48 * t76) * mrSges(7,2)) - g(3) * ((-t43 * t76 - t58 * t77) * mrSges(7,1) + (-t43 * t77 + t58 * t76) * mrSges(7,2))];
taug  = t3(:);
