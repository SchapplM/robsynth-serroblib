% Calculate Gravitation load on the joints for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:28
% EndTime: 2019-03-10 01:47:32
% DurationCPUTime: 1.54s
% Computational Cost: add. (1057->172), mult. (1810->232), div. (0->0), fcn. (2137->12), ass. (0->91)
t101 = sin(qJ(5));
t196 = mrSges(6,2) - mrSges(7,3);
t198 = t101 * t196 - mrSges(5,1);
t189 = mrSges(6,3) + mrSges(7,2);
t195 = mrSges(5,2) - t189;
t197 = mrSges(6,1) + mrSges(7,1);
t183 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t102 = sin(qJ(3));
t106 = cos(qJ(3));
t100 = sin(pkin(6));
t103 = sin(qJ(2));
t154 = t100 * t103;
t156 = cos(pkin(6));
t193 = -t102 * t154 + t156 * t106;
t104 = sin(qJ(1));
t152 = t100 * t106;
t107 = cos(qJ(2));
t134 = t104 * t156;
t177 = cos(qJ(1));
t77 = -t103 * t134 + t107 * t177;
t41 = -t102 * t77 + t104 * t152;
t99 = qJ(3) + qJ(4);
t96 = sin(t99);
t97 = cos(t99);
t192 = -m(4) * pkin(2) - mrSges(4,1) * t106 - mrSges(5,1) * t97 + mrSges(4,2) * t102 + mrSges(5,2) * t96 - mrSges(3,1);
t105 = cos(qJ(5));
t126 = t156 * t177;
t74 = t103 * t104 - t107 * t126;
t158 = t74 * t105;
t137 = t100 * t177;
t75 = t103 * t126 + t104 * t107;
t36 = -t96 * t137 + t75 * t97;
t1 = t101 * t36 - t158;
t191 = t101 * t74 + t105 * t36;
t190 = m(6) + m(7);
t155 = qJ(6) * t101;
t35 = -t97 * t137 - t75 * t96;
t161 = t105 * t35;
t188 = pkin(5) * t161 + t35 * t155;
t153 = t100 * t104;
t39 = -t153 * t97 + t77 * t96;
t160 = t105 * t39;
t187 = -pkin(5) * t160 - t39 * t155;
t63 = -t154 * t96 + t156 * t97;
t159 = t105 * t63;
t186 = pkin(5) * t159 + t63 * t155;
t64 = t154 * t97 + t156 * t96;
t181 = -t197 * t159 + t195 * t64 + t198 * t63;
t40 = t153 * t96 + t77 * t97;
t180 = t197 * t160 + t195 * t40 - t198 * t39;
t179 = -t197 * t161 + t195 * t36 + t198 * t35;
t130 = m(7) * pkin(5) + t197;
t128 = -m(7) * qJ(6) + t196;
t178 = pkin(4) * t97;
t175 = t74 * t96;
t76 = t103 * t177 + t107 * t134;
t174 = t76 * t96;
t108 = -pkin(10) - pkin(9);
t95 = pkin(3) * t106 + pkin(2);
t171 = -t75 * t108 - t74 * t95;
t170 = -t77 * t108 - t76 * t95;
t168 = t177 * pkin(1) + pkin(8) * t153;
t163 = t101 * t97;
t157 = t76 * t105;
t151 = t100 * t107;
t150 = t105 * t107;
t146 = t96 * t151;
t144 = t102 * t153;
t141 = t101 * t151;
t140 = t35 * pkin(4) + t36 * pkin(11);
t139 = -t39 * pkin(4) + pkin(11) * t40;
t138 = t63 * pkin(4) + pkin(11) * t64;
t136 = -pkin(1) * t104 + pkin(8) * t137;
t90 = t102 * t137;
t135 = -t106 * t75 + t90;
t129 = t41 * pkin(3);
t127 = pkin(3) * t144 - t76 * t108 + t77 * t95 + t168;
t119 = t193 * pkin(3);
t118 = pkin(3) * t90 + t74 * t108 - t75 * t95 + t136;
t117 = -t190 * pkin(11) + t195;
t113 = t129 + t139;
t112 = t75 * t102 + t106 * t137;
t111 = t119 + t138;
t110 = t112 * pkin(3);
t109 = -t110 + t140;
t78 = t95 * t151;
t42 = t106 * t77 + t144;
t33 = t100 * t150 + t101 * t64;
t6 = t101 * t76 + t105 * t40;
t5 = t101 * t40 - t157;
t2 = [(-t177 * mrSges(2,1) - m(3) * t168 - t77 * mrSges(3,1) - m(4) * (pkin(2) * t77 + t168) - t42 * mrSges(4,1) - t41 * mrSges(4,2) - m(5) * t127 - t40 * mrSges(5,1) - t130 * t6 + t128 * t5 + (-mrSges(3,3) * t100 + mrSges(2,2)) * t104 + t183 * t76 + t117 * t39 - t190 * (t40 * pkin(4) + t127)) * g(2) + (t104 * mrSges(2,1) + t177 * mrSges(2,2) - m(3) * t136 + t75 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * (-pkin(2) * t75 + t136) - t135 * mrSges(4,1) - t112 * mrSges(4,2) - m(5) * t118 + t36 * mrSges(5,1) + t130 * t191 - t128 * t1 - t183 * t74 + t117 * t35 + t190 * (pkin(4) * t36 - t118)) * g(1) (-m(5) * t171 - t190 * (-pkin(11) * t175 - t74 * t178 + t171) - t130 * (t101 * t75 - t158 * t97) + t128 * (-t75 * t105 - t163 * t74) + t189 * t175 + t183 * t75 - t192 * t74) * g(2) + (-m(5) * t170 - t190 * (-pkin(11) * t174 - t76 * t178 + t170) + t128 * (-t77 * t105 - t163 * t76) + t189 * t174 - t130 * (t101 * t77 - t157 * t97) + t183 * t77 - t192 * t76) * g(1) + (-m(5) * t78 - t190 * (pkin(11) * t146 - t108 * t154 + t151 * t178 + t78) + t128 * (-t105 * t154 + t141 * t97) - t189 * t146 + (-t130 * t150 * t97 + t192 * t107 + (m(5) * t108 - t130 * t101 + t183) * t103) * t100) * g(3) (-t193 * mrSges(4,1) - (-t102 * t156 - t103 * t152) * mrSges(4,2) - m(5) * t119 - m(6) * t111 - m(7) * (t111 + t186) + t181) * g(3) + (mrSges(4,1) * t112 - mrSges(4,2) * t135 + m(5) * t110 - m(6) * t109 - m(7) * (t109 + t188) + t179) * g(2) + (-mrSges(4,1) * t41 + mrSges(4,2) * t42 - m(5) * t129 - m(6) * t113 - m(7) * (t113 + t187) + t180) * g(1) (-m(6) * t138 - m(7) * (t138 + t186) + t181) * g(3) + (-m(6) * t140 - m(7) * (t140 + t188) + t179) * g(2) + (-m(6) * t139 - m(7) * (t139 + t187) + t180) * g(1) (t128 * (t105 * t64 - t141) + t130 * t33) * g(3) + (t130 * t1 + t128 * t191) * g(2) + (t128 * t6 + t130 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t33) * m(7)];
taug  = t2(:);
