% Calculate Gravitation load on the joints for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:59
% EndTime: 2019-03-09 19:52:04
% DurationCPUTime: 1.92s
% Computational Cost: add. (1374->185), mult. (3274->271), div. (0->0), fcn. (4099->16), ass. (0->84)
t90 = sin(qJ(6));
t92 = cos(qJ(6));
t166 = m(7) * pkin(5) + t92 * mrSges(7,1) - t90 * mrSges(7,2) + mrSges(6,1);
t165 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t85 = sin(pkin(13));
t88 = cos(pkin(13));
t158 = -m(5) * pkin(3) - t88 * mrSges(5,1) + t85 * mrSges(5,2) - mrSges(4,1);
t137 = cos(pkin(6));
t154 = cos(qJ(2));
t117 = t137 * t154;
t151 = sin(qJ(2));
t152 = sin(qJ(1));
t93 = cos(qJ(1));
t102 = t117 * t152 + t151 * t93;
t136 = cos(pkin(7));
t87 = sin(pkin(6));
t128 = t87 * t136;
t86 = sin(pkin(7));
t171 = t102 * t86 + t152 * t128;
t63 = -t117 * t93 + t151 * t152;
t170 = t93 * t128 - t63 * t86;
t91 = sin(qJ(3));
t127 = t91 * t136;
t143 = t87 * t93;
t135 = t86 * t143;
t153 = cos(qJ(3));
t116 = t137 * t151;
t64 = t116 * t93 + t152 * t154;
t24 = -t127 * t63 - t91 * t135 + t153 * t64;
t84 = pkin(13) + qJ(5);
t81 = sin(t84);
t82 = cos(t84);
t4 = -t170 * t81 + t24 * t82;
t169 = -t170 * t82 - t24 * t81;
t168 = -m(4) - m(5);
t167 = -m(6) - m(7);
t132 = t87 * t152;
t164 = -t102 * t136 + t86 * t132;
t110 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t163 = -t90 * mrSges(7,1) - t92 * mrSges(7,2) + t110;
t162 = -t165 * t81 + t166 * t82 - t158;
t159 = mrSges(3,2) + (-t85 * mrSges(5,1) - t88 * mrSges(5,2) - mrSges(4,3)) * t86;
t155 = pkin(10) * t86;
t150 = t170 * t85;
t148 = t81 * t86;
t147 = t82 * t86;
t131 = t87 * t151;
t125 = t86 * t131;
t133 = t87 * t154;
t139 = pkin(2) * t133 + pkin(10) * t125;
t138 = t93 * pkin(1) + pkin(9) * t132;
t129 = t86 * t137;
t120 = (pkin(4) * t85 + pkin(10)) * t86;
t119 = -pkin(1) * t152 + pkin(9) * t143;
t115 = t136 * t153;
t114 = t136 * t151;
t111 = t85 * t125;
t103 = -t64 * pkin(2) + pkin(10) * t170 + t119;
t23 = t115 * t63 + t135 * t153 + t64 * t91;
t65 = -t116 * t152 + t154 * t93;
t97 = t65 * pkin(2) + pkin(10) * t171 + t138;
t27 = -t164 * t153 + t65 * t91;
t28 = t65 * t153 + t164 * t91;
t79 = pkin(4) * t88 + pkin(3);
t89 = -pkin(11) - qJ(4);
t94 = t171 * t85;
t95 = pkin(4) * t94 - t27 * t89 + t28 * t79 + t97;
t62 = -t133 * t86 + t136 * t137;
t60 = t102 * pkin(2);
t58 = t63 * pkin(2);
t57 = (-t114 * t91 + t153 * t154) * t87;
t56 = (t114 * t153 + t154 * t91) * t87;
t46 = t91 * t129 + (t127 * t154 + t151 * t153) * t87;
t45 = -t115 * t133 - t129 * t153 + t131 * t91;
t38 = -t102 * t153 - t127 * t65;
t37 = -t102 * t91 + t115 * t65;
t36 = -t127 * t64 - t153 * t63;
t35 = t115 * t64 - t63 * t91;
t18 = t46 * t82 + t62 * t81;
t8 = t171 * t81 + t28 * t82;
t7 = -t171 * t82 + t28 * t81;
t2 = t27 * t90 + t8 * t92;
t1 = t27 * t92 - t8 * t90;
t3 = [(-t93 * mrSges(2,1) + t152 * mrSges(2,2) - m(3) * t138 - t65 * mrSges(3,1) + t102 * mrSges(3,2) - mrSges(3,3) * t132 - m(4) * t97 - t28 * mrSges(4,1) - t171 * mrSges(4,3) - m(5) * (t28 * pkin(3) + t97) - (t28 * t88 + t94) * mrSges(5,1) - (t171 * t88 - t28 * t85) * mrSges(5,2) - m(6) * t95 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t95) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t165 * t7 + t110 * t27) * g(2) + (t152 * mrSges(2,1) - m(3) * t119 + t64 * mrSges(3,1) - t63 * mrSges(3,2) - m(4) * t103 + t24 * mrSges(4,1) - t170 * mrSges(4,3) - m(5) * (-pkin(3) * t24 + t103) - (-t24 * t88 + t150) * mrSges(5,1) - (t170 * t88 + t24 * t85) * mrSges(5,2) + (-mrSges(3,3) * t87 + mrSges(2,2)) * t93 + t165 * t169 + t166 * t4 - t163 * t23 + t167 * (pkin(4) * t150 + t23 * t89 - t24 * t79 + t103)) * g(1) (-(mrSges(3,1) * t154 - mrSges(3,2) * t151) * t87 - m(4) * t139 - t57 * mrSges(4,1) - mrSges(4,3) * t125 - m(5) * (pkin(3) * t57 + t139) - (t57 * t88 + t111) * mrSges(5,1) - (t125 * t88 - t57 * t85) * mrSges(5,2) + t167 * (pkin(4) * t111 - t56 * t89 + t57 * t79 + t139) + t165 * (-t125 * t82 + t57 * t81) - t166 * (t125 * t81 + t57 * t82) + t163 * t56) * g(3) + (mrSges(3,1) * t63 + t167 * (t120 * t64 - t35 * t89 + t36 * t79 - t58) + t165 * (-t147 * t64 + t36 * t81) + t168 * (t155 * t64 - t58) + t158 * t36 + t159 * t64 - t166 * (t148 * t64 + t36 * t82) + t163 * t35) * g(2) + (t102 * mrSges(3,1) + t168 * (t155 * t65 - t60) + t158 * t38 + t159 * t65 + t167 * (t120 * t65 - t37 * t89 + t38 * t79 - t60) + t165 * (-t147 * t65 + t38 * t81) - t166 * (t148 * t65 + t38 * t82) + t163 * t37) * g(1) (t167 * (-t45 * t79 - t46 * t89) + t163 * t46 + t162 * t45) * g(3) + (t167 * (-t23 * t79 - t24 * t89) + t163 * t24 + t162 * t23) * g(2) + (t167 * (-t27 * t79 - t28 * t89) + t163 * t28 + t162 * t27) * g(1) (m(5) - t167) * (-g(1) * t27 - g(2) * t23 - g(3) * t45) (t165 * t18 - t166 * (-t46 * t81 + t62 * t82)) * g(3) + (t165 * t4 - t166 * t169) * g(2) + (t165 * t8 + t166 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t23 * t92 - t4 * t90) * mrSges(7,1) + (-t23 * t90 - t4 * t92) * mrSges(7,2)) - g(3) * ((-t18 * t90 + t45 * t92) * mrSges(7,1) + (-t18 * t92 - t45 * t90) * mrSges(7,2))];
taug  = t3(:);
