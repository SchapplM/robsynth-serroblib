% Calculate Gravitation load on the joints for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:14
% EndTime: 2019-03-10 00:34:19
% DurationCPUTime: 2.14s
% Computational Cost: add. (1448->193), mult. (3929->267), div. (0->0), fcn. (4968->14), ass. (0->111)
t150 = cos(pkin(6));
t163 = cos(qJ(2));
t134 = t150 * t163;
t160 = sin(qJ(2));
t161 = sin(qJ(1));
t89 = cos(qJ(1));
t115 = -t134 * t89 + t161 * t160;
t147 = sin(pkin(7));
t148 = sin(pkin(6));
t124 = t148 * t147;
t149 = cos(pkin(7));
t193 = t115 * t149 + t89 * t124;
t185 = mrSges(5,2) - mrSges(6,3);
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t192 = t84 * mrSges(7,1) + t87 * mrSges(7,2) - t185;
t162 = cos(qJ(3));
t133 = t150 * t160;
t68 = t133 * t89 + t161 * t163;
t86 = sin(qJ(3));
t33 = -t162 * t68 + t193 * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t125 = t149 * t148;
t98 = t115 * t147 - t89 * t125;
t10 = t33 * t88 - t85 * t98;
t9 = t33 * t85 + t88 * t98;
t186 = m(6) + m(7);
t190 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t102 = t134 * t161 + t160 * t89;
t91 = t102 * t147 + t161 * t125;
t175 = m(7) * pkin(12) + t190;
t187 = pkin(4) * t186 + t175;
t30 = t193 * t162 + t68 * t86;
t151 = qJ(5) * t85;
t159 = t30 * t88;
t184 = -pkin(4) * t159 - t30 * t151;
t179 = t102 * t149 - t161 * t124;
t69 = -t133 * t161 + t163 * t89;
t34 = t179 * t162 + t69 * t86;
t158 = t34 * t88;
t183 = -pkin(4) * t158 - t34 * t151;
t129 = t148 * t160;
t178 = t163 * t125 + t150 * t147;
t51 = t129 * t86 - t178 * t162;
t157 = t51 * t88;
t182 = -pkin(4) * t157 - t51 * t151;
t181 = -t147 * mrSges(4,3) + mrSges(3,2);
t174 = -m(7) * (pkin(5) + pkin(11)) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t173 = -t186 * qJ(5) - t192;
t171 = t190 * t88 + t192 * t85 + mrSges(4,1);
t170 = -t87 * mrSges(7,1) + t84 * mrSges(7,2) + t174;
t113 = t160 * t125;
t131 = t163 * t148;
t61 = t113 * t162 + t131 * t86;
t168 = pkin(11) * t61;
t167 = t30 * pkin(11);
t166 = t34 * pkin(11);
t132 = t149 * t162;
t38 = -t115 * t86 + t132 * t68;
t165 = t38 * pkin(11);
t40 = -t102 * t86 + t132 * t69;
t164 = t40 * pkin(11);
t111 = t160 * t124;
t153 = pkin(2) * t131 + pkin(10) * t111;
t130 = t148 * t161;
t152 = t89 * pkin(1) + pkin(9) * t130;
t62 = -t113 * t86 + t131 * t162;
t146 = t62 * pkin(3) + t153;
t24 = t30 * pkin(3);
t145 = -pkin(11) * t33 - t24;
t26 = t34 * pkin(3);
t35 = t69 * t162 - t179 * t86;
t144 = pkin(11) * t35 - t26;
t50 = t51 * pkin(3);
t52 = t162 * t129 + t178 * t86;
t143 = pkin(11) * t52 - t50;
t142 = pkin(10) * t147;
t141 = t85 * t147;
t140 = t86 * t149;
t139 = t88 * t147;
t138 = t89 * t148;
t136 = -pkin(1) * t161 + pkin(9) * t138;
t123 = -t115 * pkin(2) + t142 * t68;
t122 = -t102 * pkin(2) + t142 * t69;
t43 = -t111 * t88 + t62 * t85;
t44 = t111 * t85 + t62 * t88;
t119 = t44 * pkin(4) + qJ(5) * t43 + t146;
t39 = -t115 * t162 - t140 * t68;
t118 = t39 * pkin(3) + t123;
t41 = -t102 * t162 - t140 * t69;
t117 = t41 * pkin(3) + t122;
t15 = -t139 * t68 + t39 * t85;
t16 = t141 * t68 + t39 * t88;
t108 = t16 * pkin(4) + t15 * qJ(5) + t118;
t17 = -t139 * t69 + t41 * t85;
t18 = t141 * t69 + t41 * t88;
t107 = t18 * pkin(4) + t17 * qJ(5) + t117;
t101 = -t124 * t163 + t149 * t150;
t96 = -t68 * pkin(2) - t98 * pkin(10) + t136;
t95 = t33 * pkin(3) + t96;
t94 = t69 * pkin(2) + t91 * pkin(10) + t152;
t93 = t35 * pkin(3) + t94;
t92 = t10 * pkin(4) + t9 * qJ(5) + t95;
t11 = t35 * t85 - t88 * t91;
t12 = t35 * t88 + t85 * t91;
t90 = t12 * pkin(4) + t11 * qJ(5) + t93;
t28 = -t101 * t88 + t52 * t85;
t2 = t11 * t84 + t34 * t87;
t1 = t11 * t87 - t34 * t84;
t3 = [(-t89 * mrSges(2,1) + t161 * mrSges(2,2) - m(3) * t152 - t69 * mrSges(3,1) + t102 * mrSges(3,2) - mrSges(3,3) * t130 - m(4) * t94 - t35 * mrSges(4,1) - t91 * mrSges(4,3) - m(5) * (t93 + t166) - m(6) * (t90 + t166) - m(7) * t90 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t185 * t11 + t174 * t34 - t175 * t12) * g(2) + (t161 * mrSges(2,1) + t89 * mrSges(2,2) - m(3) * t136 + t68 * mrSges(3,1) - t115 * mrSges(3,2) - mrSges(3,3) * t138 - m(4) * t96 - t33 * mrSges(4,1) + t98 * mrSges(4,3) - m(5) * (t95 - t167) - m(6) * (t92 - t167) - m(7) * t92 - t192 * t9 - t170 * t30 - t175 * t10) * g(1) (-mrSges(3,1) * t131 + mrSges(3,2) * t129 - m(4) * t153 - t62 * mrSges(4,1) - mrSges(4,3) * t111 - m(5) * (t146 + t168) - m(6) * (t119 + t168) - m(7) * t119 - t192 * t43 + t170 * t61 - t175 * t44) * g(3) + (t115 * mrSges(3,1) - m(4) * t123 - t39 * mrSges(4,1) - m(5) * (t118 + t165) - m(6) * (t108 + t165) - m(7) * t108 + t181 * t68 - t192 * t15 + t170 * t38 - t175 * t16) * g(2) + (mrSges(3,1) * t102 - m(4) * t122 - t41 * mrSges(4,1) - m(5) * (t117 + t164) - m(6) * (t107 + t164) - m(7) * t107 + t181 * t69 - t192 * t17 + t170 * t40 - t175 * t18) * g(1) (-m(5) * t143 - m(6) * (t143 + t182) - m(7) * (-pkin(12) * t157 + t182 - t50) + t170 * t52 + t171 * t51) * g(3) + (-m(5) * t145 - m(6) * (t145 + t184) - m(7) * (-pkin(12) * t159 + t184 - t24) - t170 * t33 + t171 * t30) * g(2) + (-m(5) * t144 - m(6) * (t144 + t183) - m(7) * (-pkin(12) * t158 + t183 - t26) + t170 * t35 + t171 * t34) * g(1) (t173 * (t101 * t85 + t52 * t88) + t187 * t28) * g(3) + (-t173 * t10 - t187 * t9) * g(2) + (t187 * t11 + t173 * t12) * g(1), t186 * (-g(1) * t11 + g(2) * t9 - g(3) * t28) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t30 * t84 - t87 * t9) * mrSges(7,1) + (-t30 * t87 + t84 * t9) * mrSges(7,2)) - g(3) * ((t28 * t87 - t51 * t84) * mrSges(7,1) + (-t28 * t84 - t51 * t87) * mrSges(7,2))];
taug  = t3(:);
