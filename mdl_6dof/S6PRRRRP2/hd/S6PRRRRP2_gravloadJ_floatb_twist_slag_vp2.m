% Calculate Gravitation load on the joints for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:28:15
% EndTime: 2018-11-23 15:28:16
% DurationCPUTime: 0.95s
% Computational Cost: add. (1524->120), mult. (1595->156), div. (0->0), fcn. (1583->16), ass. (0->73)
t174 = mrSges(6,2) - mrSges(7,3);
t91 = sin(qJ(5));
t176 = t91 * t174 - mrSges(5,1);
t175 = -mrSges(6,1) - mrSges(7,1);
t173 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t112 = -m(7) * qJ(6) + t174;
t115 = m(7) * pkin(5) - t175;
t163 = -m(6) - m(7);
t87 = qJ(3) + qJ(4);
t85 = sin(t87);
t86 = cos(t87);
t92 = sin(qJ(3));
t94 = cos(qJ(5));
t95 = cos(qJ(3));
t171 = t163 * (pkin(4) * t86 + pkin(10) * t85) - m(4) * pkin(2) - t95 * mrSges(4,1) + t92 * mrSges(4,2) + t173 * t85 - mrSges(3,1) - (-t112 * t91 + t115 * t94 + mrSges(5,1)) * t86;
t170 = -m(5) + t163;
t131 = cos(pkin(11));
t129 = pkin(6) + qJ(2);
t81 = sin(t129) / 0.2e1;
t130 = pkin(6) - qJ(2);
t84 = sin(t130);
t76 = t81 - t84 / 0.2e1;
t88 = sin(pkin(11));
t96 = cos(qJ(2));
t105 = t131 * t96 - t88 * t76;
t89 = sin(pkin(6));
t140 = t88 * t89;
t169 = -t105 * t92 + t95 * t140;
t111 = cos(t129) / 0.2e1;
t116 = cos(t130);
t77 = t111 - t116 / 0.2e1;
t90 = cos(pkin(6));
t168 = t77 * t92 + t90 * t95;
t166 = m(4) * pkin(8) + t112 * t94 + t115 * t91 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t132 = qJ(6) * t91;
t106 = t131 * t76 + t88 * t96;
t120 = t89 * t131;
t31 = -t106 * t85 - t120 * t86;
t151 = t31 * t94;
t162 = pkin(5) * t151 + t31 * t132;
t33 = -t105 * t85 + t140 * t86;
t149 = t33 * t94;
t161 = pkin(5) * t149 + t33 * t132;
t53 = t77 * t85 + t86 * t90;
t147 = t53 * t94;
t160 = pkin(5) * t147 + t53 * t132;
t54 = -t77 * t86 + t85 * t90;
t158 = t175 * t147 + t173 * t54 + t176 * t53;
t34 = t105 * t86 + t140 * t85;
t157 = t175 * t149 + t173 * t34 + t176 * t33;
t32 = t106 * t86 - t120 * t85;
t156 = t175 * t151 + t173 * t32 + t176 * t31;
t123 = t31 * pkin(4) + t32 * pkin(10);
t122 = t33 * pkin(4) + pkin(10) * t34;
t121 = t53 * pkin(4) + pkin(10) * t54;
t114 = t169 * pkin(3);
t113 = t168 * pkin(3);
t103 = t114 + t122;
t102 = t113 + t121;
t101 = -t106 * t92 - t120 * t95;
t100 = t101 * pkin(3);
t99 = t116 / 0.2e1 + t111;
t98 = t100 + t123;
t97 = -pkin(9) - pkin(8);
t93 = sin(qJ(2));
t83 = pkin(3) * t95 + pkin(2);
t75 = t81 + t84 / 0.2e1;
t64 = t131 * t93 + t88 * t99;
t61 = -t131 * t99 + t88 * t93;
t9 = t54 * t91 + t75 * t94;
t3 = t34 * t91 - t64 * t94;
t1 = t32 * t91 - t61 * t94;
t2 = [(-m(2) - m(3) - m(4) + t170) * g(3) (t170 * (t75 * t83 + t77 * t97) + t166 * t77 + t171 * t75) * g(3) + (t170 * (-t106 * t97 - t61 * t83) - t166 * t106 - t171 * t61) * g(2) + (t170 * (-t105 * t97 - t64 * t83) - t166 * t105 - t171 * t64) * g(1) (-t168 * mrSges(4,1) - (t77 * t95 - t90 * t92) * mrSges(4,2) - m(5) * t113 - m(6) * t102 - m(7) * (t102 + t160) + t158) * g(3) + (-t101 * mrSges(4,1) - (-t106 * t95 + t120 * t92) * mrSges(4,2) - m(5) * t100 - m(6) * t98 - m(7) * (t98 + t162) + t156) * g(2) + (-t169 * mrSges(4,1) - (-t105 * t95 - t140 * t92) * mrSges(4,2) - m(5) * t114 - m(6) * t103 - m(7) * (t103 + t161) + t157) * g(1) (-m(6) * t121 - m(7) * (t121 + t160) + t158) * g(3) + (-m(6) * t123 - m(7) * (t123 + t162) + t156) * g(2) + (-m(6) * t122 - m(7) * (t122 + t161) + t157) * g(1) (t115 * t9 + t112 * (t54 * t94 - t75 * t91)) * g(3) + (t112 * (t32 * t94 + t61 * t91) + t115 * t1) * g(2) + (t112 * (t34 * t94 + t64 * t91) + t115 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
