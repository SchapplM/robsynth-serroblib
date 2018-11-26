% Calculate Gravitation load on the joints for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2018-11-23 17:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:45:27
% EndTime: 2018-11-23 17:45:28
% DurationCPUTime: 1.43s
% Computational Cost: add. (1710->137), mult. (1842->180), div. (0->0), fcn. (1833->16), ass. (0->67)
t118 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t166 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t160 = t118 * t89 + t166 * t85 + mrSges(5,1);
t155 = m(6) + m(7);
t161 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t81 = qJ(3) + pkin(11);
t78 = sin(t81);
t79 = cos(t81);
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t156 = t155 * (pkin(4) * t79 + pkin(10) * t78) + m(4) * pkin(2) + t90 * mrSges(4,1) - t86 * mrSges(4,2) - t161 * t78 + mrSges(3,1) + t160 * t79;
t165 = m(5) + t155;
t150 = cos(qJ(1));
t131 = pkin(6) + qJ(2);
t111 = sin(t131) / 0.2e1;
t132 = pkin(6) - qJ(2);
t119 = sin(t132);
t62 = t111 - t119 / 0.2e1;
t88 = sin(qJ(1));
t91 = cos(qJ(2));
t102 = t150 * t91 - t88 * t62;
t82 = sin(pkin(6));
t140 = t82 * t88;
t25 = -t102 * t86 + t90 * t140;
t112 = cos(t131) / 0.2e1;
t120 = cos(t132);
t63 = t112 - t120 / 0.2e1;
t83 = cos(pkin(6));
t163 = t63 * t86 + t83 * t90;
t103 = t150 * t62 + t88 * t91;
t127 = t82 * t150;
t18 = t103 * t79 - t78 * t127;
t87 = sin(qJ(2));
t93 = t120 / 0.2e1 + t112;
t47 = -t150 * t93 + t87 * t88;
t1 = t18 * t85 - t47 * t89;
t159 = t18 * t89 + t47 * t85;
t154 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t158 = t118 * t85 - t166 * t89 - t154;
t97 = t103 * t86 + t90 * t127;
t94 = t97 * pkin(3);
t133 = t150 * pkin(1) + pkin(8) * t140;
t130 = t86 * t140;
t126 = -pkin(1) * t88 + pkin(8) * t127;
t17 = -t103 * t78 - t79 * t127;
t70 = t86 * t127;
t124 = -t103 * t90 + t70;
t116 = t25 * pkin(3);
t115 = t163 * pkin(3);
t50 = t150 * t87 + t88 * t93;
t77 = pkin(3) * t90 + pkin(2);
t84 = -qJ(4) - pkin(9);
t113 = pkin(3) * t130 + t102 * t77 - t50 * t84 + t133;
t104 = pkin(3) * t70 - t103 * t77 + t47 * t84 + t126;
t101 = -t155 * pkin(10) + t161;
t61 = t111 + t119 / 0.2e1;
t37 = -t63 * t79 + t78 * t83;
t36 = t63 * t78 + t79 * t83;
t26 = t102 * t90 + t130;
t22 = t102 * t79 + t140 * t78;
t21 = t102 * t78 - t140 * t79;
t11 = t37 * t85 + t61 * t89;
t6 = t22 * t89 + t50 * t85;
t5 = t22 * t85 - t50 * t89;
t2 = [(-t150 * mrSges(2,1) - m(3) * t133 - t102 * mrSges(3,1) - m(4) * (pkin(2) * t102 + t133) - t26 * mrSges(4,1) - t25 * mrSges(4,2) - m(5) * t113 - t22 * mrSges(5,1) + (-mrSges(3,3) * t82 + mrSges(2,2)) * t88 - t118 * t6 - t166 * t5 + t154 * t50 + t101 * t21 - t155 * (t22 * pkin(4) + t113)) * g(2) + (t88 * mrSges(2,1) + t150 * mrSges(2,2) - m(3) * t126 + t103 * mrSges(3,1) - mrSges(3,3) * t127 - m(4) * (-pkin(2) * t103 + t126) - t124 * mrSges(4,1) - t97 * mrSges(4,2) - m(5) * t104 + t18 * mrSges(5,1) + t118 * t159 + t166 * t1 - t154 * t47 + t101 * t17 + t155 * (pkin(4) * t18 - t104)) * g(1) (-t165 * (t61 * t77 + t63 * t84) + t158 * t63 - t156 * t61) * g(3) + (-t165 * (-t103 * t84 - t47 * t77) - t158 * t103 + t156 * t47) * g(2) + (-t165 * (-t102 * t84 - t50 * t77) - t158 * t102 + t156 * t50) * g(1) (-t163 * mrSges(4,1) - (t63 * t90 - t83 * t86) * mrSges(4,2) - m(5) * t115 - t155 * (t36 * pkin(4) + pkin(10) * t37 + t115) + t161 * t37 - t160 * t36) * g(3) + (m(5) * t94 + mrSges(4,1) * t97 - mrSges(4,2) * t124 - t160 * t17 + t161 * t18 - t155 * (t17 * pkin(4) + t18 * pkin(10) - t94)) * g(2) + (-m(5) * t116 - mrSges(4,1) * t25 + mrSges(4,2) * t26 - t155 * (-t21 * pkin(4) + pkin(10) * t22 + t116) + t161 * t22 + t160 * t21) * g(1), t165 * (-g(1) * t50 - g(2) * t47 + g(3) * t61) (-t166 * (t37 * t89 - t61 * t85) + t118 * t11) * g(3) + (t118 * t1 - t159 * t166) * g(2) + (t118 * t5 - t166 * t6) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
