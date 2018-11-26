% Calculate Gravitation load on the joints for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:30:28
% EndTime: 2018-11-23 16:30:30
% DurationCPUTime: 1.71s
% Computational Cost: add. (3906->155), mult. (4181->209), div. (0->0), fcn. (4170->22), ass. (0->98)
t149 = cos(qJ(4));
t135 = pkin(7) - qJ(3);
t123 = cos(t135);
t115 = t123 / 0.2e1;
t134 = pkin(7) + qJ(3);
t122 = cos(t134);
t109 = t115 - t122 / 0.2e1;
t138 = sin(pkin(6));
t102 = t109 * t138;
t150 = cos(qJ(1));
t113 = sin(t134) / 0.2e1;
t121 = sin(t135);
t159 = t113 - t121 / 0.2e1;
t148 = sin(qJ(1));
t132 = pkin(6) + pkin(12);
t110 = sin(t132) / 0.2e1;
t133 = pkin(6) - pkin(12);
t119 = sin(t133);
t58 = t110 - t119 / 0.2e1;
t69 = cos(pkin(12));
t52 = t148 * t69 + t150 * t58;
t75 = cos(qJ(3));
t136 = sin(pkin(12));
t111 = cos(t133) / 0.2e1;
t120 = cos(t132);
t96 = t111 + t120 / 0.2e1;
t90 = t148 * t136 - t150 * t96;
t163 = t159 * t90 - t52 * t75;
t31 = t102 * t150 + t163;
t72 = sin(qJ(4));
t139 = cos(pkin(7));
t112 = t139 * t138;
t137 = sin(pkin(7));
t83 = -t150 * t112 + t90 * t137;
t16 = -t149 * t31 + t72 * t83;
t114 = t122 / 0.2e1;
t101 = t115 + t114;
t73 = sin(qJ(3));
t98 = t113 + t121 / 0.2e1;
t92 = t98 * t138;
t27 = t101 * t90 + t150 * t92 + t52 * t73;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t175 = t16 * t71 - t27 * t74;
t174 = t16 * t74 + t27 * t71;
t173 = t149 * t83 + t31 * t72;
t169 = mrSges(6,1) + mrSges(7,1);
t160 = mrSges(6,2) + mrSges(7,2);
t67 = pkin(5) * t74 + pkin(4);
t167 = m(6) * pkin(4) + m(7) * t67 + mrSges(5,1);
t166 = m(6) * pkin(11) - m(7) * (-qJ(6) - pkin(11)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t87 = t136 * t150 + t148 * t96;
t76 = t148 * t112 + t87 * t137;
t59 = t111 - t120 / 0.2e1;
t95 = t110 + t119 / 0.2e1;
t165 = t159 * t95 + t59 * t75;
t53 = -t148 * t58 + t150 * t69;
t164 = -t159 * t87 + t53 * t75;
t151 = m(7) * pkin(5);
t161 = mrSges(4,2) - mrSges(5,3);
t158 = -m(5) - m(6) - m(7);
t157 = -t151 - t169;
t156 = -t160 * t71 + t169 * t74 + t167;
t155 = t167 * t149 + t166 * t72 + mrSges(4,1);
t153 = -m(7) * (t71 * pkin(5) + pkin(10)) + t161;
t100 = t114 - t123 / 0.2e1;
t93 = t100 * t138;
t29 = t150 * t93 - t163;
t147 = t29 * t71;
t34 = -t148 * t93 + t164;
t146 = t34 * t71;
t140 = cos(pkin(6));
t38 = -t100 * t140 + t165;
t145 = t38 * t71;
t116 = t138 * t148;
t141 = t150 * pkin(1) + qJ(2) * t116;
t129 = t71 * t149;
t128 = t74 * t149;
t117 = t150 * t138;
t118 = -pkin(1) * t148 + qJ(2) * t117;
t33 = t102 * t148 + t164;
t20 = t149 * t33 + t72 * t76;
t32 = t101 * t87 - t148 * t92 + t53 * t73;
t5 = -t20 * t71 + t32 * t74;
t86 = -t137 * t95 + t139 * t140;
t84 = -t52 * pkin(2) - t83 * pkin(9) + t118;
t81 = t31 * pkin(3) + t84;
t80 = t53 * pkin(2) + t76 * pkin(9) + t141;
t79 = t33 * pkin(3) + t80;
t78 = -pkin(10) * t27 + t81;
t77 = t32 * pkin(10) + t79;
t37 = t109 * t140 + t165;
t36 = -t101 * t95 - t140 * t98 + t59 * t73;
t22 = t149 * t37 + t72 * t86;
t21 = -t149 * t86 + t37 * t72;
t19 = -t149 * t76 + t33 * t72;
t6 = t20 * t74 + t32 * t71;
t1 = [(-t150 * mrSges(2,1) + t148 * mrSges(2,2) - m(3) * t141 - t53 * mrSges(3,1) + t87 * mrSges(3,2) - mrSges(3,3) * t116 - m(4) * t80 - t33 * mrSges(4,1) - t76 * mrSges(4,3) - m(5) * t77 - t20 * mrSges(5,1) - m(6) * (t20 * pkin(4) + t77) - m(7) * (t20 * t67 + t79) - t169 * t6 - t160 * t5 + t153 * t32 - t166 * t19) * g(2) + (t148 * mrSges(2,1) + t150 * mrSges(2,2) - m(3) * t118 + t52 * mrSges(3,1) - t90 * mrSges(3,2) - mrSges(3,3) * t117 - m(4) * t84 - t31 * mrSges(4,1) + t83 * mrSges(4,3) - m(5) * t78 + t16 * mrSges(5,1) - m(6) * (-pkin(4) * t16 + t78) - m(7) * (-t16 * t67 + t81) + t169 * t174 - t153 * t27 - t160 * t175 - t166 * t173) * g(1) (-g(1) * t116 + g(2) * t117 - g(3) * t140) * (m(3) + m(4) - t158) (-t145 * t151 + t161 * t38 - t169 * (-t128 * t36 + t145) - t160 * (t129 * t36 + t38 * t74) + t158 * (-t36 * pkin(3) + t38 * pkin(10)) + t155 * t36) * g(3) + (-t147 * t151 - t169 * (-t128 * t27 + t147) - t160 * (t129 * t27 + t29 * t74) + t161 * t29 + t158 * (-t27 * pkin(3) + t29 * pkin(10)) + t155 * t27) * g(2) + (-t146 * t151 - t160 * (t129 * t32 + t34 * t74) + t161 * t34 + t158 * (-t32 * pkin(3) + t34 * pkin(10)) - t169 * (-t128 * t32 + t146) + t155 * t32) * g(1) (t156 * t21 - t166 * t22) * g(3) + (-t156 * t173 - t16 * t166) * g(2) + (t156 * t19 - t166 * t20) * g(1) (-t160 * (-t22 * t74 - t36 * t71) + t157 * (-t22 * t71 + t36 * t74)) * g(3) + (-t157 * t175 + t160 * t174) * g(2) + (t157 * t5 + t160 * t6) * g(1) (-g(1) * t19 + g(2) * t173 - g(3) * t21) * m(7)];
taug  = t1(:);
