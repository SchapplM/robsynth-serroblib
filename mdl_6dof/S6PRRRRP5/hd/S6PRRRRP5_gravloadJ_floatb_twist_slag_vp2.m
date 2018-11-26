% Calculate Gravitation load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:30:31
% EndTime: 2018-11-23 15:30:32
% DurationCPUTime: 1.52s
% Computational Cost: add. (4009->171), mult. (4168->243), div. (0->0), fcn. (4076->22), ass. (0->103)
t169 = mrSges(6,1) + mrSges(7,1);
t145 = -mrSges(6,2) - mrSges(7,2);
t87 = cos(qJ(5));
t81 = pkin(5) * t87 + pkin(4);
t168 = m(6) * pkin(4) + m(7) * t81 + mrSges(5,1);
t167 = m(6) * pkin(11) - m(7) * (-qJ(6) - pkin(11)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t132 = pkin(7) + qJ(3);
t116 = sin(t132) / 0.2e1;
t133 = pkin(7) - qJ(3);
t121 = sin(t133);
t162 = t116 - t121 / 0.2e1;
t137 = sin(pkin(12));
t139 = cos(pkin(12));
t134 = pkin(6) + qJ(2);
t117 = sin(t134) / 0.2e1;
t135 = pkin(6) - qJ(2);
t122 = sin(t135);
t75 = t117 - t122 / 0.2e1;
t90 = cos(qJ(2));
t68 = -t137 * t75 + t139 * t90;
t89 = cos(qJ(3));
t120 = cos(t134) / 0.2e1;
t125 = cos(t135);
t103 = t125 / 0.2e1 + t120;
t157 = sin(qJ(2));
t94 = t137 * t103 + t139 * t157;
t166 = -t162 * t94 + t68 * t89;
t66 = t137 * t90 + t139 * t75;
t93 = -t139 * t103 + t137 * t157;
t165 = -t162 * t93 + t66 * t89;
t100 = t117 + t122 / 0.2e1;
t76 = t120 - t125 / 0.2e1;
t164 = t100 * t162 - t76 * t89;
t158 = m(7) * pkin(5);
t163 = mrSges(4,2) - mrSges(5,3);
t136 = -m(5) - m(6) - m(7);
t161 = -t158 - t169;
t84 = sin(qJ(5));
t160 = t145 * t84 + t169 * t87 + t168;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t159 = t167 * t85 + t168 * t88 + mrSges(4,1);
t123 = cos(t132);
t118 = t123 / 0.2e1;
t124 = cos(t133);
t101 = t118 - t124 / 0.2e1;
t138 = sin(pkin(6));
t97 = t101 * t138;
t33 = t139 * t97 + t165;
t156 = t33 * t84;
t36 = -t137 * t97 + t166;
t155 = t36 * t84;
t141 = cos(pkin(6));
t48 = -t141 * t101 + t164;
t154 = t48 * t84;
t82 = sin(pkin(7));
t150 = t82 * t85;
t149 = t82 * t88;
t148 = t84 * t88;
t147 = t87 * t88;
t42 = -t162 * t66 - t93 * t89;
t64 = t93 * pkin(2);
t144 = t42 * pkin(3) - t64;
t44 = -t162 * t68 - t94 * t89;
t65 = t94 * pkin(2);
t143 = t44 * pkin(3) - t65;
t51 = t100 * t89 + t162 * t76;
t74 = t100 * pkin(2);
t142 = t51 * pkin(3) + t74;
t140 = cos(pkin(7));
t131 = m(4) - t136;
t119 = t124 / 0.2e1;
t114 = t140 * t138;
t113 = t119 - t123 / 0.2e1;
t108 = t113 * t138;
t107 = -mrSges(3,2) + (t131 * pkin(9) + mrSges(4,3)) * t82;
t105 = t136 * pkin(10) - t84 * t158 + t163;
t102 = t119 + t118;
t98 = t116 + t121 / 0.2e1;
t96 = t98 * t138;
t95 = -t100 * t82 + t141 * t140;
t92 = t137 * t114 + t94 * t82;
t91 = -t139 * t114 + t93 * t82;
t86 = sin(qJ(3));
t50 = t100 * t86 - t76 * t102;
t47 = t141 * t113 + t164;
t46 = -t100 * t102 - t141 * t98 - t76 * t86;
t43 = t102 * t68 - t94 * t86;
t41 = t102 * t66 - t93 * t86;
t40 = -t76 * t150 + t51 * t88;
t35 = t137 * t108 + t166;
t34 = t94 * t102 - t137 * t96 + t68 * t86;
t32 = -t139 * t108 + t165;
t31 = t93 * t102 + t139 * t96 + t66 * t86;
t28 = t47 * t88 + t95 * t85;
t27 = t47 * t85 - t95 * t88;
t26 = t150 * t68 + t44 * t88;
t24 = t150 * t66 + t42 * t88;
t18 = t35 * t88 + t92 * t85;
t17 = t35 * t85 - t92 * t88;
t16 = t32 * t88 + t91 * t85;
t15 = t32 * t85 - t91 * t88;
t1 = [(-m(2) - m(3) - t131) * g(3) (-t100 * mrSges(3,1) - m(4) * t74 - t51 * mrSges(4,1) - m(5) * t142 - t40 * mrSges(5,1) - m(6) * (pkin(4) * t40 + t142) - m(7) * (t40 * t81 + t142) + t107 * t76 - t169 * (t40 * t87 + t50 * t84) + t145 * (-t40 * t84 + t50 * t87) + t105 * t50 - t167 * (t76 * t149 + t51 * t85)) * g(3) + (t93 * mrSges(3,1) + m(4) * t64 - t42 * mrSges(4,1) - m(5) * t144 - t24 * mrSges(5,1) - m(6) * (pkin(4) * t24 + t144) - m(7) * (t24 * t81 + t144) + t145 * (-t24 * t84 + t41 * t87) - t107 * t66 - t169 * (t24 * t87 + t41 * t84) + t105 * t41 - t167 * (-t149 * t66 + t42 * t85)) * g(2) + (t94 * mrSges(3,1) + m(4) * t65 - t44 * mrSges(4,1) - m(5) * t143 - t26 * mrSges(5,1) - m(6) * (pkin(4) * t26 + t143) - m(7) * (t26 * t81 + t143) - t107 * t68 - t169 * (t26 * t87 + t43 * t84) + t145 * (-t26 * t84 + t43 * t87) + t105 * t43 - t167 * (-t149 * t68 + t44 * t85)) * g(1) (-t154 * t158 + t163 * t48 - t169 * (-t46 * t147 + t154) + t145 * (t46 * t148 + t48 * t87) + t136 * (-t46 * pkin(3) + pkin(10) * t48) + t159 * t46) * g(3) + (-t156 * t158 - t169 * (-t31 * t147 + t156) + t145 * (t31 * t148 + t33 * t87) + t163 * t33 + t136 * (-t31 * pkin(3) + pkin(10) * t33) + t159 * t31) * g(2) + (-t155 * t158 - t169 * (-t34 * t147 + t155) + t145 * (t34 * t148 + t36 * t87) + t163 * t36 + t136 * (-t34 * pkin(3) + pkin(10) * t36) + t159 * t34) * g(1) (t160 * t27 - t167 * t28) * g(3) + (t160 * t15 - t16 * t167) * g(2) + (t160 * t17 - t167 * t18) * g(1) (t145 * (-t28 * t87 - t46 * t84) + t161 * (-t28 * t84 + t46 * t87)) * g(3) + (t145 * (-t16 * t87 - t31 * t84) + t161 * (-t16 * t84 + t31 * t87)) * g(2) + (t145 * (-t18 * t87 - t34 * t84) + t161 * (-t18 * t84 + t34 * t87)) * g(1) (-g(1) * t17 - g(2) * t15 - g(3) * t27) * m(7)];
taug  = t1(:);
