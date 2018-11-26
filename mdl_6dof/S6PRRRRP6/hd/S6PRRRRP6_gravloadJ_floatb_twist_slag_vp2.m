% Calculate Gravitation load on the joints for
% S6PRRRRP6
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
% Datum: 2018-11-23 15:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:31:10
% EndTime: 2018-11-23 15:31:11
% DurationCPUTime: 1.44s
% Computational Cost: add. (4573->176), mult. (4763->255), div. (0->0), fcn. (4689->22), ass. (0->109)
t139 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t138 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t109 = cos(qJ(3));
t158 = pkin(6) + qJ(2);
t137 = cos(t158) / 0.2e1;
t159 = pkin(6) - qJ(2);
t144 = cos(t159);
t121 = t144 / 0.2e1 + t137;
t165 = sin(pkin(12));
t167 = cos(pkin(12));
t179 = sin(qJ(2));
t111 = -t167 * t121 + t165 * t179;
t156 = pkin(7) + qJ(3);
t133 = sin(t156) / 0.2e1;
t157 = pkin(7) - qJ(3);
t140 = sin(t157);
t183 = t133 - t140 / 0.2e1;
t110 = cos(qJ(2));
t134 = sin(t158) / 0.2e1;
t141 = sin(t159);
t96 = t134 - t141 / 0.2e1;
t87 = t165 * t110 + t167 * t96;
t189 = t87 * t109 - t111 * t183;
t112 = t165 * t121 + t167 * t179;
t89 = t167 * t110 - t165 * t96;
t188 = t89 * t109 - t112 * t183;
t118 = t134 + t141 / 0.2e1;
t97 = t137 - t144 / 0.2e1;
t187 = -t97 * t109 + t118 * t183;
t180 = -m(6) - m(7);
t186 = mrSges(4,2) - mrSges(5,3);
t185 = mrSges(6,3) + mrSges(7,2);
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t184 = t108 * mrSges(5,1) - mrSges(5,2) * t105 + mrSges(4,1);
t182 = mrSges(5,2) - t185;
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t181 = t138 * t104 - t139 * t107 - mrSges(5,1);
t178 = pkin(4) * t108;
t58 = -t111 * t109 - t183 * t87;
t84 = t111 * pkin(2);
t177 = t58 * pkin(3) - t84;
t60 = -t112 * t109 - t183 * t89;
t85 = t112 * pkin(2);
t176 = t60 * pkin(3) - t85;
t72 = t118 * t109 + t183 * t97;
t95 = t118 * pkin(2);
t175 = t72 * pkin(3) + t95;
t106 = sin(qJ(3));
t116 = t133 + t140 / 0.2e1;
t166 = sin(pkin(6));
t113 = t116 * t166;
t142 = cos(t156);
t135 = t142 / 0.2e1;
t143 = cos(t157);
t136 = t143 / 0.2e1;
t120 = t136 + t135;
t44 = t87 * t106 + t111 * t120 + t167 * t113;
t174 = t105 * t44;
t47 = t106 * t89 + t112 * t120 - t165 * t113;
t173 = t105 * t47;
t168 = cos(pkin(6));
t66 = -t106 * t97 - t168 * t116 - t118 * t120;
t172 = t105 * t66;
t102 = sin(pkin(7));
t164 = t102 * t105;
t163 = t102 * t108;
t162 = t104 * t108;
t161 = t107 * t108;
t160 = -m(5) + t180;
t152 = m(4) - t160;
t119 = t135 - t143 / 0.2e1;
t114 = t119 * t166;
t46 = t167 * t114 + t189;
t151 = -t44 * pkin(3) + pkin(10) * t46;
t49 = -t165 * t114 + t188;
t150 = -t47 * pkin(3) + pkin(10) * t49;
t68 = -t168 * t119 + t187;
t149 = -t66 * pkin(3) + pkin(10) * t68;
t103 = cos(pkin(7));
t145 = t103 * t166;
t126 = t136 - t142 / 0.2e1;
t125 = t160 * pkin(10) + t186;
t124 = t180 * pkin(11) + t182;
t123 = t126 * t166;
t122 = -mrSges(3,2) + (t152 * pkin(9) + mrSges(4,3)) * t102;
t86 = -t118 * t102 + t168 * t103;
t74 = t112 * t102 + t165 * t145;
t73 = t111 * t102 - t167 * t145;
t71 = t118 * t106 - t97 * t120;
t67 = t168 * t126 + t187;
t59 = -t112 * t106 + t120 * t89;
t57 = -t111 * t106 + t120 * t87;
t56 = t108 * t72 - t97 * t164;
t48 = t165 * t123 + t188;
t45 = -t167 * t123 + t189;
t33 = t105 * t86 + t108 * t67;
t32 = -t105 * t67 + t108 * t86;
t30 = t108 * t60 + t164 * t89;
t28 = t108 * t58 + t164 * t87;
t20 = t105 * t74 + t108 * t48;
t19 = -t105 * t48 + t108 * t74;
t18 = t105 * t73 + t108 * t45;
t17 = -t105 * t45 + t108 * t73;
t13 = t104 * t33 - t66 * t107;
t3 = t104 * t20 - t47 * t107;
t1 = t104 * t18 - t44 * t107;
t2 = [(-m(2) - m(3) - t152) * g(3) (-t118 * mrSges(3,1) - m(4) * t95 - t72 * mrSges(4,1) - m(5) * t175 - t56 * mrSges(5,1) + t122 * t97 + t125 * t71 - t139 * (t104 * t71 + t107 * t56) + t138 * (t104 * t56 - t71 * t107) + t124 * (t105 * t72 + t97 * t163) + t180 * (t56 * pkin(4) + t175)) * g(3) + (t111 * mrSges(3,1) + m(4) * t84 - t58 * mrSges(4,1) - m(5) * t177 - t28 * mrSges(5,1) + t138 * (t104 * t28 - t57 * t107) - t122 * t87 + t125 * t57 - t139 * (t104 * t57 + t107 * t28) + t124 * (t105 * t58 - t163 * t87) + t180 * (t28 * pkin(4) + t177)) * g(2) + (t112 * mrSges(3,1) + m(4) * t85 - t60 * mrSges(4,1) - m(5) * t176 - t30 * mrSges(5,1) - t122 * t89 + t125 * t59 - t139 * (t104 * t59 + t107 * t30) + t138 * (t104 * t30 - t59 * t107) + t124 * (t105 * t60 - t163 * t89) + t180 * (t30 * pkin(4) + t176)) * g(1) (-m(5) * t149 + t186 * t68 + t184 * t66 + t180 * (-pkin(11) * t172 - t66 * t178 + t149) - t139 * (t104 * t68 - t66 * t161) + t138 * (-t68 * t107 - t66 * t162) + t185 * t172) * g(3) + (-m(5) * t151 + t180 * (-pkin(11) * t174 - t44 * t178 + t151) - t139 * (t104 * t46 - t44 * t161) + t138 * (-t46 * t107 - t44 * t162) + t186 * t46 + t184 * t44 + t185 * t174) * g(2) + (-m(5) * t150 + t180 * (-pkin(11) * t173 - t47 * t178 + t150) - t139 * (t104 * t49 - t47 * t161) + t138 * (-t49 * t107 - t47 * t162) + t186 * t49 + t184 * t47 + t185 * t173) * g(1) (t180 * (t32 * pkin(4) + pkin(11) * t33) + t182 * t33 + t181 * t32) * g(3) + (t180 * (t17 * pkin(4) + pkin(11) * t18) + t182 * t18 + t181 * t17) * g(2) + (t180 * (t19 * pkin(4) + pkin(11) * t20) + t182 * t20 + t181 * t19) * g(1) (t138 * (t104 * t66 + t107 * t33) + t139 * t13) * g(3) + (t138 * (t104 * t44 + t107 * t18) + t139 * t1) * g(2) + (t138 * (t104 * t47 + t107 * t20) + t139 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t13) * m(7)];
taug  = t2(:);
