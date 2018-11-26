% Calculate Gravitation load on the joints for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:10:41
% EndTime: 2018-11-23 18:10:43
% DurationCPUTime: 1.50s
% Computational Cost: add. (1971->167), mult. (2463->212), div. (0->0), fcn. (2527->14), ass. (0->97)
t167 = -mrSges(5,2) + mrSges(7,2) + mrSges(6,3);
t100 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t164 = pkin(4) * t84 + qJ(5) * t80;
t163 = mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t162 = t100 * t84 + t167 * t80 + mrSges(4,1);
t129 = sin(pkin(6));
t144 = cos(qJ(1));
t110 = t144 * t129;
t128 = pkin(6) - qJ(2);
t113 = sin(t128);
t127 = pkin(6) + qJ(2);
t112 = sin(t127);
t75 = t112 / 0.2e1;
t132 = t75 - t113 / 0.2e1;
t83 = sin(qJ(1));
t86 = cos(qJ(2));
t134 = t83 * t86;
t56 = t132 * t144 + t134;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t30 = -t110 * t81 + t56 * t85;
t82 = sin(qJ(2));
t109 = cos(t127) / 0.2e1;
t114 = cos(t128);
t89 = t114 / 0.2e1 + t109;
t55 = -t144 * t89 + t82 * t83;
t5 = t30 * t80 - t55 * t84;
t6 = t30 * t84 + t55 * t80;
t161 = m(6) + m(7);
t160 = mrSges(3,2) - mrSges(4,3);
t116 = -t110 * t85 - t56 * t81;
t159 = t164 * t116;
t115 = t83 * t129;
t125 = t144 * t86;
t59 = -t132 * t83 + t125;
t33 = -t115 * t85 + t59 * t81;
t158 = t164 * t33;
t69 = t109 - t114 / 0.2e1;
t79 = cos(pkin(6));
t53 = t69 * t81 + t79 * t85;
t157 = t164 * t53;
t156 = pkin(4) * t161 + t100;
t153 = -m(7) * (pkin(5) + pkin(10)) - mrSges(7,1) - t163;
t151 = mrSges(4,1) * t85 + t163 * t81 + mrSges(3,1);
t148 = pkin(3) * t85;
t146 = pkin(10) * t116;
t145 = pkin(10) * t33;
t142 = t55 * t81;
t58 = t144 * t82 + t83 * t89;
t140 = t58 * t81;
t108 = t113 / 0.2e1;
t68 = t75 + t108;
t139 = t68 * t81;
t137 = t80 * t85;
t133 = t84 * t85;
t131 = pkin(1) * t144 + pkin(8) * t115;
t124 = -pkin(1) * t83 + pkin(8) * t110;
t88 = t108 - t112 / 0.2e1;
t57 = -t144 * t88 + t134;
t123 = -pkin(2) * t55 + pkin(9) * t57;
t60 = t83 * t88 + t125;
t122 = -pkin(2) * t58 + pkin(9) * t60;
t121 = pkin(2) * t68 - pkin(9) * t69;
t24 = t116 * pkin(3);
t120 = pkin(10) * t30 + t24;
t26 = t33 * pkin(3);
t34 = t115 * t81 + t59 * t85;
t119 = pkin(10) * t34 - t26;
t48 = t53 * pkin(3);
t54 = -t69 * t85 + t79 * t81;
t118 = pkin(10) * t54 + t48;
t104 = -pkin(10) * t142 - t148 * t55 + t123;
t103 = -pkin(10) * t140 - t148 * t58 + t122;
t102 = pkin(2) * t59 + pkin(9) * t58 + t131;
t101 = pkin(10) * t139 + t148 * t68 + t121;
t99 = pkin(3) * t34 + t102;
t98 = -pkin(2) * t56 - pkin(9) * t55 + t124;
t97 = -qJ(5) * t161 - t167;
t95 = -pkin(3) * t30 + t98;
t13 = -t137 * t55 - t57 * t84;
t14 = -t133 * t55 + t57 * t80;
t93 = pkin(4) * t14 + qJ(5) * t13 + t104;
t15 = -t137 * t58 - t60 * t84;
t16 = -t133 * t58 + t60 * t80;
t92 = pkin(4) * t16 + qJ(5) * t15 + t103;
t35 = t137 * t68 + t69 * t84;
t36 = t133 * t68 - t69 * t80;
t91 = pkin(4) * t36 + qJ(5) * t35 + t101;
t10 = t34 * t84 + t58 * t80;
t9 = t34 * t80 - t58 * t84;
t90 = pkin(4) * t10 + qJ(5) * t9 + t99;
t87 = -pkin(4) * t6 - qJ(5) * t5 + t95;
t19 = t54 * t84 - t68 * t80;
t18 = t54 * t80 + t68 * t84;
t1 = [(-t144 * mrSges(2,1) + t83 * mrSges(2,2) - m(3) * t131 - t59 * mrSges(3,1) - mrSges(3,3) * t115 - m(4) * t102 - t34 * mrSges(4,1) - m(5) * (t99 + t145) - m(6) * (t90 + t145) - m(7) * t90 + t160 * t58 - t167 * t9 - t100 * t10 + t153 * t33) * g(2) + (t83 * mrSges(2,1) + t144 * mrSges(2,2) - m(3) * t124 + t56 * mrSges(3,1) - mrSges(3,3) * t110 - m(4) * t98 + t30 * mrSges(4,1) - m(5) * (t95 + t146) - m(6) * (t87 + t146) - m(7) * t87 - t160 * t55 + t100 * t6 + t167 * t5 + t153 * t116) * g(1) (-m(4) * t121 - m(5) * t101 - m(6) * t91 - m(7) * (pkin(5) * t139 + t91) - mrSges(7,1) * t139 - t160 * t69 - t100 * t36 - t167 * t35 - t151 * t68) * g(3) + (-m(4) * t123 - m(5) * t104 - m(6) * t93 - m(7) * (-pkin(5) * t142 + t93) + mrSges(7,1) * t142 + t160 * t57 - t100 * t14 - t167 * t13 + t151 * t55) * g(2) + (-m(4) * t122 - m(5) * t103 - m(6) * t92 - m(7) * (-pkin(5) * t140 + t92) + mrSges(7,1) * t140 + t160 * t60 - t100 * t16 - t167 * t15 + t151 * t58) * g(1) (-m(5) * t118 - m(6) * (t118 + t157) - m(7) * (t48 + t157) + t153 * t54 - t162 * t53) * g(3) + (-m(5) * t120 - m(6) * (t120 + t159) - m(7) * (t24 + t159) + t153 * t30 - t162 * t116) * g(2) + (-m(5) * t119 - m(6) * (t119 - t158) - m(7) * (-t26 - t158) + t153 * t34 + t162 * t33) * g(1) (t156 * t18 + t97 * t19) * g(3) + (t156 * t5 + t97 * t6) * g(2) + (t10 * t97 + t156 * t9) * g(1), t161 * (-g(1) * t9 - g(2) * t5 - g(3) * t18) (-g(1) * t10 - g(2) * t6 - g(3) * t19) * m(7)];
taug  = t1(:);
