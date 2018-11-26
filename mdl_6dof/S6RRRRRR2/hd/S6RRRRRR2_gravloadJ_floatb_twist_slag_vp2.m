% Calculate Gravitation load on the joints for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:38:08
% EndTime: 2018-11-23 18:38:09
% DurationCPUTime: 0.97s
% Computational Cost: add. (699->135), mult. (598->147), div. (0->0), fcn. (521->12), ass. (0->87)
t49 = qJ(5) + qJ(6);
t42 = sin(t49);
t51 = sin(qJ(5));
t155 = -t51 * mrSges(6,2) - t42 * mrSges(7,2);
t44 = cos(t49);
t54 = cos(qJ(5));
t154 = mrSges(6,1) * t54 + t44 * mrSges(7,1);
t153 = -mrSges(6,3) - mrSges(7,3);
t50 = qJ(2) + qJ(3);
t46 = qJ(4) + t50;
t38 = sin(t46);
t152 = t154 * t38;
t151 = t155 * t38;
t39 = cos(t46);
t40 = pkin(5) * t54 + pkin(4);
t57 = -pkin(11) - pkin(10);
t73 = -t38 * t40 - t39 * t57;
t150 = -m(7) * t73 + t152;
t128 = m(7) * pkin(5);
t43 = sin(t50);
t45 = cos(t50);
t76 = mrSges(5,1) * t38 + mrSges(5,2) * t39;
t149 = mrSges(4,1) * t43 + mrSges(4,2) * t45 + t76;
t56 = cos(qJ(1));
t114 = t39 * t56;
t148 = t153 * t114 + t151 * t56;
t53 = sin(qJ(1));
t115 = t39 * t53;
t147 = t153 * t115 + t151 * t53;
t146 = g(1) * t56 + g(2) * t53;
t145 = -t39 * mrSges(5,1) + (mrSges(5,2) + t153) * t38;
t136 = -t38 * t57 + t39 * t40;
t142 = t39 * pkin(4) + t38 * pkin(10);
t144 = -m(6) * t142 - m(7) * t136;
t143 = t51 * t128;
t141 = t150 * t56 + t148;
t140 = t150 * t53 + t147;
t123 = pkin(4) * t38;
t124 = pkin(3) * t43;
t139 = -m(7) * (t73 - t124) - m(6) * (-t123 - t124) + t152;
t137 = mrSges(6,1) + t128;
t88 = t45 * mrSges(4,1) - t43 * mrSges(4,2);
t134 = m(5) + m(6) + m(7);
t133 = t145 + (-t154 - t155) * t39;
t58 = -pkin(8) - pkin(7);
t131 = -m(3) * pkin(7) + m(4) * t58 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t130 = -t88 + t133;
t55 = cos(qJ(2));
t47 = t55 * pkin(2);
t52 = sin(qJ(2));
t80 = t55 * mrSges(3,1) - t52 * mrSges(3,2);
t129 = mrSges(2,1) + m(4) * (t47 + pkin(1)) + t88 + m(3) * pkin(1) + t80 - t145;
t107 = t44 * t56;
t112 = t42 * t53;
t5 = t112 * t39 + t107;
t108 = t44 * t53;
t111 = t42 * t56;
t6 = -t108 * t39 + t111;
t127 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t111 * t39 + t108;
t8 = t107 * t39 + t112;
t126 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t125 = pkin(2) * t52;
t37 = pkin(3) * t45;
t119 = g(3) * t38;
t105 = t51 * t53;
t104 = t51 * t56;
t103 = t53 * t54;
t102 = t54 * t56;
t96 = t37 + t47;
t91 = t37 + t142;
t29 = pkin(10) * t115;
t84 = -t123 * t53 + t29;
t30 = pkin(10) * t114;
t83 = -t123 * t56 + t30;
t82 = t37 + t136;
t75 = -mrSges(7,1) * t42 - mrSges(7,2) * t44;
t11 = -t104 * t39 + t103;
t9 = t105 * t39 + t102;
t48 = -pkin(9) + t58;
t23 = -t124 - t125;
t22 = pkin(1) + t96;
t17 = t56 * t23;
t16 = t53 * t23;
t12 = t102 * t39 + t105;
t10 = -t103 * t39 + t104;
t1 = [(-t105 * t128 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t134 * (t56 * t22 - t48 * t53) + t131 * t53 + (-t129 + t144) * t56) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t134 * t48 + t131 - t143) * t56 + (m(5) * t22 - m(6) * (-t22 - t142) - m(7) * (-t22 - t136) + t129) * t53) * g(1) (-m(6) * (t16 + t84) - m(7) * t16 + t140) * g(2) + (-m(6) * (t17 + t83) - m(7) * t17 + t141) * g(1) + (-t80 - m(4) * t47 - m(5) * t96 - m(6) * (t47 + t91) - m(7) * (t47 + t82) + t130) * g(3) + t146 * (m(4) * t125 - m(5) * t23 + mrSges(3,1) * t52 + mrSges(3,2) * t55 + t149) (-m(6) * t29 + t139 * t53 + t147) * g(2) + (-m(6) * t30 + t139 * t56 + t148) * g(1) + (-m(5) * t37 - m(6) * t91 - m(7) * t82 + t130) * g(3) + t146 * (m(5) * t124 + t149) t146 * t76 + (-m(6) * t84 + t140) * g(2) + (-m(6) * t83 + t141) * g(1) + (t133 + t144) * g(3) (mrSges(6,1) * t51 + mrSges(6,2) * t54 + t143 - t75) * t119 + (-t10 * mrSges(6,2) + t137 * t9 - t127) * g(2) + (t12 * mrSges(6,2) - t137 * t11 - t126) * g(1), -g(1) * t126 - g(2) * t127 - t119 * t75];
taug  = t1(:);
