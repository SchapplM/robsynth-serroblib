% Calculate Gravitation load on the joints for
% S6RPRRRP12
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:04
% EndTime: 2019-03-09 06:44:08
% DurationCPUTime: 1.43s
% Computational Cost: add. (1372->119), mult. (3742->173), div. (0->0), fcn. (4786->14), ass. (0->66)
t131 = sin(pkin(6));
t133 = cos(pkin(7));
t106 = t133 * t131;
t130 = sin(pkin(7));
t142 = cos(qJ(1));
t132 = cos(pkin(12));
t134 = cos(pkin(6));
t109 = t134 * t132;
t129 = sin(pkin(12));
t140 = sin(qJ(1));
t96 = -t142 * t109 + t140 * t129;
t156 = -t142 * t106 + t96 * t130;
t141 = cos(qJ(3));
t105 = t131 * t130;
t161 = t142 * t105 + t96 * t133;
t107 = t134 * t129;
t63 = t142 * t107 + t140 * t132;
t77 = sin(qJ(3));
t44 = -t63 * t141 + t161 * t77;
t76 = sin(qJ(4));
t79 = cos(qJ(4));
t20 = -t156 * t76 + t44 * t79;
t41 = t141 * t161 + t63 * t77;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t165 = t20 * t75 + t41 * t78;
t164 = t20 * t78 - t41 * t75;
t145 = -m(6) - m(7);
t149 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t102 = t145 * pkin(11) + t149;
t17 = t156 * t79 + t44 * t76;
t121 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t122 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t148 = t121 * t75 - t122 * t78 - mrSges(5,1);
t158 = -m(5) + t145;
t90 = t140 * t109 + t142 * t129;
t157 = t140 * t106 + t90 * t130;
t153 = mrSges(4,2) - mrSges(5,3);
t155 = -t121 * t78 - t122 * t75 + t153;
t154 = mrSges(4,1) - t102 * t76 + (-t145 * pkin(4) - t148) * t79;
t152 = -t140 * t105 + t90 * t133;
t150 = t132 * t106 + t134 * t130;
t117 = t131 * t140;
t135 = t142 * pkin(1) + qJ(2) * t117;
t118 = t142 * t131;
t120 = -t140 * pkin(1) + qJ(2) * t118;
t104 = t131 * t129;
t103 = t158 * pkin(10) + t153;
t87 = -t63 * pkin(2) - t156 * pkin(9) + t120;
t85 = t44 * pkin(3) + t87;
t64 = -t140 * t107 + t142 * t132;
t84 = t64 * pkin(2) + t157 * pkin(9) + t135;
t46 = t64 * t141 - t152 * t77;
t82 = t46 * pkin(3) + t84;
t62 = -t132 * t105 + t134 * t133;
t55 = t141 * t104 + t150 * t77;
t54 = t77 * t104 - t150 * t141;
t45 = t152 * t141 + t64 * t77;
t40 = t55 * t79 + t62 * t76;
t39 = -t55 * t76 + t62 * t79;
t22 = t157 * t76 + t46 * t79;
t21 = -t157 * t79 + t46 * t76;
t11 = t40 * t75 - t54 * t78;
t6 = t22 * t78 + t45 * t75;
t5 = t22 * t75 - t45 * t78;
t1 = [(-m(3) * t135 - m(4) * t84 - m(5) * t82 - t142 * mrSges(2,1) - t64 * mrSges(3,1) - t46 * mrSges(4,1) - t22 * mrSges(5,1) + t140 * mrSges(2,2) + t90 * mrSges(3,2) - mrSges(3,3) * t117 - t157 * mrSges(4,3) + t102 * t21 + t103 * t45 + t121 * t5 - t122 * t6 + t145 * (t22 * pkin(4) + t82)) * g(2) + (t140 * mrSges(2,1) + t142 * mrSges(2,2) - m(3) * t120 + t63 * mrSges(3,1) - t96 * mrSges(3,2) - mrSges(3,3) * t118 - m(4) * t87 - t44 * mrSges(4,1) + t156 * mrSges(4,3) - m(5) * t85 - t20 * mrSges(5,1) - t103 * t41 - t122 * t164 + t121 * t165 + t102 * t17 + t145 * (t20 * pkin(4) + t85)) * g(1) (-g(1) * t117 + g(2) * t118 - g(3) * t134) * (m(3) + m(4) - t158) (t158 * (-t54 * pkin(3) + pkin(10) * t55) + t155 * t55 + t154 * t54) * g(3) + (t158 * (-t41 * pkin(3) - pkin(10) * t44) - t155 * t44 + t154 * t41) * g(2) + (t158 * (-t45 * pkin(3) + pkin(10) * t46) + t155 * t46 + t154 * t45) * g(1) (t145 * (t39 * pkin(4) + pkin(11) * t40) + t149 * t40 + t148 * t39) * g(3) + (t145 * (t17 * pkin(4) - pkin(11) * t20) - t149 * t20 + t148 * t17) * g(2) + (t145 * (-t21 * pkin(4) + pkin(11) * t22) + t149 * t22 - t148 * t21) * g(1) (t121 * (t40 * t78 + t54 * t75) + t122 * t11) * g(3) + (-t121 * t164 - t122 * t165) * g(2) + (t121 * t6 + t122 * t5) * g(1) (-g(1) * t5 + g(2) * t165 - g(3) * t11) * m(7)];
taug  = t1(:);
