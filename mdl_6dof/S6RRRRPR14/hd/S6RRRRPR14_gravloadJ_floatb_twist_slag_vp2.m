% Calculate Gravitation load on the joints for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:45
% EndTime: 2019-03-10 00:08:51
% DurationCPUTime: 1.99s
% Computational Cost: add. (1553->160), mult. (4056->233), div. (0->0), fcn. (5136->16), ass. (0->85)
t72 = sin(pkin(13));
t73 = cos(pkin(13));
t152 = -m(7) * (pkin(5) * t73 + pkin(4)) - m(6) * pkin(4) - t73 * mrSges(6,1) + t72 * mrSges(6,2) - mrSges(5,1);
t71 = pkin(13) + qJ(6);
t68 = sin(t71);
t69 = cos(t71);
t150 = -t69 * mrSges(7,1) + t68 * mrSges(7,2) + t152;
t142 = cos(pkin(6));
t148 = cos(qJ(2));
t127 = t142 * t148;
t145 = sin(qJ(2));
t146 = sin(qJ(1));
t149 = cos(qJ(1));
t100 = -t127 * t149 + t146 * t145;
t139 = sin(pkin(7));
t140 = sin(pkin(6));
t116 = t140 * t139;
t141 = cos(pkin(7));
t180 = t100 * t141 + t149 * t116;
t147 = cos(qJ(3));
t126 = t142 * t145;
t52 = t126 * t149 + t146 * t148;
t76 = sin(qJ(3));
t22 = -t147 * t52 + t180 * t76;
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t117 = t141 * t140;
t84 = -t100 * t139 + t149 * t117;
t179 = t22 * t75 - t77 * t84;
t178 = t22 * t77 + t75 * t84;
t172 = -m(5) - m(6);
t154 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(12) - qJ(5)) - mrSges(7,3);
t92 = t146 * t127 + t145 * t149;
t78 = t146 * t117 + t92 * t139;
t171 = m(6) + m(7);
t174 = mrSges(4,1) - t150 * t77 - t154 * t75 + pkin(3) * (m(5) + t171);
t153 = -t72 * mrSges(6,1) - t73 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3) - m(7) * (pkin(5) * t72 + pkin(11));
t151 = -t68 * mrSges(7,1) - t69 * mrSges(7,2) + t153;
t19 = t180 * t147 + t52 * t76;
t163 = -t139 * mrSges(4,3) + mrSges(3,2);
t162 = t146 * t116 - t92 * t141;
t160 = t148 * t117 + t142 * t139;
t155 = t172 * pkin(11) + t151;
t103 = t145 * t116;
t123 = t148 * t140;
t144 = pkin(2) * t123 + pkin(10) * t103;
t122 = t140 * t146;
t143 = t149 * pkin(1) + pkin(9) * t122;
t107 = t145 * t117;
t46 = -t107 * t76 + t123 * t147;
t138 = t46 * pkin(3) + t144;
t135 = pkin(10) * t139;
t133 = t75 * t139;
t132 = t76 * t141;
t131 = t77 * t139;
t124 = t149 * t140;
t129 = -pkin(1) * t146 + pkin(9) * t124;
t125 = t141 * t147;
t121 = t140 * t145;
t115 = -t100 * pkin(2) + t135 * t52;
t53 = -t146 * t126 + t148 * t149;
t114 = -t92 * pkin(2) + t135 * t53;
t28 = -t100 * t147 - t132 * t52;
t112 = t28 * pkin(3) + t115;
t30 = -t132 * t53 - t147 * t92;
t111 = t30 * pkin(3) + t114;
t91 = -t116 * t148 + t141 * t142;
t86 = -t52 * pkin(2) + t84 * pkin(10) + t129;
t83 = t22 * pkin(3) + t86;
t82 = t53 * pkin(2) + t78 * pkin(10) + t143;
t24 = t53 * t147 + t162 * t76;
t81 = t24 * pkin(3) + t82;
t45 = t107 * t147 + t123 * t76;
t38 = t147 * t121 + t160 * t76;
t37 = t121 * t76 - t160 * t147;
t29 = t125 * t53 - t76 * t92;
t27 = -t100 * t76 + t125 * t52;
t23 = -t162 * t147 + t53 * t76;
t18 = t38 * t77 + t75 * t91;
t17 = t38 * t75 - t77 * t91;
t8 = t24 * t77 + t75 * t78;
t7 = t24 * t75 - t77 * t78;
t2 = t23 * t68 + t69 * t8;
t1 = t23 * t69 - t68 * t8;
t3 = [(-m(3) * t143 - m(4) * t82 - m(7) * t81 - t149 * mrSges(2,1) - t53 * mrSges(3,1) - t24 * mrSges(4,1) - t2 * mrSges(7,1) + t146 * mrSges(2,2) + t92 * mrSges(3,2) - t1 * mrSges(7,2) - mrSges(3,3) * t122 - t78 * mrSges(4,3) + t172 * (t23 * pkin(11) + t81) + t152 * t8 + t153 * t23 + t154 * t7) * g(2) + (-m(3) * t129 - m(4) * t86 - m(7) * t83 + t146 * mrSges(2,1) + t52 * mrSges(3,1) - t22 * mrSges(4,1) + t149 * mrSges(2,2) - t100 * mrSges(3,2) - mrSges(3,3) * t124 - t84 * mrSges(4,3) + t172 * (-pkin(11) * t19 + t83) + t150 * t178 - t151 * t19 + t154 * t179) * g(1) (-m(4) * t144 - m(7) * t138 - mrSges(3,1) * t123 - t46 * mrSges(4,1) + mrSges(3,2) * t121 - mrSges(4,3) * t103 + t172 * (pkin(11) * t45 + t138) + t150 * (t103 * t75 + t46 * t77) + t151 * t45 + t154 * (-t103 * t77 + t46 * t75)) * g(3) + (-m(4) * t115 - m(7) * t112 + mrSges(3,1) * t100 - t28 * mrSges(4,1) + t172 * (t27 * pkin(11) + t112) + t163 * t52 + t150 * (t133 * t52 + t28 * t77) + t151 * t27 + t154 * (-t131 * t52 + t28 * t75)) * g(2) + (-m(4) * t114 - m(7) * t111 + mrSges(3,1) * t92 - t30 * mrSges(4,1) + t172 * (t29 * pkin(11) + t111) + t163 * t53 + t150 * (t133 * t53 + t30 * t77) + t151 * t29 + t154 * (-t131 * t53 + t30 * t75)) * g(1) (t155 * t38 + t174 * t37) * g(3) + (-t155 * t22 + t174 * t19) * g(2) + (t155 * t24 + t174 * t23) * g(1) (-t150 * t17 + t154 * t18) * g(3) + (t150 * t179 - t154 * t178) * g(2) + (-t150 * t7 + t154 * t8) * g(1), t171 * (-g(1) * t7 + g(2) * t179 - g(3) * t17) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t178 * t68 + t19 * t69) * mrSges(7,1) + (t178 * t69 - t19 * t68) * mrSges(7,2)) - g(3) * ((-t18 * t68 + t37 * t69) * mrSges(7,1) + (-t18 * t69 - t37 * t68) * mrSges(7,2))];
taug  = t3(:);
