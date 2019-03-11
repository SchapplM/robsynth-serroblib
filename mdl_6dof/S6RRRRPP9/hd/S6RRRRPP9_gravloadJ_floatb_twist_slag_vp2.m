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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:53
% EndTime: 2019-03-09 21:41:56
% DurationCPUTime: 1.29s
% Computational Cost: add. (837->156), mult. (2085->209), div. (0->0), fcn. (2527->10), ass. (0->86)
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t153 = pkin(4) * t77 + qJ(5) * t74;
t152 = -mrSges(7,1) - mrSges(6,1) - mrSges(5,3);
t143 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t91 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t133 = cos(qJ(1));
t73 = sin(pkin(6));
t111 = t73 * t133;
t132 = sin(qJ(1));
t76 = sin(qJ(2));
t79 = cos(qJ(2));
t115 = cos(pkin(6));
t99 = t115 * t133;
t56 = t132 * t79 + t76 * t99;
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t29 = -t75 * t111 + t56 * t78;
t55 = t132 * t76 - t79 * t99;
t5 = t29 * t74 - t55 * t77;
t6 = t29 * t77 + t55 * t74;
t151 = m(6) + m(7);
t150 = -t78 * mrSges(4,1) + t75 * mrSges(4,2) - mrSges(3,1);
t148 = mrSges(3,2) - mrSges(4,3);
t103 = -t78 * t111 - t56 * t75;
t147 = t153 * t103;
t110 = t73 * t132;
t98 = t115 * t132;
t58 = t133 * t79 - t76 * t98;
t32 = -t78 * t110 + t58 * t75;
t146 = t153 * t32;
t127 = t73 * t76;
t53 = t115 * t78 - t75 * t127;
t145 = t153 * t53;
t144 = t151 * pkin(4) + t91;
t142 = t143 * t74 - t91 * t77 - mrSges(4,1);
t141 = mrSges(4,2) - m(7) * (pkin(5) + pkin(10)) + t152;
t140 = -t152 * t75 - t150;
t137 = pkin(3) * t78;
t135 = pkin(10) * t32;
t134 = t103 * pkin(10);
t130 = t55 * t75;
t57 = t133 * t76 + t79 * t98;
t128 = t57 * t75;
t126 = t73 * t79;
t124 = t74 * t78;
t120 = t77 * t78;
t119 = t78 * t79;
t118 = pkin(2) * t126 + pkin(9) * t127;
t117 = t133 * pkin(1) + pkin(8) * t110;
t114 = t75 * t126;
t113 = t74 * t126;
t109 = -t55 * pkin(2) + pkin(9) * t56;
t108 = -t57 * pkin(2) + pkin(9) * t58;
t22 = t103 * pkin(3);
t107 = pkin(10) * t29 + t22;
t24 = t32 * pkin(3);
t33 = t75 * t110 + t58 * t78;
t106 = pkin(10) * t33 - t24;
t48 = t53 * pkin(3);
t54 = t115 * t75 + t78 * t127;
t105 = pkin(10) * t54 + t48;
t102 = t73 * pkin(3) * t119 + pkin(10) * t114 + t118;
t100 = -t132 * pkin(1) + pkin(8) * t111;
t94 = -pkin(10) * t130 - t55 * t137 + t109;
t93 = -pkin(10) * t128 - t57 * t137 + t108;
t92 = t58 * pkin(2) + pkin(9) * t57 + t117;
t90 = t33 * pkin(3) + t92;
t89 = -t151 * qJ(5) + t143;
t86 = -t56 * pkin(2) - t55 * pkin(9) + t100;
t13 = -t55 * t124 - t56 * t77;
t14 = -t55 * t120 + t56 * t74;
t84 = t14 * pkin(4) + qJ(5) * t13 + t94;
t15 = -t57 * t124 - t58 * t77;
t16 = -t57 * t120 + t58 * t74;
t83 = t16 * pkin(4) + qJ(5) * t15 + t93;
t82 = -pkin(3) * t29 + t86;
t10 = t33 * t77 + t57 * t74;
t9 = t33 * t74 - t57 * t77;
t81 = t10 * pkin(4) + qJ(5) * t9 + t90;
t80 = -pkin(4) * t6 - qJ(5) * t5 + t82;
t36 = (t77 * t119 + t74 * t76) * t73;
t35 = t78 * t113 - t77 * t127;
t27 = t54 * t77 - t113;
t26 = t77 * t126 + t54 * t74;
t1 = [(-t133 * mrSges(2,1) + t132 * mrSges(2,2) - m(3) * t117 - t58 * mrSges(3,1) - mrSges(3,3) * t110 - m(4) * t92 - t33 * mrSges(4,1) - m(5) * (t90 + t135) - m(6) * (t81 + t135) - m(7) * t81 + t148 * t57 + t143 * t9 - t91 * t10 + t141 * t32) * g(2) + (t132 * mrSges(2,1) + t133 * mrSges(2,2) - m(3) * t100 + t56 * mrSges(3,1) - mrSges(3,3) * t111 - m(4) * t86 + t29 * mrSges(4,1) - m(5) * (t82 + t134) - m(6) * (t80 + t134) - m(7) * t80 - t148 * t55 + t91 * t6 - t143 * t5 + t141 * t103) * g(1) (-m(4) * t118 - m(5) * t102 - t151 * (t36 * pkin(4) + t35 * qJ(5) + t102) + (t148 * t76 + t150 * t79) * t73 - t91 * t36 + t143 * t35 + (-m(7) * pkin(5) + t152) * t114) * g(3) + (-m(4) * t109 - m(5) * t94 - m(6) * t84 - m(7) * (-pkin(5) * t130 + t84) + t148 * t56 - t91 * t14 + t143 * t13 + t140 * t55) * g(2) + (-m(4) * t108 - m(5) * t93 - m(6) * t83 - m(7) * (-pkin(5) * t128 + t83) + t148 * t58 - t91 * t16 + t143 * t15 + t140 * t57) * g(1) (-m(5) * t105 - m(6) * (t105 + t145) - m(7) * (t48 + t145) + t141 * t54 + t142 * t53) * g(3) + (-m(5) * t107 - m(6) * (t107 + t147) - m(7) * (t22 + t147) + t141 * t29 + t142 * t103) * g(2) + (-m(5) * t106 - m(6) * (t106 - t146) - m(7) * (-t24 - t146) + t141 * t33 - t142 * t32) * g(1) (t144 * t26 + t89 * t27) * g(3) + (t144 * t5 + t89 * t6) * g(2) + (t89 * t10 + t144 * t9) * g(1), t151 * (-g(1) * t9 - g(2) * t5 - g(3) * t26) (-g(1) * t10 - g(2) * t6 - g(3) * t27) * m(7)];
taug  = t1(:);
