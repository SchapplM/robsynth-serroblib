% Calculate Gravitation load on the joints for
% S6PRRRRP1
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
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:27:32
% EndTime: 2018-11-23 15:27:33
% DurationCPUTime: 1.07s
% Computational Cost: add. (1351->128), mult. (1410->170), div. (0->0), fcn. (1378->16), ass. (0->69)
t167 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t158 = -mrSges(6,2) - mrSges(7,2);
t168 = mrSges(6,1) + mrSges(7,1);
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t169 = t158 * t78 + t81 * t168 + mrSges(5,1);
t113 = sin(pkin(6));
t76 = sin(pkin(11));
t105 = t76 * t113;
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t114 = cos(pkin(11));
t112 = pkin(6) - qJ(2);
t102 = sin(t112);
t111 = pkin(6) + qJ(2);
t96 = sin(t111) / 0.2e1;
t61 = t96 - t102 / 0.2e1;
t83 = cos(qJ(2));
t91 = t114 * t83 - t76 * t61;
t165 = t82 * t105 - t91 * t79;
t115 = cos(pkin(6));
t103 = cos(t112);
t97 = cos(t111) / 0.2e1;
t62 = t97 - t103 / 0.2e1;
t164 = t115 * t82 + t62 * t79;
t75 = qJ(3) + qJ(4);
t73 = sin(t75);
t74 = cos(t75);
t92 = t114 * t61 + t76 * t83;
t93 = t114 * t113;
t29 = t73 * t92 + t74 * t93;
t30 = -t73 * t93 + t74 * t92;
t163 = t167 * t30 + t169 * t29;
t31 = -t105 * t74 + t73 * t91;
t32 = t105 * t73 + t74 * t91;
t162 = t167 * t32 + t169 * t31;
t46 = -t115 * t74 - t62 * t73;
t47 = t115 * t73 - t62 * t74;
t161 = t167 * t47 + t169 * t46;
t71 = pkin(5) * t81 + pkin(4);
t77 = -qJ(6) - pkin(10);
t160 = -m(4) * pkin(2) - t82 * mrSges(4,1) + t79 * mrSges(4,2) - mrSges(3,1) + (-m(6) * pkin(4) - m(7) * t71 - mrSges(5,1)) * t74 + (-m(6) * pkin(10) + m(7) * t77 + t167) * t73;
t146 = m(7) * pkin(5);
t154 = -m(5) - m(6) - m(7);
t152 = -t146 - t168;
t151 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t145 = -t29 * t71 - t30 * t77;
t138 = t92 * t78;
t136 = t91 * t78;
t135 = t62 * t78;
t131 = t74 * t78;
t130 = t74 * t81;
t129 = -t31 * t71 - t32 * t77;
t122 = -t46 * t71 - t47 * t77;
t109 = -t29 * pkin(4) + t30 * pkin(10);
t108 = -t31 * pkin(4) + pkin(10) * t32;
t107 = -t46 * pkin(4) + pkin(10) * t47;
t101 = t165 * pkin(3);
t100 = t164 * pkin(3);
t88 = -t79 * t92 - t82 * t93;
t86 = t103 / 0.2e1 + t97;
t85 = t88 * pkin(3);
t84 = -pkin(9) - pkin(8);
t80 = sin(qJ(2));
t72 = pkin(3) * t82 + pkin(2);
t60 = t96 + t102 / 0.2e1;
t53 = t114 * t80 + t76 * t86;
t50 = -t114 * t86 + t76 * t80;
t1 = [(-m(2) - m(3) - m(4) + t154) * g(3) (t135 * t146 - t168 * (t130 * t60 - t135) + t158 * (-t131 * t60 - t62 * t81) + t154 * (t60 * t72 + t62 * t84) - t151 * t62 + t160 * t60) * g(3) + (-t138 * t146 - t168 * (-t130 * t50 + t138) + t158 * (t131 * t50 + t81 * t92) + t154 * (-t50 * t72 - t84 * t92) + t151 * t92 - t160 * t50) * g(2) + (-t136 * t146 - t168 * (-t130 * t53 + t136) + t158 * (t131 * t53 + t81 * t91) + t154 * (-t53 * t72 - t84 * t91) + t151 * t91 - t160 * t53) * g(1) (-t164 * mrSges(4,1) - (-t115 * t79 + t62 * t82) * mrSges(4,2) - m(5) * t100 - m(6) * (t100 + t107) - m(7) * (t100 + t122) + t161) * g(3) + (-t88 * mrSges(4,1) - (t79 * t93 - t82 * t92) * mrSges(4,2) - m(5) * t85 - m(6) * (t109 + t85) - m(7) * (t85 + t145) + t163) * g(2) + (-t165 * mrSges(4,1) - (-t105 * t79 - t82 * t91) * mrSges(4,2) - m(5) * t101 - m(6) * (t101 + t108) - m(7) * (t101 + t129) + t162) * g(1) (-m(6) * t107 - m(7) * t122 + t161) * g(3) + (-m(6) * t109 - m(7) * t145 + t163) * g(2) + (-m(6) * t108 - m(7) * t129 + t162) * g(1) (t158 * (-t47 * t81 + t60 * t78) + t152 * (-t47 * t78 - t60 * t81)) * g(3) + (t158 * (-t30 * t81 - t50 * t78) + t152 * (-t30 * t78 + t50 * t81)) * g(2) + (t158 * (-t32 * t81 - t53 * t78) + t152 * (-t32 * t78 + t53 * t81)) * g(1) (-g(1) * t31 - g(2) * t29 - g(3) * t46) * m(7)];
taug  = t1(:);
