% Calculate Gravitation load on the joints for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:19
% EndTime: 2019-03-08 19:07:21
% DurationCPUTime: 1.13s
% Computational Cost: add. (1746->113), mult. (5008->190), div. (0->0), fcn. (6532->18), ass. (0->76)
t128 = m(6) + m(7);
t98 = -m(5) - t128;
t134 = m(3) + m(4) - t98;
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t133 = -t51 * mrSges(7,1) - t54 * mrSges(7,2) - pkin(11) * t128 + mrSges(5,2) - mrSges(6,3);
t132 = m(7) * pkin(5) + t54 * mrSges(7,1) - t51 * mrSges(7,2) + mrSges(6,1);
t95 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t130 = pkin(4) * t128 + t132 * t55 - t95 * t52 + mrSges(5,1);
t105 = cos(pkin(8));
t50 = sin(pkin(8));
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t104 = cos(pkin(13));
t106 = cos(pkin(7));
t107 = cos(pkin(6));
t100 = sin(pkin(13));
t99 = sin(pkin(14));
t83 = t100 * t99;
t103 = cos(pkin(14));
t90 = t104 * t103;
t71 = t107 * t90 - t83;
t101 = sin(pkin(7));
t102 = sin(pkin(6));
t87 = t102 * t101;
t123 = -t104 * t87 + t71 * t106;
t84 = t100 * t103;
t88 = t104 * t99;
t72 = t107 * t88 + t84;
t58 = t72 * t113 - t123 * t115;
t89 = t104 * t102;
t64 = -t71 * t101 - t106 * t89;
t126 = t58 * t105 - t64 * t50;
t73 = -t107 * t84 - t88;
t86 = t102 * t100;
t122 = t101 * t86 + t73 * t106;
t74 = -t107 * t83 + t90;
t59 = t74 * t113 - t122 * t115;
t65 = -t73 * t101 + t106 * t86;
t125 = t59 * t105 - t65 * t50;
t121 = t106 * t103 * t102 + t107 * t101;
t85 = t102 * t99;
t63 = t113 * t85 - t121 * t115;
t70 = -t103 * t87 + t107 * t106;
t124 = t63 * t105 - t70 * t50;
t114 = cos(qJ(4));
t112 = t50 * t52;
t111 = t50 * t55;
t53 = sin(qJ(4));
t96 = t53 * t105;
t94 = t105 * t114;
t76 = mrSges(4,2) + (t98 * pkin(10) - mrSges(5,3)) * t50;
t47 = t121 * t113 + t115 * t85;
t46 = t63 * pkin(3);
t41 = t70 * t105 + t63 * t50;
t40 = t122 * t113 + t74 * t115;
t39 = t123 * t113 + t72 * t115;
t38 = t59 * pkin(3);
t37 = t58 * pkin(3);
t32 = -t63 * t114 - t47 * t96;
t29 = t65 * t105 + t59 * t50;
t28 = t64 * t105 + t58 * t50;
t27 = t47 * t114 - t124 * t53;
t26 = t124 * t114 + t47 * t53;
t22 = -t59 * t114 - t40 * t96;
t20 = -t58 * t114 - t39 * t96;
t16 = t27 * t55 + t41 * t52;
t14 = t40 * t114 - t125 * t53;
t13 = t125 * t114 + t40 * t53;
t12 = t39 * t114 - t126 * t53;
t11 = t126 * t114 + t39 * t53;
t4 = t14 * t55 + t29 * t52;
t2 = t12 * t55 + t28 * t52;
t1 = [(-m(2) - t134) * g(3) (-t86 * g(1) + t89 * g(2) - t107 * g(3)) * t134 (t63 * mrSges(4,1) + m(5) * t46 - t32 * mrSges(5,1) + t76 * t47 - t132 * (t47 * t112 + t32 * t55) + t133 * (t47 * t94 - t63 * t53) + t95 * (-t47 * t111 + t32 * t52) - t128 * (t32 * pkin(4) - t46)) * g(3) + (t58 * mrSges(4,1) + m(5) * t37 - t20 * mrSges(5,1) + t95 * (-t39 * t111 + t20 * t52) + t76 * t39 - t132 * (t39 * t112 + t20 * t55) + t133 * (t39 * t94 - t58 * t53) - t128 * (t20 * pkin(4) - t37)) * g(2) + (t59 * mrSges(4,1) + m(5) * t38 - t22 * mrSges(5,1) + t95 * (-t40 * t111 + t22 * t52) + t76 * t40 - t132 * (t40 * t112 + t22 * t55) + t133 * (t40 * t94 - t59 * t53) - t128 * (t22 * pkin(4) - t38)) * g(1) (t130 * t26 + t133 * t27) * g(3) + (t130 * t11 + t12 * t133) * g(2) + (t130 * t13 + t133 * t14) * g(1) (t95 * t16 - t132 * (-t27 * t52 + t41 * t55)) * g(3) + (t95 * t2 - t132 * (-t12 * t52 + t28 * t55)) * g(2) + (t95 * t4 - t132 * (-t14 * t52 + t29 * t55)) * g(1), -g(1) * ((t13 * t54 - t4 * t51) * mrSges(7,1) + (-t13 * t51 - t4 * t54) * mrSges(7,2)) - g(2) * ((t11 * t54 - t2 * t51) * mrSges(7,1) + (-t11 * t51 - t2 * t54) * mrSges(7,2)) - g(3) * ((-t16 * t51 + t26 * t54) * mrSges(7,1) + (-t16 * t54 - t26 * t51) * mrSges(7,2))];
taug  = t1(:);
