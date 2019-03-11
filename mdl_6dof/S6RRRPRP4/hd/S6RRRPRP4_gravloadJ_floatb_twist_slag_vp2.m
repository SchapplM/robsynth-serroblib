% Calculate Gravitation load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:36
% EndTime: 2019-03-09 16:43:39
% DurationCPUTime: 0.96s
% Computational Cost: add. (472->115), mult. (596->133), div. (0->0), fcn. (531->8), ass. (0->64)
t122 = -m(6) - m(7);
t126 = -mrSges(6,3) - mrSges(7,2);
t128 = -t126 + t122 * (-pkin(3) - pkin(9));
t127 = -mrSges(6,1) - mrSges(7,1);
t109 = m(7) * pkin(5) - t127;
t41 = qJ(2) + qJ(3);
t37 = cos(t41);
t125 = t37 * t109;
t36 = sin(t41);
t124 = t128 * t36;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t117 = g(1) * t47 + g(2) * t44;
t123 = (-mrSges(4,1) + mrSges(5,2)) * t37 + (mrSges(4,2) - mrSges(5,3)) * t36;
t113 = -m(7) * qJ(6) - mrSges(7,3);
t45 = cos(qJ(5));
t96 = t37 * t45;
t120 = -t113 * t96 + t124;
t116 = m(5) - t122;
t86 = qJ(4) * t37;
t16 = t47 * t86;
t83 = mrSges(6,2) * t96;
t42 = sin(qJ(5));
t93 = t42 * t47;
t95 = t37 * t47;
t98 = t36 * t47;
t115 = -mrSges(5,2) * t98 - mrSges(5,3) * t95 + t122 * t16 - t93 * t125 - t47 * t83;
t14 = t44 * t86;
t94 = t42 * t44;
t99 = t36 * t44;
t114 = -mrSges(5,2) * t99 + t122 * t14 - t94 * t125 + (-t37 * mrSges(5,3) - t83) * t44;
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t64 = t46 * mrSges(3,1) - t43 * mrSges(3,2);
t112 = -m(3) * pkin(1) - mrSges(2,1) + t123 - t64;
t111 = -m(3) * pkin(7) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t100 = t36 * t42;
t91 = t45 * mrSges(7,3);
t110 = t123 + t126 * t37 + (-mrSges(6,2) * t45 + t91) * t36 + t127 * t100;
t105 = pkin(2) * t43;
t85 = qJ(6) * t45;
t108 = -m(5) * (-pkin(3) * t36 - t105) - m(7) * (-t37 * t85 - t105) + t37 * t91 + m(6) * t105 + t124;
t107 = -mrSges(6,2) - t113;
t102 = g(3) * t37;
t34 = t37 * pkin(3);
t39 = t46 * pkin(2);
t92 = t44 * t45;
t90 = t45 * t47;
t28 = t36 * qJ(4);
t87 = t34 + t28;
t33 = t37 * pkin(9);
t79 = t33 + t87;
t78 = t39 + t87;
t35 = t39 + pkin(1);
t48 = -pkin(8) - pkin(7);
t74 = t47 * t35 - t44 * t48;
t61 = mrSges(4,1) * t36 + mrSges(4,2) * t37;
t58 = pkin(3) * t95 + t47 * t28 + t74;
t53 = pkin(5) * t100 - t36 * t85 + t79;
t4 = -t36 * t94 + t90;
t3 = t36 * t92 + t93;
t2 = t36 * t93 + t92;
t1 = -t36 * t90 + t94;
t5 = [(-m(4) * t74 - m(5) * t58 + t126 * t95 + t122 * (t44 * pkin(4) + pkin(9) * t95 + t58) - t109 * t2 - t107 * t1 + t111 * t44) * g(2) + (-t109 * t4 - t107 * t3 + (m(4) * t35 + m(5) * t34 - t116 * (-t35 - t28) + t128 * t37 - t112) * t44) * g(1) + (t112 * g(2) + (t122 * (pkin(4) - t48) + (m(4) + m(5)) * t48 + t111) * g(1)) * t47 (-m(5) * t14 + t108 * t44 + t114) * g(2) + (-m(5) * t16 + t108 * t47 + t115) * g(1) + (-t64 - m(4) * t39 - m(5) * t78 - m(6) * (t33 + t78) - m(7) * (t39 + t53) + t110) * g(3) + t117 * (m(4) * t105 + mrSges(3,1) * t43 + mrSges(3,2) * t46 + t61) t117 * t61 + (-m(5) * (-pkin(3) * t99 + t14) + t120 * t44 + t114) * g(2) + (-m(5) * (-pkin(3) * t98 + t16) + t120 * t47 + t115) * g(1) + (-m(5) * t87 - m(6) * t79 - m(7) * t53 + t110) * g(3) (-t117 * t36 + t102) * t116 (t107 * t42 + t109 * t45) * t102 + (t107 * t4 - t109 * t3) * g(2) + (t109 * t1 - t107 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - g(3) * t96) * m(7)];
taug  = t5(:);
