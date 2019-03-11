% Calculate Gravitation load on the joints for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:41:07
% EndTime: 2019-03-08 20:41:09
% DurationCPUTime: 0.82s
% Computational Cost: add. (463->100), mult. (855->141), div. (0->0), fcn. (961->12), ass. (0->57)
t110 = mrSges(6,2) - mrSges(7,3);
t48 = sin(qJ(6));
t51 = cos(qJ(6));
t109 = -t51 * mrSges(7,1) + t48 * mrSges(7,2) - mrSges(6,1);
t102 = m(6) + m(7);
t103 = m(4) + m(5);
t112 = -t103 - t102;
t44 = qJ(4) + qJ(5);
t42 = sin(t44);
t43 = cos(t44);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t111 = -t49 * mrSges(5,1) - t52 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(7) * pkin(5) + t109) * t42 + (m(7) * pkin(10) - t110) * t43;
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t79 = cos(pkin(6));
t68 = t53 * t79;
t31 = t45 * t50 - t47 * t68;
t46 = sin(pkin(6));
t88 = t46 * t49;
t107 = t31 * t52 + t47 * t88;
t33 = t45 * t68 + t47 * t50;
t106 = t33 * t52 - t45 * t88;
t85 = t46 * t53;
t25 = -t42 * t79 - t43 * t85;
t26 = -t42 * t85 + t43 * t79;
t100 = t109 * t25 + t110 * t26;
t89 = t46 * t47;
t13 = t31 * t43 + t42 * t89;
t14 = -t31 * t42 + t43 * t89;
t99 = t109 * t13 - t110 * t14;
t90 = t45 * t46;
t11 = t33 * t43 - t42 * t90;
t12 = t33 * t42 + t43 * t90;
t98 = t109 * t11 + t110 * t12;
t97 = m(5) * pkin(8) + t48 * mrSges(7,1) + t51 * mrSges(7,2) + mrSges(3,1) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t96 = t112 * qJ(3) + t111;
t93 = pkin(4) * t49;
t87 = t46 * t50;
t86 = t46 * t52;
t81 = t107 * pkin(4);
t80 = pkin(2) * t85 + qJ(3) * t87;
t72 = t11 * pkin(5) + pkin(10) * t12;
t71 = t13 * pkin(5) - pkin(10) * t14;
t70 = t25 * pkin(5) + t26 * pkin(10);
t69 = t50 * t79;
t62 = t106 * pkin(4);
t58 = -t49 * t79 - t52 * t85;
t56 = t58 * pkin(4);
t54 = -pkin(9) - pkin(8);
t34 = -t45 * t69 + t47 * t53;
t32 = t45 * t53 + t47 * t69;
t30 = t33 * pkin(2);
t29 = t31 * pkin(2);
t1 = [(-m(2) - m(3) + t112) * g(3) (-t103 * t80 - t102 * (t87 * t93 + t80) + ((t102 * t54 - t97) * t53 + t111 * t50) * t46) * g(3) + (-t102 * (t31 * t54 + t32 * t93 - t29) + t103 * t29 + t96 * t32 + t97 * t31) * g(2) + (-t102 * (t33 * t54 + t34 * t93 - t30) + t103 * t30 + t96 * t34 + t97 * t33) * g(1) -(-g(1) * t33 - g(2) * t31 + g(3) * t85) * t112 (-t58 * mrSges(5,1) - (t49 * t85 - t52 * t79) * mrSges(5,2) - m(6) * t56 - m(7) * (t56 + t70) + t100) * g(3) + (-t107 * mrSges(5,1) - (-t31 * t49 + t47 * t86) * mrSges(5,2) - m(6) * t81 - m(7) * (t71 + t81) + t99) * g(2) + (-t106 * mrSges(5,1) - (-t33 * t49 - t45 * t86) * mrSges(5,2) - m(6) * t62 - m(7) * (t62 + t72) + t98) * g(1) (-m(7) * t70 + t100) * g(3) + (-m(7) * t71 + t99) * g(2) + (-m(7) * t72 + t98) * g(1), -g(1) * ((-t12 * t48 + t34 * t51) * mrSges(7,1) + (-t12 * t51 - t34 * t48) * mrSges(7,2)) - g(2) * ((t14 * t48 + t32 * t51) * mrSges(7,1) + (t14 * t51 - t32 * t48) * mrSges(7,2)) - g(3) * ((-t26 * t48 + t51 * t87) * mrSges(7,1) + (-t26 * t51 - t48 * t87) * mrSges(7,2))];
taug  = t1(:);
