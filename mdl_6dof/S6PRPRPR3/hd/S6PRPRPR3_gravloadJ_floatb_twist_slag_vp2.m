% Calculate Gravitation load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:30
% EndTime: 2019-03-08 19:34:32
% DurationCPUTime: 0.80s
% Computational Cost: add. (520->96), mult. (1328->140), div. (0->0), fcn. (1624->12), ass. (0->56)
t97 = m(6) + m(7);
t103 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t45 = sin(qJ(6));
t48 = cos(qJ(6));
t102 = t45 * mrSges(7,1) + t48 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t40 = sin(pkin(11));
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t69 = cos(pkin(11));
t30 = -t50 * t40 - t47 * t69;
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t44 = cos(pkin(6));
t78 = t44 * t50;
t99 = -t41 * t47 + t43 * t78;
t98 = m(7) * pkin(9) + pkin(4) * t97 + t103;
t58 = -t47 * t40 + t50 * t69;
t53 = t44 * t58;
t14 = t41 * t30 + t43 * t53;
t46 = sin(qJ(4));
t70 = qJ(5) * t46;
t49 = cos(qJ(4));
t87 = t14 * t49;
t96 = pkin(4) * t87 + t14 * t70;
t17 = t30 * t43 - t41 * t53;
t86 = t17 * t49;
t95 = pkin(4) * t86 + t17 * t70;
t42 = sin(pkin(6));
t27 = t58 * t42;
t85 = t27 * t49;
t94 = pkin(4) * t85 + t27 * t70;
t72 = t30 * t44;
t18 = t41 * t72 + t43 * t58;
t13 = -t41 * t58 + t43 * t72;
t91 = -t97 * qJ(5) - t102;
t90 = t102 * t46 + t103 * t49 + mrSges(4,1);
t89 = m(7) * (pkin(5) + pkin(8)) + t48 * mrSges(7,1) - t45 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t82 = t42 * t46;
t81 = t42 * t49;
t79 = t44 * t47;
t37 = t42 * t50 * pkin(2);
t71 = t27 * pkin(3) + t37;
t67 = m(4) + m(5) + t97;
t64 = t99 * pkin(2);
t28 = t30 * t42;
t63 = -pkin(8) * t28 + t71;
t60 = t14 * pkin(3) + t64;
t59 = -t41 * t78 - t43 * t47;
t56 = t59 * pkin(2);
t54 = -pkin(8) * t13 + t60;
t52 = t17 * pkin(3) + t56;
t51 = pkin(8) * t18 + t52;
t20 = -t28 * t46 - t44 * t49;
t5 = t18 * t46 - t41 * t81;
t3 = -t13 * t46 + t43 * t81;
t1 = [(-m(2) - m(3) - t67) * g(3) (-(t50 * mrSges(3,1) - t47 * mrSges(3,2)) * t42 - m(4) * t37 - m(5) * t63 - m(6) * (t63 + t94) - m(7) * (pkin(9) * t85 + t71 + t94) + t89 * t28 - t90 * t27) * g(3) + (-t99 * mrSges(3,1) - (-t41 * t50 - t43 * t79) * mrSges(3,2) - m(4) * t64 - m(5) * t54 - m(6) * (t54 + t96) - m(7) * (pkin(9) * t87 + t60 + t96) + t89 * t13 - t90 * t14) * g(2) + (-t59 * mrSges(3,1) - (t41 * t79 - t43 * t50) * mrSges(3,2) - m(4) * t56 - m(5) * t51 - m(6) * (t51 + t95) - m(7) * (pkin(9) * t86 + t52 + t95) - t89 * t18 - t90 * t17) * g(1) (-g(3) * t44 + (-g(1) * t41 + g(2) * t43) * t42) * t67 (t91 * (-t28 * t49 + t44 * t46) + t98 * t20) * g(3) + (t91 * (-t13 * t49 - t43 * t82) + t98 * t3) * g(2) + (t91 * (t18 * t49 + t41 * t82) + t98 * t5) * g(1), t97 * (-g(1) * t5 - g(2) * t3 - g(3) * t20) -g(1) * ((t17 * t45 + t48 * t5) * mrSges(7,1) + (t17 * t48 - t45 * t5) * mrSges(7,2)) - g(2) * ((t14 * t45 + t3 * t48) * mrSges(7,1) + (t14 * t48 - t3 * t45) * mrSges(7,2)) - g(3) * ((t20 * t48 + t27 * t45) * mrSges(7,1) + (-t20 * t45 + t27 * t48) * mrSges(7,2))];
taug  = t1(:);
