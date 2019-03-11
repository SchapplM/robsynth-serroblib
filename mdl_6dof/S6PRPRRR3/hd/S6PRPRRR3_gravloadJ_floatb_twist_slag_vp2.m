% Calculate Gravitation load on the joints for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:47
% EndTime: 2019-03-08 20:31:49
% DurationCPUTime: 0.75s
% Computational Cost: add. (620->98), mult. (830->137), div. (0->0), fcn. (928->14), ass. (0->54)
t109 = mrSges(6,2) - mrSges(7,3);
t51 = sin(qJ(6));
t53 = cos(qJ(6));
t108 = -t53 * mrSges(7,1) + t51 * mrSges(7,2) - mrSges(6,1);
t45 = pkin(12) + qJ(4);
t40 = sin(t45);
t41 = cos(t45);
t82 = cos(pkin(6));
t48 = sin(pkin(6));
t52 = sin(qJ(2));
t89 = t48 * t52;
t107 = -t40 * t89 + t82 * t41;
t54 = cos(qJ(2));
t47 = sin(pkin(11));
t71 = t47 * t82;
t81 = cos(pkin(11));
t30 = -t52 * t71 + t81 * t54;
t90 = t47 * t48;
t106 = -t30 * t40 + t41 * t90;
t50 = -pkin(8) - qJ(3);
t95 = -m(4) * qJ(3) + m(5) * t50 - t51 * mrSges(7,1) - t53 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t42 = qJ(5) + t45;
t37 = sin(t42);
t38 = cos(t42);
t49 = cos(pkin(12));
t39 = t49 * pkin(3) + pkin(2);
t104 = -m(4) * pkin(2) - m(5) * t39 - t49 * mrSges(4,1) - t41 * mrSges(5,1) + sin(pkin(12)) * mrSges(4,2) + t40 * mrSges(5,2) - mrSges(3,1) + (-m(7) * pkin(5) + t108) * t38 + (-m(7) * pkin(10) + t109) * t37;
t103 = m(6) + m(7);
t22 = -t37 * t89 + t82 * t38;
t23 = t82 * t37 + t38 * t89;
t101 = t108 * t22 + t109 * t23;
t13 = -t30 * t37 + t38 * t90;
t14 = t30 * t38 + t37 * t90;
t100 = t108 * t13 + t109 * t14;
t62 = t82 * t81;
t28 = t47 * t54 + t52 * t62;
t70 = t48 * t81;
t11 = -t28 * t37 - t38 * t70;
t12 = t28 * t38 - t37 * t70;
t99 = t108 * t11 + t109 * t12;
t88 = t48 * t54;
t78 = m(4) + m(5) + t103;
t76 = t11 * pkin(5) + t12 * pkin(10);
t74 = t13 * pkin(5) + pkin(10) * t14;
t73 = t22 * pkin(5) + pkin(10) * t23;
t68 = t106 * pkin(4);
t63 = t107 * pkin(4);
t57 = -t28 * t40 - t41 * t70;
t56 = t57 * pkin(4);
t44 = -pkin(9) + t50;
t32 = pkin(4) * t41 + t39;
t29 = t81 * t52 + t54 * t71;
t27 = t47 * t52 - t54 * t62;
t1 = [(-m(2) - m(3) - t78) * g(3) (-t103 * (-t27 * t32 - t28 * t44) + t95 * t28 - t104 * t27) * g(2) + (-t103 * (-t29 * t32 - t30 * t44) + t95 * t30 - t104 * t29) * g(1) + (-t103 * t32 * t88 + (t104 * t54 + (t103 * t44 + t95) * t52) * t48) * g(3) (-g(1) * t29 - g(2) * t27 + g(3) * t88) * t78 (-t107 * mrSges(5,1) - (-t82 * t40 - t41 * t89) * mrSges(5,2) - m(6) * t63 - m(7) * (t63 + t73) + t101) * g(3) + (-t57 * mrSges(5,1) - (-t28 * t41 + t40 * t70) * mrSges(5,2) - m(6) * t56 - m(7) * (t56 + t76) + t99) * g(2) + (-t106 * mrSges(5,1) - (-t30 * t41 - t40 * t90) * mrSges(5,2) - m(6) * t68 - m(7) * (t68 + t74) + t100) * g(1) (-m(7) * t73 + t101) * g(3) + (-m(7) * t76 + t99) * g(2) + (-m(7) * t74 + t100) * g(1), -g(1) * ((-t14 * t51 + t29 * t53) * mrSges(7,1) + (-t14 * t53 - t29 * t51) * mrSges(7,2)) - g(2) * ((-t12 * t51 + t27 * t53) * mrSges(7,1) + (-t12 * t53 - t27 * t51) * mrSges(7,2)) - g(3) * ((-t23 * t51 - t53 * t88) * mrSges(7,1) + (-t23 * t53 + t51 * t88) * mrSges(7,2))];
taug  = t1(:);
