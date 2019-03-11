% Calculate Gravitation load on the joints for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:26
% DurationCPUTime: 0.53s
% Computational Cost: add. (298->81), mult. (362->93), div. (0->0), fcn. (315->10), ass. (0->48)
t76 = mrSges(5,2) - m(6) * pkin(8) + m(7) * (-pkin(9) - pkin(8)) - mrSges(6,3) - mrSges(7,3);
t24 = pkin(10) + qJ(4);
t17 = sin(t24);
t18 = cos(t24);
t75 = -t17 * mrSges(5,1) - t76 * t18;
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t74 = g(1) * t30 - g(2) * t32;
t31 = cos(qJ(5));
t73 = -m(6) * pkin(4) - m(7) * (t31 * pkin(5) + pkin(4));
t64 = m(7) * pkin(5);
t70 = mrSges(6,1) + t64;
t25 = qJ(5) + qJ(6);
t19 = sin(t25);
t20 = cos(t25);
t29 = sin(qJ(5));
t69 = t31 * mrSges(6,1) + t20 * mrSges(7,1) - t29 * mrSges(6,2) - t19 * mrSges(7,2) - t73;
t67 = -m(5) - m(6) - m(7);
t66 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t26 = sin(pkin(10));
t65 = -t26 * mrSges(4,1) - cos(pkin(10)) * mrSges(4,2) + mrSges(2,2) - mrSges(3,3) + t73 * t17 + t75;
t52 = t32 * t20;
t57 = t30 * t19;
t5 = -t17 * t57 + t52;
t53 = t32 * t19;
t56 = t30 * t20;
t6 = t17 * t56 + t53;
t63 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t17 * t53 + t56;
t8 = t17 * t52 - t57;
t62 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t61 = pkin(3) * t26;
t58 = g(3) * t18;
t55 = t30 * t29;
t54 = t30 * t31;
t51 = t32 * t29;
t50 = t32 * t31;
t49 = t32 * pkin(1) + t30 * qJ(2);
t48 = -m(4) + t67;
t22 = t32 * qJ(2);
t45 = -t30 * pkin(1) + t22;
t40 = -mrSges(7,1) * t19 - mrSges(7,2) * t20;
t11 = t17 * t51 + t54;
t9 = -t17 * t55 + t50;
t28 = -pkin(7) - qJ(3);
t12 = t17 * t50 - t55;
t10 = t17 * t54 + t51;
t1 = [(-t51 * t64 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (-m(3) - m(4)) * t49 + t67 * (-t32 * t28 + t30 * t61 + t49) + (-m(4) * qJ(3) - t66) * t32 + t65 * t30) * g(2) + (t55 * t64 - m(3) * t45 - m(4) * t22 - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + t67 * (t30 * t28 + t32 * t61 + t45) + (-m(4) * (-pkin(1) - qJ(3)) + t66) * t30 + t65 * t32) * g(1), -t74 * (m(3) - t48) (g(1) * t32 + g(2) * t30) * t48 (t69 * t17 - t75) * g(3) + t74 * ((-mrSges(5,1) - t69) * t18 + t76 * t17) (mrSges(6,2) * t31 + t70 * t29 - t40) * t58 + (-t12 * mrSges(6,2) - t70 * t11 - t62) * g(2) + (t10 * mrSges(6,2) - t70 * t9 - t63) * g(1), -g(1) * t63 - g(2) * t62 - t40 * t58];
taug  = t1(:);
