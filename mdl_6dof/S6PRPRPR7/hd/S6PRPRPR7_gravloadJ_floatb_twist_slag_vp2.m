% Calculate Gravitation load on the joints for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:13
% EndTime: 2019-03-08 19:51:15
% DurationCPUTime: 0.81s
% Computational Cost: add. (318->79), mult. (791->109), div. (0->0), fcn. (890->10), ass. (0->42)
t69 = m(6) + m(7);
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t61 = t29 * mrSges(7,1) + t32 * mrSges(7,2) + t69 * qJ(5) - mrSges(5,2) + mrSges(6,3);
t80 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t79 = t80 * t30 + t61 * t33 + mrSges(3,2) - mrSges(4,3);
t50 = m(4) + m(5) + t69;
t71 = pkin(4) * t69 - t80;
t72 = t32 * mrSges(7,1) - t29 * mrSges(7,2) + mrSges(3,1) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t63 = -m(7) * (-pkin(5) - pkin(8)) + t72;
t60 = -t50 * qJ(3) + t79;
t59 = pkin(4) * t30;
t27 = sin(pkin(10));
t28 = sin(pkin(6));
t58 = t27 * t28;
t31 = sin(qJ(2));
t57 = t28 * t31;
t34 = cos(qJ(2));
t56 = t28 * t34;
t54 = pkin(2) * t56 + qJ(3) * t57;
t53 = cos(pkin(6));
t52 = cos(pkin(10));
t51 = pkin(8) * t56 + t54;
t42 = t53 * t52;
t14 = t27 * t31 - t34 * t42;
t11 = t14 * pkin(2);
t49 = -pkin(8) * t14 - t11;
t47 = t27 * t53;
t16 = t52 * t31 + t34 * t47;
t12 = t16 * pkin(2);
t48 = -pkin(8) * t16 - t12;
t46 = t28 * t52;
t18 = t53 * t30 + t33 * t56;
t17 = -t31 * t47 + t52 * t34;
t15 = t27 * t34 + t31 * t42;
t8 = t17 * t59;
t7 = t15 * t59;
t5 = t14 * t33 + t30 * t46;
t3 = -t16 * t33 + t30 * t58;
t1 = [(-m(2) - m(3) - t50) * g(3) (-m(4) * t54 - m(5) * t51 - t69 * (t57 * t59 + t51) + ((-m(7) * pkin(5) - t72) * t34 + t79 * t31) * t28) * g(3) + (m(4) * t11 - m(5) * t49 - m(6) * (t49 + t7) - m(7) * (-t11 + t7) + t60 * t15 + t63 * t14) * g(2) + (m(4) * t12 - m(5) * t48 - m(6) * (t48 + t8) - m(7) * (-t12 + t8) + t60 * t17 + t63 * t16) * g(1) (-g(1) * t16 - g(2) * t14 + g(3) * t56) * t50 (-t61 * (-t30 * t56 + t53 * t33) + t71 * t18) * g(3) + (t61 * (-t14 * t30 + t33 * t46) - t71 * t5) * g(2) + (-t61 * (t16 * t30 + t33 * t58) + t71 * t3) * g(1), t69 * (-g(1) * t3 + g(2) * t5 - g(3) * t18) -g(1) * ((-t17 * t29 + t3 * t32) * mrSges(7,1) + (-t17 * t32 - t29 * t3) * mrSges(7,2)) - g(2) * ((-t15 * t29 - t32 * t5) * mrSges(7,1) + (-t15 * t32 + t29 * t5) * mrSges(7,2)) - g(3) * ((t18 * t32 - t29 * t57) * mrSges(7,1) + (-t18 * t29 - t32 * t57) * mrSges(7,2))];
taug  = t1(:);
