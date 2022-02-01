% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:19
% DurationCPUTime: 0.47s
% Computational Cost: add. (216->81), mult. (262->94), div. (0->0), fcn. (239->10), ass. (0->48)
t30 = sin(pkin(8));
t32 = cos(pkin(8));
t41 = qJ(3) + pkin(6);
t31 = cos(pkin(9));
t57 = t31 * pkin(3) + pkin(2);
t71 = m(4) * (pkin(2) * t32 + pkin(1)) + m(5) * (t57 * t32 + pkin(1)) + mrSges(2,1) + t32 * mrSges(3,1) + (m(4) * qJ(3) + m(5) * t41 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t30;
t70 = -m(3) - m(6);
t68 = -m(4) - m(6);
t67 = m(6) * pkin(4) + mrSges(5,1);
t29 = sin(pkin(9));
t58 = t29 * pkin(3);
t28 = pkin(9) + qJ(4);
t20 = sin(t28);
t60 = pkin(4) * t20;
t65 = -m(5) * (qJ(2) + t58) - m(6) * (t58 + t60) + mrSges(2,2) - mrSges(3,3);
t22 = qJ(5) + t28;
t19 = cos(t22);
t34 = cos(qJ(1));
t47 = t34 * t19;
t18 = sin(t22);
t33 = sin(qJ(1));
t55 = t33 * t18;
t5 = t32 * t55 + t47;
t48 = t34 * t18;
t54 = t33 * t19;
t6 = -t32 * t54 + t48;
t62 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t32 * t48 + t54;
t8 = t32 * t47 + t55;
t61 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t59 = g(3) * t30;
t53 = t33 * t20;
t21 = cos(t28);
t52 = t33 * t21;
t51 = t33 * t29;
t49 = t33 * t31;
t46 = t34 * t20;
t45 = t34 * t21;
t44 = t34 * t29;
t42 = t34 * t31;
t39 = m(5) - t68;
t36 = -mrSges(6,1) * t18 - mrSges(6,2) * t19;
t35 = (pkin(4) * t21 + t57) * t32 - (-pkin(7) - t41) * t30;
t11 = -t32 * t46 + t52;
t9 = t32 * t53 + t45;
t12 = t32 * t45 + t53;
t10 = -t32 * t52 + t46;
t1 = [(-(t32 * t42 + t51) * mrSges(4,1) - (-t32 * t44 + t49) * mrSges(4,2) - t12 * mrSges(5,1) - t11 * mrSges(5,2) - t8 * mrSges(6,1) - t7 * mrSges(6,2) + (t65 + (-m(4) + t70) * qJ(2)) * t33) * g(2) + (-(-t32 * t49 + t44) * mrSges(4,1) - (t32 * t51 + t42) * mrSges(4,2) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - t6 * mrSges(6,1) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(6) * (-pkin(1) - t35) + t71) * t33) * g(1) + ((-m(6) * t35 + t70 * pkin(1) - t71) * g(2) + ((-m(3) + t68) * qJ(2) + t65) * g(1)) * t34, (-g(1) * t33 + g(2) * t34) * (m(3) + t39), (g(3) * t32 + (-g(1) * t34 - g(2) * t33) * t30) * t39, (m(6) * t60 + mrSges(5,1) * t20 + mrSges(5,2) * t21 - t36) * t59 + (-t10 * mrSges(5,2) + t67 * t9 - t62) * g(2) + (t12 * mrSges(5,2) - t67 * t11 - t61) * g(1), -g(1) * t61 - g(2) * t62 - t36 * t59];
taug = t1(:);
