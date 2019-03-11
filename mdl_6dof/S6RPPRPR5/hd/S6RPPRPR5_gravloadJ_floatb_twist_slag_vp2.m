% Calculate Gravitation load on the joints for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:32
% EndTime: 2019-03-09 01:48:32
% DurationCPUTime: 0.47s
% Computational Cost: add. (180->67), mult. (305->76), div. (0->0), fcn. (253->8), ass. (0->31)
t13 = sin(pkin(9));
t14 = cos(pkin(9));
t59 = m(6) * pkin(4) + t14 * mrSges(6,1) - t13 * mrSges(6,2) + mrSges(5,1);
t58 = m(7) * (t14 * pkin(5) + pkin(4));
t56 = m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5) - mrSges(6,3);
t16 = sin(qJ(4));
t18 = cos(qJ(4));
t57 = t18 * mrSges(5,2) + t59 * t16;
t51 = m(6) + m(7);
t35 = -t51 - m(4) - m(5);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t49 = g(1) * t19 + g(2) * t17;
t55 = -t16 * t58 - t56 * t18 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - t57;
t53 = -m(4) - m(6);
t52 = -m(5) - m(7);
t46 = m(6) * pkin(7) + t14 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t13;
t41 = t19 * pkin(1) + t17 * qJ(2);
t40 = t16 * t17;
t39 = t16 * t19;
t36 = t19 * qJ(3) + t41;
t12 = pkin(9) + qJ(6);
t6 = sin(t12);
t7 = cos(t12);
t25 = t7 * mrSges(7,1) - t6 * mrSges(7,2) + t58;
t10 = t19 * qJ(2);
t4 = -t17 * t6 + t7 * t39;
t3 = -t17 * t7 - t6 * t39;
t2 = -t19 * t6 - t7 * t40;
t1 = -t19 * t7 + t6 * t40;
t5 = [(-m(3) * t41 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + t53 * t36 + t52 * (-t17 * pkin(7) + t36) + t46 * t17 + t55 * t19) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + t52 * (-t19 * pkin(7) + t10) + (-m(3) + t53) * t10 + t46 * t19 + (m(3) * pkin(1) + t35 * (-pkin(1) - qJ(3)) - t55) * t17) * g(1) (-g(1) * t17 + g(2) * t19) * (m(3) - t35) t49 * t35, t57 * g(3) + (t56 * g(3) + t49 * (-t25 - t59)) * t18 + (t25 * g(3) + t49 * (mrSges(5,2) + t56)) * t16 (-g(3) * t16 + t49 * t18) * t51, -g(1) * (t3 * mrSges(7,1) - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(3) * (-mrSges(7,1) * t6 - mrSges(7,2) * t7) * t18];
taug  = t5(:);
