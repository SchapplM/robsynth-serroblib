% Calculate Gravitation load on the joints for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:43
% EndTime: 2019-03-09 01:50:44
% DurationCPUTime: 0.45s
% Computational Cost: add. (144->77), mult. (286->90), div. (0->0), fcn. (230->6), ass. (0->34)
t55 = m(6) + m(7);
t35 = -m(4) - m(5) - t55;
t21 = cos(qJ(1));
t58 = g(1) * t21;
t17 = sin(qJ(4));
t20 = cos(qJ(4));
t43 = t17 * mrSges(6,2);
t57 = -(m(7) * (-pkin(4) - pkin(8)) - mrSges(7,3)) * t17 - t20 * mrSges(6,3) - t43;
t16 = sin(qJ(6));
t19 = cos(qJ(6));
t28 = mrSges(7,1) * t16 + mrSges(7,2) * t19;
t54 = (-m(7) * pkin(8) + mrSges(6,2) - mrSges(7,3)) * t20 + (-t28 - mrSges(6,3)) * t17;
t18 = sin(qJ(1));
t53 = -g(2) * t18 - t58;
t29 = t17 * mrSges(5,1) + t20 * mrSges(5,2);
t52 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3) - t29;
t51 = m(7) * pkin(5) + mrSges(6,1) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t45 = g(3) * t17;
t44 = t17 * pkin(4);
t42 = t18 * t20;
t41 = t21 * t16;
t40 = t21 * t19;
t38 = pkin(1) * t21 + qJ(2) * t18;
t37 = qJ(5) * t17;
t12 = t20 * qJ(5);
t36 = qJ(3) * t21 + t38;
t14 = t21 * qJ(2);
t34 = -t21 * pkin(7) + t14;
t25 = -t18 * pkin(7) + t36;
t4 = t16 * t42 - t40;
t3 = t19 * t42 + t41;
t2 = t18 * t19 + t20 * t41;
t1 = t16 * t18 - t20 * t40;
t5 = [(-m(3) * t38 - m(4) * t36 - m(5) * t25 + t2 * mrSges(7,1) - t1 * mrSges(7,2) - t55 * (t21 * t44 + t25) + t51 * t18 + (t43 - (-m(6) * qJ(5) - mrSges(6,3)) * t20 - m(7) * (pkin(8) * t17 - t12) - t17 * mrSges(7,3) + t52) * t21) * g(2) + (-m(5) * t34 - t4 * mrSges(7,1) - t3 * mrSges(7,2) - t55 * (t18 * t12 + t34) + (-m(3) - m(4)) * t14 + t51 * t21 + (m(3) * pkin(1) + m(6) * t44 + t35 * (-pkin(1) - qJ(3)) - t52 + t57) * t18) * g(1) (-g(1) * t18 + g(2) * t21) * (m(3) - t35) -t53 * t35, t53 * (mrSges(5,1) * t20 - mrSges(5,2) * t17) + (-t55 * (pkin(4) * t42 + t18 * t37) + t54 * t18) * g(2) + (-t55 * (pkin(4) * t20 + t37) + t54) * t58 + (t29 - m(6) * (t12 - t44) - m(7) * t12 - t20 * t28 + t57) * g(3) (-t20 * t53 - t45) * t55, -g(1) * (mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(2) * (-mrSges(7,1) * t3 + mrSges(7,2) * t4) - (mrSges(7,1) * t19 - mrSges(7,2) * t16) * t45];
taug  = t5(:);
