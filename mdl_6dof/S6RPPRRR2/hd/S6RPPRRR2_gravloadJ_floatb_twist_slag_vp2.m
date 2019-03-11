% Calculate Gravitation load on the joints for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:19
% EndTime: 2019-03-09 02:20:20
% DurationCPUTime: 0.50s
% Computational Cost: add. (413->91), mult. (339->103), div. (0->0), fcn. (297->12), ass. (0->54)
t74 = mrSges(6,3) + mrSges(7,3);
t31 = cos(qJ(5));
t15 = pkin(5) * t31 + pkin(4);
t25 = qJ(5) + qJ(6);
t20 = sin(t25);
t21 = cos(t25);
t29 = sin(qJ(5));
t73 = -m(6) * pkin(4) - m(7) * t15 - mrSges(6,1) * t31 - mrSges(7,1) * t21 + mrSges(6,2) * t29 + mrSges(7,2) * t20;
t33 = -pkin(9) - pkin(8);
t72 = -m(6) * pkin(8) + m(7) * t33 - t74;
t65 = m(7) * pkin(5);
t71 = m(3) + m(4);
t70 = mrSges(6,1) + t65;
t69 = -m(5) - m(6) - m(7);
t23 = pkin(11) + qJ(4);
t16 = sin(t23);
t27 = cos(pkin(11));
t18 = cos(t23);
t42 = t18 * mrSges(5,1) - t16 * mrSges(5,2);
t67 = m(4) * pkin(2) + t27 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + mrSges(3,1) + t42 + t74 * t16;
t66 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t24 = qJ(1) + pkin(10);
t19 = cos(t24);
t51 = t19 * t21;
t17 = sin(t24);
t56 = t17 * t20;
t5 = t18 * t56 + t51;
t52 = t19 * t20;
t55 = t17 * t21;
t6 = -t18 * t55 + t52;
t64 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t7 = -t18 * t52 + t55;
t8 = t18 * t51 + t56;
t63 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t60 = g(3) * t16;
t30 = sin(qJ(1));
t59 = t30 * pkin(1);
t32 = cos(qJ(1));
t22 = t32 * pkin(1);
t54 = t17 * t29;
t53 = t17 * t31;
t50 = t19 * t29;
t49 = t19 * t31;
t48 = m(4) - t69;
t44 = pkin(4) * t18 + pkin(8) * t16;
t40 = -mrSges(7,1) * t20 - mrSges(7,2) * t21;
t39 = t15 * t18 - t16 * t33;
t11 = -t18 * t50 + t53;
t9 = t18 * t54 + t49;
t28 = -pkin(7) - qJ(3);
t14 = pkin(3) * t27 + pkin(2);
t12 = t18 * t49 + t54;
t10 = -t18 * t53 + t50;
t1 = [(-t54 * t65 - t32 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t30 * mrSges(2,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + t69 * (t14 * t19 - t17 * t28 + t22) - t71 * t22 + t66 * t17 + (-m(6) * t44 - m(7) * t39 - t67) * t19) * g(2) + (-t50 * t65 + t30 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t32 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t71 * t59 + t69 * (-t19 * t28 - t59) + t66 * t19 + (m(5) * t14 - m(6) * (-t14 - t44) - m(7) * (-t14 - t39) + t67) * t17) * g(1) (-m(3) - t48) * g(3) (-g(1) * t17 + g(2) * t19) * t48 (t16 * t72 + t18 * t73 - t42) * g(3) + (g(1) * t19 + g(2) * t17) * ((mrSges(5,2) + t72) * t18 + (mrSges(5,1) - t73) * t16) (mrSges(6,2) * t31 + t29 * t70 - t40) * t60 + (-t10 * mrSges(6,2) + t70 * t9 - t64) * g(2) + (t12 * mrSges(6,2) - t11 * t70 - t63) * g(1), -g(1) * t63 - g(2) * t64 - t40 * t60];
taug  = t1(:);
