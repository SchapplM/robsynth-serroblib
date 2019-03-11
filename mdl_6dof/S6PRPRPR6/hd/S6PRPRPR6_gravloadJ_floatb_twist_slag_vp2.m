% Calculate Gravitation load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:42
% EndTime: 2019-03-08 19:46:43
% DurationCPUTime: 0.64s
% Computational Cost: add. (365->70), mult. (822->99), div. (0->0), fcn. (929->12), ass. (0->40)
t70 = m(6) + m(7);
t76 = m(5) + t70;
t52 = m(4) + t76;
t64 = mrSges(5,2) - m(6) * qJ(5) + m(7) * (-pkin(9) - qJ(5)) - mrSges(6,3) - mrSges(7,3);
t24 = pkin(11) + qJ(6);
t22 = sin(t24);
t23 = cos(t24);
t25 = sin(pkin(11));
t28 = cos(pkin(11));
t65 = mrSges(5,1) + m(7) * (pkin(5) * t28 + pkin(4)) + t23 * mrSges(7,1) - t22 * mrSges(7,2) + m(6) * pkin(4) + t28 * mrSges(6,1) - t25 * mrSges(6,2);
t75 = -t22 * mrSges(7,1) - t28 * mrSges(6,2) - t23 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t25;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t74 = -t65 * t30 - t64 * t32 + mrSges(3,2) - mrSges(4,3);
t73 = pkin(2) * t52 + t76 * pkin(8) - t75;
t61 = -t52 * qJ(3) + t74;
t26 = sin(pkin(10));
t27 = sin(pkin(6));
t59 = t26 * t27;
t31 = sin(qJ(2));
t58 = t27 * t31;
t33 = cos(qJ(2));
t57 = t27 * t33;
t56 = pkin(2) * t57 + qJ(3) * t58;
t55 = cos(pkin(6));
t54 = cos(pkin(10));
t49 = t26 * t55;
t48 = t27 * t54;
t44 = t55 * t54;
t14 = -t30 * t57 + t55 * t32;
t13 = t55 * t30 + t32 * t57;
t12 = -t31 * t49 + t54 * t33;
t11 = t54 * t31 + t33 * t49;
t10 = t26 * t33 + t31 * t44;
t9 = t26 * t31 - t33 * t44;
t4 = -t9 * t30 + t32 * t48;
t3 = t30 * t48 + t9 * t32;
t2 = t11 * t30 + t32 * t59;
t1 = -t11 * t32 + t30 * t59;
t5 = [(-m(2) - m(3) - t52) * g(3) (-m(4) * t56 - t76 * (pkin(8) * t57 + t56) + (t74 * t31 + t75 * t33) * t27) * g(3) + (t61 * t10 + t73 * t9) * g(2) + (t73 * t11 + t61 * t12) * g(1) (-g(1) * t11 - g(2) * t9 + g(3) * t57) * t52 (t65 * t13 + t64 * t14) * g(3) + (-t65 * t3 - t64 * t4) * g(2) + (t65 * t1 + t64 * t2) * g(1), t70 * (-g(1) * t1 + g(2) * t3 - g(3) * t13) -g(1) * ((t12 * t23 - t2 * t22) * mrSges(7,1) + (-t12 * t22 - t2 * t23) * mrSges(7,2)) - g(2) * ((t10 * t23 + t22 * t4) * mrSges(7,1) + (-t10 * t22 + t23 * t4) * mrSges(7,2)) - g(3) * ((-t14 * t22 + t23 * t58) * mrSges(7,1) + (-t14 * t23 - t22 * t58) * mrSges(7,2))];
taug  = t5(:);
