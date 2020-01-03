% Calculate Gravitation load on the joints for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:38:10
% DurationCPUTime: 0.59s
% Computational Cost: add. (196->86), mult. (351->102), div. (0->0), fcn. (322->8), ass. (0->45)
t65 = m(5) + m(6);
t63 = m(4) + t65;
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t70 = (mrSges(3,1) + mrSges(4,1)) * t30 + (-mrSges(3,2) + mrSges(4,3)) * t28;
t25 = sin(pkin(8));
t54 = t30 * t25;
t26 = cos(pkin(8));
t16 = t26 * pkin(4) + pkin(3);
t57 = -pkin(2) - t16;
t60 = -pkin(2) - pkin(3);
t69 = -m(6) * pkin(4) * t54 + (m(4) * pkin(2) - m(5) * t60 - m(6) * t57 + mrSges(4,1)) * t28 + (-qJ(3) * t63 - mrSges(4,3)) * t30;
t67 = mrSges(2,1) + t70;
t66 = qJ(4) * m(5) - (-pkin(7) - qJ(4)) * m(6) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t64 = g(1) * t31 + g(2) * t29;
t21 = t30 * pkin(2);
t24 = pkin(8) + qJ(5);
t18 = cos(t24);
t56 = t28 * t18;
t55 = t28 * t25;
t53 = t31 * t30;
t19 = t28 * qJ(3);
t52 = t21 + t19;
t51 = t31 * pkin(1) + t29 * pkin(6);
t47 = -pkin(1) - t19;
t17 = sin(t24);
t38 = t30 * t17 - t56;
t1 = t38 * t29;
t37 = t28 * t17 + t30 * t18;
t2 = t37 * t29;
t45 = -t1 * mrSges(6,1) - t2 * mrSges(6,2);
t3 = t17 * t53 - t31 * t56;
t4 = t37 * t31;
t44 = -t3 * mrSges(6,1) - t4 * mrSges(6,2);
t42 = -mrSges(6,1) * t37 + mrSges(6,2) * t38;
t36 = -t28 * t26 + t54;
t35 = t30 * t26 + t55;
t34 = pkin(4) * t55 + t30 * t16;
t8 = t35 * t31;
t7 = t36 * t31;
t6 = t35 * t29;
t5 = t36 * t29;
t9 = [(-m(5) * pkin(3) * t53 - m(3) * t51 - t8 * mrSges(5,1) - t4 * mrSges(6,1) + t7 * mrSges(5,2) + t3 * mrSges(6,2) - t63 * (pkin(2) * t53 + t31 * t19 + t51) + (-m(6) * t34 - t67) * t31 + t66 * t29) * g(2) + (t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t1 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * (t47 - t21) - m(5) * (t30 * t60 + t47) - m(6) * (-pkin(1) + t57 * t30 + (-pkin(4) * t25 - qJ(3)) * t28) + t67) * t29 + ((-m(3) - t63) * pkin(6) + t66) * t31) * g(1), t64 * (mrSges(3,1) * t28 + mrSges(3,2) * t30) + (-t5 * mrSges(5,1) - t6 * mrSges(5,2) + t69 * t29 + t45) * g(2) + (-t7 * mrSges(5,1) - t8 * mrSges(5,2) + t69 * t31 + t44) * g(1) + (-m(4) * t52 - m(5) * (t30 * pkin(3) + t52) - t35 * mrSges(5,1) + t36 * mrSges(5,2) - m(6) * (t34 + t52) + t42 - t70) * g(3), (g(3) * t30 - t64 * t28) * t63, t65 * (g(1) * t29 - g(2) * t31), -g(1) * t44 - g(2) * t45 - g(3) * t42];
taug = t9(:);
