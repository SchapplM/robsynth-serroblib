% Calculate Gravitation load on the joints for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:07
% EndTime: 2019-12-05 16:30:10
% DurationCPUTime: 0.61s
% Computational Cost: add. (310->62), mult. (693->92), div. (0->0), fcn. (797->12), ass. (0->38)
t23 = pkin(10) + qJ(5);
t21 = sin(t23);
t22 = cos(t23);
t24 = sin(pkin(10));
t26 = cos(pkin(10));
t59 = mrSges(4,1) + m(6) * (pkin(4) * t26 + pkin(3)) + t22 * mrSges(6,1) - t21 * mrSges(6,2) + m(5) * pkin(3) + t26 * mrSges(5,1) - t24 * mrSges(5,2);
t58 = mrSges(4,2) - m(5) * qJ(4) - mrSges(5,3) + m(6) * (-pkin(8) - qJ(4)) - mrSges(6,3);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t71 = t28 * t58 - t30 * t59 - mrSges(3,1);
t65 = m(5) + m(6);
t60 = m(4) + t65;
t68 = -t21 * mrSges(6,1) - t26 * mrSges(5,2) - t22 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t24;
t67 = pkin(2) * t60 - t71;
t55 = -pkin(7) * t60 + t68;
t25 = sin(pkin(5));
t29 = sin(qJ(2));
t53 = t25 * t29;
t31 = cos(qJ(2));
t52 = t25 * t31;
t50 = cos(pkin(5));
t49 = cos(pkin(9));
t48 = sin(pkin(9));
t45 = t25 * t49;
t44 = t25 * t48;
t40 = t50 * t49;
t39 = t50 * t48;
t12 = t28 * t50 + t30 * t53;
t11 = t28 * t53 - t30 * t50;
t10 = -t29 * t39 + t31 * t49;
t9 = t29 * t49 + t31 * t39;
t8 = t29 * t40 + t31 * t48;
t7 = t29 * t48 - t31 * t40;
t4 = t10 * t30 + t28 * t44;
t3 = t10 * t28 - t30 * t44;
t2 = -t28 * t45 + t30 * t8;
t1 = t28 * t8 + t30 * t45;
t5 = [(-m(2) - m(3) - t60) * g(3), (-t60 * (pkin(2) * t52 + pkin(7) * t53) + (t68 * t29 + t31 * t71) * t25) * g(3) + (t55 * t8 + t67 * t7) * g(2) + (t55 * t10 + t67 * t9) * g(1), (t11 * t59 + t12 * t58) * g(3) + (t1 * t59 + t2 * t58) * g(2) + (t3 * t59 + t4 * t58) * g(1), t65 * (-g(1) * t3 - g(2) * t1 - g(3) * t11), -g(1) * ((-t21 * t4 + t22 * t9) * mrSges(6,1) + (-t21 * t9 - t22 * t4) * mrSges(6,2)) - g(2) * ((-t2 * t21 + t22 * t7) * mrSges(6,1) + (-t2 * t22 - t21 * t7) * mrSges(6,2)) - g(3) * ((-t12 * t21 - t22 * t52) * mrSges(6,1) + (-t12 * t22 + t21 * t52) * mrSges(6,2))];
taug = t5(:);
