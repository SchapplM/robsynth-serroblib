% Calculate Gravitation load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:22:26
% DurationCPUTime: 0.26s
% Computational Cost: add. (199->48), mult. (165->47), div. (0->0), fcn. (121->10), ass. (0->27)
t16 = qJ(3) + pkin(9);
t10 = cos(t16);
t19 = sin(qJ(3));
t12 = qJ(5) + t16;
t5 = sin(t12);
t6 = cos(t12);
t30 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t8 = sin(t16);
t44 = t10 * mrSges(5,1) - t19 * mrSges(4,2) - t8 * mrSges(5,2) + t30;
t42 = m(5) + m(6);
t40 = m(3) + m(4) + t42;
t43 = pkin(1) * t40 + mrSges(2,1);
t17 = qJ(1) + pkin(8);
t11 = cos(t17);
t9 = sin(t17);
t41 = g(1) * t11 + g(2) * t9;
t21 = cos(qJ(3));
t13 = t21 * pkin(3);
t33 = pkin(4) * t10 + t13;
t39 = m(4) * pkin(2) + t21 * mrSges(4,1) + mrSges(3,1) + m(5) * (t13 + pkin(2)) + m(6) * (pkin(2) + t33) + t44;
t18 = -qJ(4) - pkin(6);
t38 = -m(4) * pkin(6) + m(5) * t18 - m(6) * (pkin(7) - t18) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t31 = m(5) * pkin(3) + mrSges(4,1);
t28 = mrSges(6,1) * t5 + mrSges(6,2) * t6;
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = [(t20 * mrSges(2,2) - t39 * t11 - t43 * t22 + t38 * t9) * g(2) + (t22 * mrSges(2,2) + t38 * t11 + t43 * t20 + t39 * t9) * g(1), -t40 * g(3), (-m(6) * t33 - t21 * t31 - t44) * g(3) + t41 * (-m(6) * (-t19 * pkin(3) - pkin(4) * t8) + mrSges(5,1) * t8 + mrSges(4,2) * t21 + mrSges(5,2) * t10 + t31 * t19 + t28), t42 * (-g(1) * t9 + g(2) * t11), -g(3) * t30 + t41 * t28];
taug = t1(:);
