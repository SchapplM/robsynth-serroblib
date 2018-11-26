% Calculate Gravitation load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:05
% EndTime: 2018-11-23 15:38:06
% DurationCPUTime: 0.38s
% Computational Cost: add. (242->69), mult. (375->82), div. (0->0), fcn. (399->10), ass. (0->36)
t21 = cos(pkin(10));
t56 = -mrSges(4,1) - m(5) * pkin(3) - t21 * mrSges(5,1) + sin(pkin(10)) * mrSges(5,2);
t55 = -m(4) - m(5);
t54 = -m(6) - m(7);
t37 = m(5) - t54;
t53 = m(4) + t37;
t52 = mrSges(2,1) + mrSges(3,1);
t51 = mrSges(2,2) - mrSges(3,3);
t48 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t19 = pkin(10) + qJ(5);
t13 = sin(t19);
t14 = cos(t19);
t23 = sin(qJ(6));
t24 = cos(qJ(6));
t28 = m(7) * pkin(5) + t24 * mrSges(7,1) - t23 * mrSges(7,2);
t32 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t35 = m(7) * pkin(8) + mrSges(7,3);
t47 = t35 * t13 + t28 * t14 + t32;
t44 = cos(qJ(1));
t43 = sin(qJ(1));
t42 = t14 * t23;
t41 = t14 * t24;
t40 = t44 * pkin(1) + t43 * qJ(2);
t39 = cos(pkin(9));
t38 = sin(pkin(9));
t36 = t44 * pkin(2) + t40;
t33 = -t43 * pkin(1) + t44 * qJ(2);
t30 = mrSges(7,1) * t23 + mrSges(7,2) * t24;
t27 = -t43 * pkin(2) + t33;
t22 = -pkin(7) - qJ(4);
t12 = t21 * pkin(4) + pkin(3);
t8 = t44 * t38 - t43 * t39;
t7 = -t43 * t38 - t44 * t39;
t2 = t8 * t23 - t7 * t41;
t1 = t8 * t24 + t7 * t42;
t3 = [(-m(3) * t40 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t52 * t44 + t51 * t43 + t55 * t36 + t54 * (-t7 * t12 - t8 * t22 + t36) + t48 * t8 + (t32 - m(7) * (-t14 * pkin(5) - t13 * pkin(8)) + t13 * mrSges(7,3) - t56) * t7) * g(2) + (-m(3) * t33 + t51 * t44 + t52 * t43 + t55 * t27 + t54 * (t8 * t12 - t7 * t22 + t27) + (-t47 + t56) * t8 + (-t30 + t48) * t7) * g(1) (-t43 * g(1) + t44 * g(2)) * (m(3) + t53) t53 * g(3) (-g(1) * t8 + g(2) * t7) * t37, t47 * g(3) + ((mrSges(6,2) - t35) * t14 + (mrSges(6,1) + t28) * t13) * (-g(1) * t7 - g(2) * t8) -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * ((-t7 * t24 + t8 * t42) * mrSges(7,1) + (t7 * t23 + t8 * t41) * mrSges(7,2)) - g(3) * t30 * t13];
taug  = t3(:);
