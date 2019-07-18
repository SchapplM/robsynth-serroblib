% Calculate Gravitation load on the joints for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:20
% EndTime: 2019-07-18 13:28:21
% DurationCPUTime: 0.24s
% Computational Cost: add. (128->44), mult. (195->55), div. (0->0), fcn. (162->8), ass. (0->32)
t13 = qJ(3) + qJ(4);
t11 = sin(t13);
t45 = (mrSges(5,2) - mrSges(6,3)) * t11;
t12 = cos(t13);
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t39 = m(5) + m(6);
t24 = pkin(2) * t39 + mrSges(4,1);
t22 = -mrSges(4,2) * t15 + t18 * t24;
t43 = mrSges(5,1) * t12 + mrSges(3,1) + t22 - t45;
t14 = sin(qJ(5));
t34 = mrSges(6,2) * t14;
t42 = mrSges(6,3) * t12 + t11 * t34;
t16 = sin(qJ(2));
t38 = t42 * t16;
t19 = cos(qJ(2));
t37 = t42 * t19;
t32 = t14 * t19;
t31 = t16 * t14;
t17 = cos(qJ(5));
t30 = t16 * t17;
t29 = t17 * t19;
t26 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t25 = mrSges(6,1) * t17 + mrSges(5,1);
t23 = mrSges(5,2) * t12 + t11 * t25;
t21 = (t25 - t34) * t12 - t45;
t20 = mrSges(4,2) * t18 + t15 * t24 + t23;
t4 = t12 * t29 + t31;
t3 = -t12 * t32 + t30;
t2 = -t12 * t30 + t32;
t1 = t12 * t31 + t29;
t5 = [(-m(2) - m(3) - m(4) - t39) * g(3), (-t4 * mrSges(6,1) - t3 * mrSges(6,2) + t26 * t16 - t43 * t19) * g(3) + (-mrSges(6,1) * t2 - mrSges(6,2) * t1 + t43 * t16 + t26 * t19) * g(1), (t21 + t22) * g(2) + (t16 * t20 - t38) * g(3) + (t19 * t20 - t37) * g(1), (t16 * t23 - t38) * g(3) + t21 * g(2) + (t19 * t23 - t37) * g(1), -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(3) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t14 + mrSges(6,2) * t17) * t11];
taug  = t5(:);
