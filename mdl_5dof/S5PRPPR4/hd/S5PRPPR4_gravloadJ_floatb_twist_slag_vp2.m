% Calculate Gravitation load on the joints for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:47
% DurationCPUTime: 0.30s
% Computational Cost: add. (142->45), mult. (156->55), div. (0->0), fcn. (132->6), ass. (0->23)
t28 = m(5) + m(6);
t24 = m(4) + t28;
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t30 = -mrSges(3,1) + (-mrSges(4,1) - mrSges(5,1)) * t14 + (mrSges(4,2) - mrSges(5,3)) * t13;
t29 = m(6) * pkin(6) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t12 = pkin(7) + qJ(2);
t10 = sin(t12);
t11 = cos(t12);
t27 = t11 * pkin(2) + t10 * qJ(3);
t26 = t11 * t14;
t25 = qJ(4) * t13;
t23 = pkin(3) * t26 + t11 * t25 + t27;
t22 = -pkin(2) - t25;
t15 = sin(qJ(5));
t16 = cos(qJ(5));
t18 = t13 * t16 - t14 * t15;
t17 = t13 * t15 + t14 * t16;
t4 = t17 * t11;
t3 = t18 * t11;
t2 = t17 * t10;
t1 = t18 * t10;
t5 = [(-m(2) - m(3) - t24) * g(3), (-m(4) * t27 - m(5) * t23 - m(6) * (pkin(4) * t26 + t23) - t4 * mrSges(6,1) - t3 * mrSges(6,2) + t30 * t11 + t29 * t10) * g(2) + (t2 * mrSges(6,1) + t1 * mrSges(6,2) + (m(4) * pkin(2) - m(5) * (-pkin(3) * t14 + t22) - m(6) * ((-pkin(3) - pkin(4)) * t14 + t22) - t30) * t10 + (-t24 * qJ(3) + t29) * t11) * g(1), (-g(1) * t10 + g(2) * t11) * t24, (t14 * g(3) + t13 * (-g(1) * t11 - g(2) * t10)) * t28, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(3) * (-t17 * mrSges(6,1) - t18 * mrSges(6,2))];
taug = t5(:);
