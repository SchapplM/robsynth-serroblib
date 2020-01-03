% Calculate Gravitation load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:07:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (103->61), mult. (225->76), div. (0->0), fcn. (203->6), ass. (0->35)
t49 = m(4) + m(5);
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t48 = (-mrSges(3,1) - mrSges(4,1)) * t21 + (mrSges(3,2) - mrSges(4,3)) * t18;
t47 = -mrSges(2,1) + t48;
t46 = pkin(6) * m(5) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t44 = g(1) * t22 + g(2) * t19;
t43 = -pkin(2) - pkin(3);
t14 = t21 * pkin(2);
t20 = cos(qJ(4));
t40 = t18 * t20;
t39 = t22 * t21;
t12 = t18 * qJ(3);
t38 = t14 + t12;
t37 = t22 * pkin(1) + t19 * pkin(5);
t36 = qJ(3) * t21;
t35 = t18 * t43;
t34 = -pkin(1) - t12;
t33 = pkin(2) * t39 + t22 * t12 + t37;
t17 = sin(qJ(4));
t25 = t17 * t21 - t40;
t1 = t25 * t19;
t24 = t17 * t18 + t20 * t21;
t2 = t24 * t19;
t32 = -t1 * mrSges(5,1) - t2 * mrSges(5,2);
t3 = t17 * t39 - t22 * t40;
t4 = t24 * t22;
t31 = -t3 * mrSges(5,1) - t4 * mrSges(5,2);
t30 = -mrSges(5,1) * t24 + mrSges(5,2) * t25;
t23 = t21 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t18;
t9 = t22 * t36;
t8 = t19 * t36;
t5 = [(-m(3) * t37 - m(4) * t33 - m(5) * (pkin(3) * t39 + t33) - t4 * mrSges(5,1) + t3 * mrSges(5,2) + t47 * t22 + t46 * t19) * g(2) + (t2 * mrSges(5,1) - t1 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * (t34 - t14) - m(5) * (t43 * t21 + t34) - t47) * t19 + ((-m(3) - t49) * pkin(5) + t46) * t22) * g(1), t44 * (mrSges(3,1) * t18 + mrSges(3,2) * t21) + (-m(4) * t8 - t23 * t19 - m(5) * (t19 * t35 + t8) + t32) * g(2) + (-m(4) * t9 - t23 * t22 - m(5) * (t22 * t35 + t9) + t31) * g(1) + (-m(4) * t38 - m(5) * (pkin(3) * t21 + t38) + t30 + t48) * g(3), (g(3) * t21 - t18 * t44) * t49, -g(1) * t31 - g(2) * t32 - g(3) * t30];
taug = t5(:);
