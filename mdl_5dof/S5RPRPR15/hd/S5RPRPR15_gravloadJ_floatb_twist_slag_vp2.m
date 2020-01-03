% Calculate Gravitation load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:36:39
% DurationCPUTime: 0.42s
% Computational Cost: add. (155->55), mult. (265->64), div. (0->0), fcn. (227->8), ass. (0->27)
t52 = mrSges(4,2) + m(6) * (-pkin(7) - qJ(4)) - mrSges(6,3) - m(5) * qJ(4) - mrSges(5,3);
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t51 = t16 * mrSges(4,1) + t52 * t18;
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t49 = -m(6) * (pkin(4) * t14 + pkin(3)) - m(5) * pkin(3) - t14 * mrSges(5,1) + t13 * mrSges(5,2);
t43 = m(5) + m(6);
t47 = -m(4) - t43;
t48 = -m(3) + t47;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t40 = -g(1) * t17 + g(2) * t19;
t46 = t49 * t16 + mrSges(2,2) - mrSges(3,3) - t51;
t44 = t14 * mrSges(5,2) + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t13;
t12 = pkin(8) + qJ(5);
t6 = sin(t12);
t7 = cos(t12);
t42 = t7 * mrSges(6,1) - t6 * mrSges(6,2) - t49;
t34 = t19 * pkin(1) + t17 * qJ(2);
t33 = t16 * t17;
t32 = t16 * t19;
t4 = -t17 * t6 + t7 * t32;
t3 = t17 * t7 + t6 * t32;
t2 = t19 * t6 + t7 * t33;
t1 = t19 * t7 - t6 * t33;
t5 = [(-m(3) * t34 - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t47 * (t19 * pkin(6) + t34) - t44 * t19 + t46 * t17) * g(2) + (-t4 * mrSges(6,1) + t3 * mrSges(6,2) + (m(3) * pkin(1) + t47 * (-pkin(1) - pkin(6)) + t44) * t17 + (t48 * qJ(2) + t46) * t19) * g(1), -t40 * t48, (t16 * t42 + t51) * g(3) + t40 * ((mrSges(4,1) + t42) * t18 - t52 * t16), (-t16 * g(3) - t40 * t18) * t43, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(3) * (-mrSges(6,1) * t6 - mrSges(6,2) * t7) * t18];
taug = t5(:);
