% Calculate Gravitation load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (129->38), mult. (115->40), div. (0->0), fcn. (80->8), ass. (0->18)
t30 = m(5) + m(6);
t10 = qJ(1) + pkin(7);
t5 = sin(t10);
t7 = cos(t10);
t29 = -g(1) * t5 + g(2) * t7;
t23 = m(4) + t30;
t28 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t9 = pkin(8) + qJ(5);
t4 = sin(t9);
t6 = cos(t9);
t18 = -t4 * mrSges(6,1) - t6 * mrSges(6,2);
t27 = -cos(pkin(8)) * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + t18 + (-m(6) * pkin(4) - mrSges(5,1)) * sin(pkin(8));
t14 = sin(qJ(1));
t24 = pkin(1) * t14;
t15 = cos(qJ(1));
t8 = t15 * pkin(1);
t13 = -pkin(6) - qJ(4);
t1 = [(-m(3) * t8 - mrSges(2,1) * t15 + t14 * mrSges(2,2) - t23 * (t7 * pkin(2) + t5 * qJ(3) + t8) + (-m(5) * qJ(4) + m(6) * t13 - t28) * t7 + t27 * t5) * g(2) + (m(3) * t24 + t14 * mrSges(2,1) + mrSges(2,2) * t15 - t23 * (t7 * qJ(3) - t24) + t27 * t7 + (m(4) * pkin(2) - m(5) * (-pkin(2) - qJ(4)) - m(6) * (-pkin(2) + t13) + t28) * t5) * g(1), (-m(3) - t23) * g(3), t29 * t23, t30 * (-g(1) * t7 - g(2) * t5), -g(3) * t18 + t29 * (mrSges(6,1) * t6 - mrSges(6,2) * t4)];
taug = t1(:);
