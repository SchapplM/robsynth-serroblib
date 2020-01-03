% Calculate Gravitation load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:45
% DurationCPUTime: 0.22s
% Computational Cost: add. (130->35), mult. (132->38), div. (0->0), fcn. (93->6), ass. (0->18)
t11 = cos(qJ(4));
t15 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t20 = mrSges(5,2) + mrSges(6,2);
t9 = sin(qJ(4));
t30 = t20 * t11 + t15 * t9;
t7 = qJ(1) + pkin(7);
t4 = sin(t7);
t5 = cos(t7);
t27 = -g(1) * t4 + g(2) * t5;
t19 = m(4) + m(5) + m(6);
t26 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t25 = mrSges(3,2) - mrSges(4,3) - t30;
t10 = sin(qJ(1));
t21 = pkin(1) * t10;
t12 = cos(qJ(1));
t6 = t12 * pkin(1);
t8 = -qJ(5) - pkin(6);
t1 = [(-m(3) * t6 - mrSges(2,1) * t12 + t10 * mrSges(2,2) - t19 * (t5 * pkin(2) + t4 * qJ(3) + t6) + (-m(5) * pkin(6) + m(6) * t8 - t26) * t5 + t25 * t4) * g(2) + (m(3) * t21 + t10 * mrSges(2,1) + mrSges(2,2) * t12 - t19 * (t5 * qJ(3) - t21) + t25 * t5 + (m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(6)) - m(6) * (-pkin(2) + t8) + t26) * t4) * g(1), (-m(3) - t19) * g(3), t27 * t19, t30 * g(3) + t27 * (t15 * t11 - t20 * t9), (-g(1) * t5 - g(2) * t4) * m(6)];
taug = t1(:);
