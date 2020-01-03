% Calculate Gravitation load on the joints for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.18s
% Computational Cost: add. (59->25), mult. (112->28), div. (0->0), fcn. (81->4), ass. (0->13)
t11 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t14 = mrSges(4,2) + mrSges(5,2);
t5 = sin(qJ(3));
t7 = cos(qJ(3));
t25 = t11 * t5 + t14 * t7;
t6 = sin(qJ(1));
t8 = cos(qJ(1));
t22 = -g(1) * t6 + g(2) * t8;
t21 = m(3) + m(4) + m(5);
t20 = mrSges(2,2) - mrSges(3,3) - t25;
t19 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t4 = -qJ(4) - pkin(5);
t1 = [(-t21 * (t8 * pkin(1) + t6 * qJ(2)) + (-m(4) * pkin(5) + m(5) * t4 - t19) * t8 + t20 * t6) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(5)) - m(5) * (-pkin(1) + t4) + t19) * t6 + (-t21 * qJ(2) + t20) * t8) * g(1), t22 * t21, t25 * g(3) + t22 * (t11 * t7 - t14 * t5), (-g(1) * t8 - g(2) * t6) * m(5)];
taug = t1(:);
