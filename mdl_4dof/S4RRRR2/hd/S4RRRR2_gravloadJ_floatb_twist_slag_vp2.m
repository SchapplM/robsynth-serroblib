% Calculate Gravitation load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (165->41), mult. (149->44), div. (0->0), fcn. (113->8), ass. (0->25)
t21 = qJ(3) + qJ(4);
t16 = sin(t21);
t18 = cos(t21);
t46 = t18 * mrSges(5,1) - t16 * mrSges(5,2);
t23 = sin(qJ(3));
t45 = t23 * mrSges(4,2) - t46;
t25 = cos(qJ(3));
t44 = -t25 * mrSges(4,1) - mrSges(3,1) + t45;
t42 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t22 = qJ(1) + qJ(2);
t17 = sin(t22);
t19 = cos(t22);
t41 = g(1) * t19 + g(2) * t17;
t35 = t19 * pkin(2) + t17 * pkin(6);
t34 = m(5) * pkin(3) + mrSges(4,1);
t15 = pkin(3) * t25 + pkin(2);
t27 = -pkin(7) - pkin(6);
t33 = t19 * t15 - t17 * t27;
t32 = mrSges(5,1) * t16 + mrSges(5,2) * t18;
t29 = t42 * t17 + t44 * t19;
t28 = (m(4) * pkin(2) + m(5) * t15 - t44) * t17 + (-m(4) * pkin(6) + m(5) * t27 + t42) * t19;
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t20 = t26 * pkin(1);
t1 = [(mrSges(2,2) * t24 - m(4) * (t20 + t35) - m(5) * (t20 + t33) + (-m(3) * pkin(1) - mrSges(2,1)) * t26 + t29) * g(2) + (mrSges(2,2) * t26 + (mrSges(2,1) + (m(3) + m(4) + m(5)) * pkin(1)) * t24 + t28) * g(1), (-m(4) * t35 - m(5) * t33 + t29) * g(2) + t28 * g(1), (-t34 * t25 + t45) * g(3) + t41 * (mrSges(4,2) * t25 + t34 * t23 + t32), -g(3) * t46 + t41 * t32];
taug = t1(:);
