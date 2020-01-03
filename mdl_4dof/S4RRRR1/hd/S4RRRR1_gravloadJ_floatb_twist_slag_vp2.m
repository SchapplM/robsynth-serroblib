% Calculate Gravitation load on the joints for
% S4RRRR1
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:08
% DurationCPUTime: 0.15s
% Computational Cost: add. (178->37), mult. (123->40), div. (0->0), fcn. (90->8), ass. (0->23)
t34 = mrSges(4,2) - mrSges(5,3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t33 = t20 * mrSges(5,1) - t18 * mrSges(5,2);
t32 = -mrSges(4,1) - t33;
t31 = m(4) + m(5);
t17 = qJ(1) + qJ(2);
t15 = qJ(3) + t17;
t11 = sin(t15);
t12 = cos(t15);
t30 = t12 * pkin(3) + t11 * pkin(7);
t14 = cos(t17);
t10 = pkin(2) * t14;
t27 = t10 + t30;
t25 = t34 * t11 + t32 * t12;
t13 = sin(t17);
t24 = -t14 * mrSges(3,1) + t13 * mrSges(3,2) + t25;
t23 = (m(5) * pkin(3) - t32) * t11 + (-m(5) * pkin(7) + t34) * t12;
t22 = mrSges(3,2) * t14 + (t31 * pkin(2) + mrSges(3,1)) * t13 + t23;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t16 = t21 * pkin(1);
t1 = [(t19 * mrSges(2,2) - m(4) * (t10 + t16) - m(5) * (t16 + t27) + (-m(3) * pkin(1) - mrSges(2,1)) * t21 + t24) * g(2) + (mrSges(2,2) * t21 + (mrSges(2,1) + (m(3) + t31) * pkin(1)) * t19 + t22) * g(1), (-m(4) * t10 - m(5) * t27 + t24) * g(2) + t22 * g(1), (-m(5) * t30 + t25) * g(2) + t23 * g(1), -g(3) * t33 + (g(1) * t12 + g(2) * t11) * (mrSges(5,1) * t18 + mrSges(5,2) * t20)];
taug = t1(:);
