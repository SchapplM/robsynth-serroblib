% Calculate Gravitation load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:37
% DurationCPUTime: 0.20s
% Computational Cost: add. (190->40), mult. (174->48), div. (0->0), fcn. (172->6), ass. (0->22)
t11 = sin(qJ(5));
t12 = cos(qJ(5));
t17 = -mrSges(6,1) * t12 + t11 * mrSges(6,2);
t32 = m(6) * pkin(4) + mrSges(5,1) - t17;
t31 = mrSges(3,1) + mrSges(4,1);
t30 = mrSges(3,2) - mrSges(4,3);
t24 = m(4) + m(5) + m(6);
t22 = pkin(8) + qJ(2);
t10 = cos(t22);
t20 = sin(t22);
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t1 = -t10 * t26 - t20 * t25;
t2 = t10 * t25 - t20 * t26;
t29 = t2 * mrSges(5,2) + t32 * t1;
t28 = t1 * mrSges(5,2) - t32 * t2;
t27 = t10 * pkin(2) + t20 * qJ(3);
t23 = t10 * pkin(3) + t27;
t21 = -m(6) * pkin(7) - mrSges(6,3);
t15 = -t20 * pkin(2) + t10 * qJ(3);
t13 = -t20 * pkin(3) + t15;
t3 = [(-m(2) - m(3) - t24) * g(3), (-m(4) * t27 - m(5) * t23 - m(6) * (pkin(7) * t2 + t23) - t2 * mrSges(6,3) + t30 * t20 - t31 * t10 + t29) * g(2) + (-m(4) * t15 - m(5) * t13 - m(6) * (t1 * pkin(7) + t13) - t1 * mrSges(6,3) + t31 * t20 + t30 * t10 + t28) * g(1), (-g(1) * t20 + g(2) * t10) * t24, (-t21 * t2 - t29) * g(2) + (-t21 * t1 - t28) * g(1), -g(3) * t17 + (-g(1) * t1 - g(2) * t2) * (mrSges(6,1) * t11 + mrSges(6,2) * t12)];
taug = t3(:);
