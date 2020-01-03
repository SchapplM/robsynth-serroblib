% Calculate Gravitation load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (52->28), mult. (124->40), div. (0->0), fcn. (100->6), ass. (0->14)
t33 = mrSges(3,1) + mrSges(5,3) - mrSges(4,2);
t7 = sin(qJ(4));
t9 = cos(qJ(4));
t32 = -t7 * mrSges(5,1) - t9 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t5 = sin(pkin(6));
t6 = cos(pkin(6));
t25 = g(1) * t6 + g(2) * t5;
t19 = m(4) + m(5);
t8 = sin(qJ(2));
t21 = t7 * t8;
t20 = t8 * t9;
t10 = cos(qJ(2));
t17 = g(3) * t10;
t1 = [(-m(2) - m(3) - t19) * g(3), (-t19 * (t10 * pkin(2) + t8 * qJ(3)) + t32 * t8 + (-m(5) * pkin(5) - t33) * t10) * g(3) + ((-m(5) * (-pkin(2) - pkin(5)) + m(4) * pkin(2) + t33) * t8 + (-t19 * qJ(3) + t32) * t10) * t25, (-t25 * t8 + t17) * t19, -g(1) * ((t6 * t20 - t5 * t7) * mrSges(5,1) + (-t6 * t21 - t5 * t9) * mrSges(5,2)) - g(2) * ((t5 * t20 + t6 * t7) * mrSges(5,1) + (-t5 * t21 + t6 * t9) * mrSges(5,2)) - (-mrSges(5,1) * t9 + mrSges(5,2) * t7) * t17];
taug = t1(:);
