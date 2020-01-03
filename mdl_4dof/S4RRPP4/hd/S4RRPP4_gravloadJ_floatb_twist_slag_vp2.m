% Calculate Gravitation load on the joints for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:58:54
% DurationCPUTime: 0.32s
% Computational Cost: add. (83->40), mult. (172->41), div. (0->0), fcn. (135->4), ass. (0->19)
t43 = -mrSges(3,1) - mrSges(4,1);
t42 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t35 = m(4) + m(5);
t41 = -m(3) - t35;
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t33 = g(1) * t13 + g(2) * t11;
t10 = sin(qJ(2));
t5 = t10 * qJ(3);
t12 = cos(qJ(2));
t7 = t12 * pkin(2);
t28 = t7 + t5;
t39 = g(1) * t11;
t32 = t42 * t10 + t43 * t12;
t30 = -mrSges(2,1) + t32;
t29 = m(5) * qJ(4) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t21 = -pkin(1) - t5;
t20 = m(5) * (-pkin(2) - pkin(3)) - mrSges(5,1);
t1 = [(t29 * t11 + t41 * (t13 * pkin(1) + t11 * pkin(5))) * g(2) + (m(3) * pkin(1) - m(4) * (t21 - t7) - m(5) * t21 - t20 * t12 - t30) * t39 + ((-t35 * t28 - (m(5) * pkin(3) + mrSges(5,1)) * t12 + t30) * g(2) + (t41 * pkin(5) + t29) * g(1)) * t13, (-m(4) * t28 - m(5) * (pkin(3) * t12 + t28) - t12 * mrSges(5,1) + t32) * g(3) + ((m(4) * pkin(2) - t20 - t43) * t10 + (-qJ(3) * t35 + t42) * t12) * t33, (g(3) * t12 - t10 * t33) * t35, (-g(2) * t13 + t39) * m(5)];
taug = t1(:);
