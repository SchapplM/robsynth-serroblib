% Calculate Gravitation load on the joints for
% S4PRPR1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (76->24), mult. (66->27), div. (0->0), fcn. (54->4), ass. (0->13)
t21 = m(4) + m(5);
t18 = cos(qJ(4));
t17 = sin(qJ(4));
t16 = mrSges(3,1) + mrSges(4,1);
t15 = mrSges(3,2) - mrSges(4,3);
t11 = pkin(6) + qJ(2);
t10 = cos(t11);
t9 = sin(t11);
t1 = -t10 * t18 - t9 * t17;
t2 = t10 * t17 - t9 * t18;
t13 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t12 = t1 * mrSges(5,1) + t2 * mrSges(5,2);
t3 = [(-m(2) - m(3) - t21) * g(3) (t15 * t9 + (-m(5) * pkin(3) - t16) * t10 + t12 - t21 * (t10 * pkin(2) + t9 * qJ(3))) * g(2) + ((m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(3)) + t16) * t9 + t13 + (-t21 * qJ(3) + t15) * t10) * g(1), t21 * (-g(1) * t9 + g(2) * t10) -g(1) * t13 - g(2) * t12];
taug  = t3(:);
