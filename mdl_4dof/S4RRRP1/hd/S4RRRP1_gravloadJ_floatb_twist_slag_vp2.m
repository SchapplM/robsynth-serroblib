% Calculate Gravitation load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:30
% EndTime: 2018-11-14 13:54:30
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->32), mult. (87->31), div. (0->0), fcn. (56->6), ass. (0->19)
t22 = -mrSges(4,1) - mrSges(5,1);
t12 = qJ(1) + qJ(2);
t9 = cos(t12);
t5 = pkin(2) * t9;
t10 = qJ(3) + t12;
t7 = cos(t10);
t3 = pkin(3) * t7;
t21 = m(4) + m(5);
t14 = cos(qJ(1));
t20 = t14 * pkin(1) + t5;
t19 = mrSges(4,2) + mrSges(5,2);
t6 = sin(t10);
t18 = t19 * t6 + t22 * t7;
t8 = sin(t12);
t17 = -t9 * mrSges(3,1) + t8 * mrSges(3,2) + t18;
t16 = t19 * t7 + (m(5) * pkin(3) - t22) * t6;
t15 = mrSges(3,2) * t9 + (t21 * pkin(2) + mrSges(3,1)) * t8 + t16;
t13 = sin(qJ(1));
t1 = [(t13 * mrSges(2,2) - m(4) * t20 - m(5) * (t3 + t20) + (-m(3) * pkin(1) - mrSges(2,1)) * t14 + t17) * g(2) + (mrSges(2,2) * t14 + (mrSges(2,1) + (m(3) + t21) * pkin(1)) * t13 + t15) * g(1) (-m(4) * t5 - m(5) * (t3 + t5) + t17) * g(2) + t15 * g(1) (-m(5) * t3 + t18) * g(2) + t16 * g(1), -g(3) * m(5)];
taug  = t1(:);
