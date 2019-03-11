% Calculate Gravitation load on the joints for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:14
% EndTime: 2019-03-08 18:25:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (100->21), mult. (54->21), div. (0->0), fcn. (32->6), ass. (0->17)
t19 = m(5) * pkin(3) + mrSges(4,1);
t16 = m(4) + m(5);
t18 = t16 * pkin(2) + mrSges(3,1);
t11 = pkin(7) + qJ(2);
t10 = qJ(3) + t11;
t7 = qJ(4) + t10;
t3 = sin(t7);
t4 = cos(t7);
t15 = t4 * mrSges(5,1) - t3 * mrSges(5,2);
t14 = -t3 * mrSges(5,1) - t4 * mrSges(5,2);
t5 = sin(t10);
t6 = cos(t10);
t13 = t5 * mrSges(4,2) - t19 * t6 - t15;
t12 = t6 * mrSges(4,2) + t19 * t5 - t14;
t9 = cos(t11);
t8 = sin(t11);
t1 = [(-m(2) - m(3) - t16) * g(3) (t8 * mrSges(3,2) - t18 * t9 + t13) * g(2) + (t9 * mrSges(3,2) + t18 * t8 + t12) * g(1), t12 * g(1) + t13 * g(2), -g(1) * t14 - g(2) * t15];
taug  = t1(:);
