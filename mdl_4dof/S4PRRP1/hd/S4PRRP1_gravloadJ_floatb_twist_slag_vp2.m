% Calculate Gravitation load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.09s
% Computational Cost: add. (94->19), mult. (56->19), div. (0->0), fcn. (34->4), ass. (0->13)
t19 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t18 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t16 = m(4) + m(5);
t17 = t16 * pkin(2) + mrSges(3,1);
t13 = pkin(6) + qJ(2);
t12 = qJ(3) + t13;
t8 = sin(t12);
t9 = cos(t12);
t15 = t18 * t8 - t19 * t9;
t14 = t18 * t9 + t19 * t8;
t11 = cos(t13);
t10 = sin(t13);
t1 = [(-m(2) - m(3) - t16) * g(3) (t10 * mrSges(3,2) - t17 * t11 + t15) * g(2) + (t11 * mrSges(3,2) + t17 * t10 + t14) * g(1), t14 * g(1) + t15 * g(2) (-g(1) * t8 + g(2) * t9) * m(5)];
taug  = t1(:);
