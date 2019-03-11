% Calculate potential energy for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:38
% EndTime: 2019-03-08 18:31:38
% DurationCPUTime: 0.13s
% Computational Cost: add. (91->47), mult. (58->38), div. (0->0), fcn. (28->8), ass. (0->20)
t14 = sin(qJ(1));
t11 = t14 * pkin(1);
t13 = qJ(1) + pkin(7);
t8 = sin(t13);
t21 = pkin(2) * t8 + t11;
t15 = cos(qJ(1));
t12 = t15 * pkin(1);
t9 = cos(t13);
t20 = pkin(2) * t9 + t12;
t19 = pkin(4) + r_base(3);
t10 = qJ(3) + t13;
t18 = qJ(2) + t19;
t17 = -m(1) - m(2) - m(3) - m(4) - m(5);
t16 = pkin(5) + t18;
t7 = qJ(4) + t10;
t6 = cos(t10);
t5 = sin(t10);
t2 = cos(t7);
t1 = sin(t7);
t3 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t19 - mrSges(2,3) - m(3) * t18 - mrSges(3,3) - m(4) * t16 - mrSges(4,3) - m(5) * (pkin(6) + t16) - mrSges(5,3)) * g(3) + (-mrSges(1,2) - t14 * mrSges(2,1) - t15 * mrSges(2,2) - m(3) * t11 - t8 * mrSges(3,1) - t9 * mrSges(3,2) - m(4) * t21 - t5 * mrSges(4,1) - t6 * mrSges(4,2) - m(5) * (pkin(3) * t5 + t21) - t1 * mrSges(5,1) - t2 * mrSges(5,2) + t17 * r_base(2)) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t15 + t14 * mrSges(2,2) - m(3) * t12 - t9 * mrSges(3,1) + t8 * mrSges(3,2) - m(4) * t20 - t6 * mrSges(4,1) + t5 * mrSges(4,2) - m(5) * (pkin(3) * t6 + t20) - t2 * mrSges(5,1) + t1 * mrSges(5,2) + t17 * r_base(1)) * g(1);
U  = t3;
