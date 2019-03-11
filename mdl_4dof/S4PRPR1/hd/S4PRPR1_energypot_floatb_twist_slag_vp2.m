% Calculate potential energy for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:51
% EndTime: 2019-03-08 18:20:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (91->44), mult. (70->36), div. (0->0), fcn. (44->6), ass. (0->18)
t23 = m(4) + m(5);
t21 = mrSges(3,2) - mrSges(4,3);
t19 = qJ(1) + r_base(3);
t18 = pkin(4) + t19;
t17 = -m(1) - m(2) - m(3) - t23;
t16 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t13 = cos(pkin(6));
t12 = sin(pkin(6));
t11 = pkin(6) + qJ(2);
t10 = t13 * pkin(1);
t9 = t12 * pkin(1);
t7 = cos(t11);
t6 = sin(t11);
t2 = -t7 * t14 + t15 * t6;
t1 = -t6 * t14 - t15 * t7;
t3 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t19 - mrSges(2,3) - mrSges(3,3) - mrSges(4,2) - m(5) * (-pkin(5) + t18) + mrSges(5,3) + (-m(3) - m(4)) * t18) * g(3) + (-mrSges(1,2) - mrSges(2,1) * t12 - t13 * mrSges(2,2) - m(3) * t9 - t2 * mrSges(5,1) - t1 * mrSges(5,2) + t17 * r_base(2) + (t23 * qJ(3) - t21) * t7 + t16 * t6 - t23 * (t6 * pkin(2) + t9)) * g(2) + (-m(3) * t10 - mrSges(2,1) * t13 + t1 * mrSges(5,1) + mrSges(2,2) * t12 - t2 * mrSges(5,2) + t16 * t7 + t17 * r_base(1) + t21 * t6 - mrSges(1,1) - t23 * (t7 * pkin(2) + t6 * qJ(3) + t10)) * g(1);
U  = t3;
