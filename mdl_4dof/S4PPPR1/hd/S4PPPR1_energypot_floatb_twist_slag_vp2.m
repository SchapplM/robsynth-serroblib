% Calculate potential energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:52
% EndTime: 2019-03-08 18:08:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (66->41), mult. (72->32), div. (0->0), fcn. (46->4), ass. (0->16)
t22 = -m(4) - m(5);
t10 = cos(pkin(5));
t9 = sin(pkin(5));
t20 = t10 * pkin(1) + t9 * qJ(2);
t19 = m(3) - t22;
t17 = qJ(1) + r_base(3);
t16 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t15 = pkin(2) + t17;
t14 = -m(1) - m(2) - t19;
t13 = m(5) * pkin(3) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t6 = t9 * pkin(1);
t2 = t10 * t11 + t12 * t9;
t1 = t10 * t12 - t9 * t11;
t3 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t15 + mrSges(4,2) - m(5) * (pkin(4) + t15) - mrSges(5,3) + (-m(2) - m(3)) * t17) * g(3) + (-mrSges(1,2) - m(3) * t6 + t1 * mrSges(5,1) - t2 * mrSges(5,2) + t14 * r_base(2) + t16 * t9 + (t19 * qJ(2) + t13) * t10 + t22 * (t9 * qJ(3) + t6)) * g(2) + (-m(3) * t20 - t2 * mrSges(5,1) - t1 * mrSges(5,2) + t16 * t10 - t13 * t9 + t14 * r_base(1) - mrSges(1,1) + t22 * (t10 * qJ(3) + t20)) * g(1);
U  = t3;
