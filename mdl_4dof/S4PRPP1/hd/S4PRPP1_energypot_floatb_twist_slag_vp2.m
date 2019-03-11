% Calculate potential energy for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:53
% EndTime: 2019-03-08 18:17:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (83->39), mult. (62->28), div. (0->0), fcn. (32->4), ass. (0->14)
t19 = m(4) + m(5);
t16 = qJ(1) + r_base(3);
t15 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t14 = pkin(4) + t16;
t13 = -m(1) - m(2) - m(3) - t19;
t12 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t11 = cos(pkin(5));
t10 = sin(pkin(5));
t9 = pkin(5) + qJ(2);
t8 = t11 * pkin(1);
t7 = t10 * pkin(1);
t5 = cos(t9);
t4 = sin(t9);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t16 - mrSges(2,3) - mrSges(3,3) - mrSges(4,1) - m(5) * (pkin(3) + t14) - mrSges(5,1) + (-m(3) - m(4)) * t14) * g(3) + (-mrSges(1,2) - t10 * mrSges(2,1) - mrSges(2,2) * t11 - m(3) * t7 + t13 * r_base(2) + (t19 * qJ(3) - t15) * t5 + t12 * t4 - t19 * (t4 * pkin(2) + t7)) * g(2) + (-m(3) * t8 - mrSges(2,1) * t11 + t10 * mrSges(2,2) + t12 * t5 + t13 * r_base(1) + t15 * t4 - mrSges(1,1) - t19 * (t5 * pkin(2) + t4 * qJ(3) + t8)) * g(1);
U  = t1;
