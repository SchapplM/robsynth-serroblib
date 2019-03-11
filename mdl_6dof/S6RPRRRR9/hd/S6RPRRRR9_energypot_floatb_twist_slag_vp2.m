% Calculate potential energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:22:56
% EndTime: 2019-03-09 07:22:57
% DurationCPUTime: 0.47s
% Computational Cost: add. (195->66), mult. (214->50), div. (0->0), fcn. (184->10), ass. (0->25)
t18 = sin(qJ(4));
t21 = cos(qJ(4));
t17 = qJ(4) + qJ(5);
t9 = qJ(6) + t17;
t3 = sin(t9);
t4 = cos(t9);
t5 = t21 * pkin(4) + pkin(3);
t6 = sin(t17);
t7 = cos(t17);
t46 = mrSges(4,1) + m(5) * pkin(3) + t21 * mrSges(5,1) - t18 * mrSges(5,2) + m(7) * (pkin(5) * t7 + t5) + t4 * mrSges(7,1) - t3 * mrSges(7,2) + m(6) * t5 + t7 * mrSges(6,1) - t6 * mrSges(6,2);
t24 = -pkin(9) - pkin(8);
t45 = mrSges(4,2) - m(5) * pkin(8) + m(6) * t24 + m(7) * (-pkin(10) + t24) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t44 = t46 * t19 + t45 * t22 - mrSges(2,2) + mrSges(3,3);
t43 = -m(1) - m(2);
t40 = -m(4) - m(5) - m(6) - m(7);
t37 = pkin(4) * t18;
t39 = -m(6) * t37 - t6 * mrSges(6,1) - t7 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t6 + t37) - t3 * mrSges(7,1) - t4 * mrSges(7,2) - t18 * mrSges(5,1) - t21 * mrSges(5,2);
t15 = pkin(6) + r_base(3);
t20 = sin(qJ(1));
t36 = t20 * pkin(1) + r_base(2);
t23 = cos(qJ(1));
t34 = t23 * pkin(1) + t20 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t15 + t40 * (pkin(2) + t15) - t46 * t22 + t45 * t19) * g(3) + (-m(3) * t36 - mrSges(1,2) + t43 * r_base(2) + t40 * (t20 * pkin(7) + t36) + ((m(3) - t40) * qJ(2) + t44) * t23 + t39 * t20) * g(2) + (-m(3) * t34 - mrSges(1,1) + t43 * r_base(1) + t40 * (t23 * pkin(7) + t34) + t39 * t23 - t44 * t20) * g(1);
U  = t1;
