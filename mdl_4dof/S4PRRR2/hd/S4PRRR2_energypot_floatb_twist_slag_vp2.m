% Calculate potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->35), mult. (48->25), div. (0->0), fcn. (18->6), ass. (0->13)
t13 = m(4) + m(5);
t8 = qJ(2) + qJ(3);
t12 = -m(2) - m(3) - t13;
t11 = -m(1) + t12;
t10 = cos(qJ(2));
t9 = sin(qJ(2));
t6 = t10 * pkin(1);
t5 = qJ(4) + t8;
t4 = cos(t8);
t3 = sin(t8);
t2 = cos(t5);
t1 = sin(t5);
t7 = (-mrSges(1,3) - mrSges(2,1) - mrSges(3,1) * t10 + t9 * mrSges(3,2) - m(4) * t6 - t4 * mrSges(4,1) + t3 * mrSges(4,2) - m(5) * (pkin(2) * t4 + t6) - t2 * mrSges(5,1) + t1 * mrSges(5,2) + t11 * r_base(3)) * g(3) + (-m(1) * r_base(2) - mrSges(1,2) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + t12 * (-qJ(1) + r_base(2))) * g(2) + (t1 * mrSges(5,1) + mrSges(3,2) * t10 + t4 * mrSges(4,2) + t2 * mrSges(5,2) - mrSges(1,1) + mrSges(2,2) + (m(5) * pkin(2) + mrSges(4,1)) * t3 + (t13 * pkin(1) + mrSges(3,1)) * t9 + t11 * r_base(1)) * g(1);
U  = t7;
