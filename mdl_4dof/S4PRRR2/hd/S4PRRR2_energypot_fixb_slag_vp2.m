% Calculate potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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

function U = S4PRRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->28), mult. (41->18), div. (0->0), fcn. (18->6), ass. (0->12)
t23 = m(4) + m(5);
t24 = t23 * pkin(1) + mrSges(3,1);
t19 = qJ(2) + qJ(3);
t22 = m(5) * pkin(2) + mrSges(4,1);
t21 = cos(qJ(2));
t20 = sin(qJ(2));
t18 = qJ(4) + t19;
t17 = cos(t19);
t16 = sin(t19);
t15 = cos(t18);
t14 = sin(t18);
t1 = (-t15 * mrSges(5,1) + t20 * mrSges(3,2) + t16 * mrSges(4,2) + t14 * mrSges(5,2) - t22 * t17 - t24 * t21 - mrSges(2,1) - mrSges(1,3)) * g(3) + (-mrSges(1,2) + mrSges(2,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + (m(2) + m(3) + t23) * qJ(1)) * g(2) + (t14 * mrSges(5,1) + t21 * mrSges(3,2) + t17 * mrSges(4,2) + t15 * mrSges(5,2) + t22 * t16 + t24 * t20 - mrSges(1,1) + mrSges(2,2)) * g(1);
U  = t1;
