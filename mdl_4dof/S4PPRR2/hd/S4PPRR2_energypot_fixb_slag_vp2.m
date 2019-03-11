% Calculate potential energy for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:48
% EndTime: 2019-03-08 18:16:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->39), mult. (43->28), div. (0->0), fcn. (18->6), ass. (0->12)
t27 = cos(pkin(6));
t29 = t27 * pkin(2) + pkin(1);
t28 = pkin(4) + qJ(2);
t25 = pkin(6) + qJ(3);
t26 = sin(pkin(6));
t23 = t26 * pkin(2);
t22 = qJ(4) + t25;
t21 = cos(t25);
t20 = sin(t25);
t19 = cos(t22);
t18 = sin(t22);
t1 = (-mrSges(1,3) + mrSges(2,2) - m(3) * qJ(2) - mrSges(3,3) - m(4) * t28 - mrSges(4,3) - m(5) * (pkin(5) + t28) - mrSges(5,3)) * g(3) + (-mrSges(1,2) - mrSges(2,3) - t26 * mrSges(3,1) - t27 * mrSges(3,2) - m(4) * t23 - t20 * mrSges(4,1) - t21 * mrSges(4,2) - m(5) * (pkin(3) * t20 + t23) - t18 * mrSges(5,1) - t19 * mrSges(5,2) + (-m(2) - m(3) - m(4) - m(5)) * qJ(1)) * g(2) + (-mrSges(1,1) - mrSges(2,1) - m(3) * pkin(1) - t27 * mrSges(3,1) + t26 * mrSges(3,2) - m(4) * t29 - t21 * mrSges(4,1) + t20 * mrSges(4,2) - m(5) * (pkin(3) * t21 + t29) - t19 * mrSges(5,1) + t18 * mrSges(5,2)) * g(1);
U  = t1;
