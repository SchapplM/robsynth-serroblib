% Calculate potential energy for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:48
% EndTime: 2019-03-08 18:18:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->32), mult. (45->20), div. (0->0), fcn. (20->4), ass. (0->9)
t29 = -m(4) - m(5);
t27 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t26 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t25 = cos(qJ(2));
t24 = sin(qJ(2));
t22 = qJ(2) + pkin(5);
t19 = cos(t22);
t18 = sin(t22);
t1 = (-m(3) * pkin(4) + mrSges(2,2) - mrSges(5,2) - mrSges(1,3) - mrSges(3,3) - mrSges(4,3) + t29 * (qJ(3) + pkin(4))) * g(3) + (-t25 * mrSges(3,2) - mrSges(1,2) - mrSges(2,3) + t26 * t19 + t27 * t18 + (-m(2) - m(3) + t29) * qJ(1) + (t29 * pkin(2) - mrSges(3,1)) * t24) * g(2) + (-m(3) * pkin(1) - t25 * mrSges(3,1) + t24 * mrSges(3,2) - t26 * t18 + t27 * t19 - mrSges(1,1) - mrSges(2,1) + t29 * (t25 * pkin(2) + pkin(1))) * g(1);
U  = t1;
