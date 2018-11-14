% Calculate potential energy for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PRPP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:08
% EndTime: 2018-11-14 14:01:08
% DurationCPUTime: 0.06s
% Computational Cost: add. (41->28), mult. (47->17), div. (0->0), fcn. (22->2), ass. (0->7)
t25 = m(4) + m(5);
t24 = -m(3) - t25;
t22 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t21 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1) - mrSges(5,1);
t20 = cos(qJ(2));
t19 = sin(qJ(2));
t1 = (m(5) * qJ(4) + t24 * pkin(4) + mrSges(2,2) - mrSges(4,2) - mrSges(1,3) - mrSges(3,3) + mrSges(5,3)) * g(3) + (-mrSges(1,2) - mrSges(2,3) + (-m(2) + t24) * qJ(1) + (t25 * qJ(3) - t22) * t20 + (-t25 * pkin(2) + t21) * t19) * g(2) + (-m(3) * pkin(1) + t22 * t19 + t21 * t20 - mrSges(1,1) - mrSges(2,1) - t25 * (t20 * pkin(2) + t19 * qJ(3) + pkin(1))) * g(1);
U  = t1;
