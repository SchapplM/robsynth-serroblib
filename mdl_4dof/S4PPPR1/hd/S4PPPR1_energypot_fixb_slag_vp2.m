% Calculate potential energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:03
% EndTime: 2018-11-14 13:38:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (51->35), mult. (67->29), div. (0->0), fcn. (46->4), ass. (0->13)
t40 = -m(4) - m(5);
t30 = sin(pkin(5));
t31 = cos(pkin(5));
t38 = t31 * pkin(1) + t30 * qJ(2);
t37 = m(3) - t40;
t36 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t34 = m(5) * pkin(3) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t33 = cos(qJ(4));
t32 = sin(qJ(4));
t28 = t30 * pkin(1);
t24 = t30 * t33 + t31 * t32;
t23 = -t30 * t32 + t31 * t33;
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * pkin(2) + mrSges(4,2) - m(5) * (pkin(4) + pkin(2)) - mrSges(5,3) + (-m(2) - t37) * qJ(1)) * g(3) + (-mrSges(1,2) - m(3) * t28 + t23 * mrSges(5,1) - t24 * mrSges(5,2) + t36 * t30 + (t37 * qJ(2) + t34) * t31 + t40 * (t30 * qJ(3) + t28)) * g(2) + (-m(3) * t38 - t24 * mrSges(5,1) - t23 * mrSges(5,2) - t34 * t30 + t36 * t31 - mrSges(1,1) + t40 * (t31 * qJ(3) + t38)) * g(1);
U  = t1;
