% Calculate potential energy for
% S4PPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR2_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:55:24
% EndTime: 2018-11-14 13:55:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->33), mult. (55->25), div. (0->0), fcn. (34->4), ass. (0->11)
t33 = m(4) + m(5);
t32 = mrSges(3,2) - mrSges(4,3);
t31 = -m(3) - t33;
t29 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t28 = cos(qJ(4));
t27 = sin(qJ(4));
t26 = cos(pkin(5));
t25 = sin(pkin(5));
t21 = t25 * t28 - t26 * t27;
t20 = -t25 * t27 - t26 * t28;
t1 = (m(5) * pkin(4) + t31 * qJ(2) + mrSges(2,2) - mrSges(4,2) - mrSges(1,3) - mrSges(3,3) + mrSges(5,3)) * g(3) + (-t21 * mrSges(5,1) - t20 * mrSges(5,2) - mrSges(1,2) - mrSges(2,3) + (t33 * qJ(3) - t32) * t26 + (-m(2) + t31) * qJ(1) + (-t33 * pkin(2) + t29) * t25) * g(2) + (-m(3) * pkin(1) + t20 * mrSges(5,1) - t21 * mrSges(5,2) + t32 * t25 + t29 * t26 - mrSges(1,1) - mrSges(2,1) - t33 * (t26 * pkin(2) + t25 * qJ(3) + pkin(1))) * g(1);
U  = t1;
