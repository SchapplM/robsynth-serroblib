% Calculate potential energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:27
% EndTime: 2018-11-14 13:56:27
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->32), mult. (35->20), div. (0->0), fcn. (10->4), ass. (0->8)
t19 = -qJ(3) - pkin(1);
t18 = -m(3) - m(4) - m(5);
t17 = cos(pkin(5));
t16 = sin(pkin(5));
t15 = pkin(5) + qJ(4);
t14 = cos(t15);
t13 = sin(t15);
t1 = (-mrSges(1,3) + mrSges(2,1) + m(3) * pkin(1) - mrSges(3,2) - m(4) * t19 + mrSges(4,3) - m(5) * (-pkin(4) + t19) + mrSges(5,3)) * g(3) + (-mrSges(1,2) - mrSges(2,3) - mrSges(3,1) - m(4) * pkin(2) - t17 * mrSges(4,1) + t16 * mrSges(4,2) - m(5) * (t17 * pkin(3) + pkin(2)) - t14 * mrSges(5,1) + t13 * mrSges(5,2) + (-m(2) + t18) * qJ(1)) * g(2) + (-t13 * mrSges(5,1) - t17 * mrSges(4,2) - t14 * mrSges(5,2) - mrSges(1,1) + mrSges(2,2) - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t16 + t18 * qJ(2)) * g(1);
U  = t1;
