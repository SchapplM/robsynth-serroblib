% Calculate potential energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:00
% EndTime: 2019-12-31 16:47:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (56->32), mult. (83->25), div. (0->0), fcn. (58->4), ass. (0->12)
t45 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t44 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t43 = -m(4) - m(5);
t42 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t41 = t45 * t30 + t44 * t32 - mrSges(2,2) + mrSges(3,3);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t38 = t33 * pkin(1) + t31 * qJ(2);
t27 = t31 * pkin(1);
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t43 * (pkin(2) + pkin(4)) - t45 * t32 + t44 * t30 + (-m(2) - m(3)) * pkin(4)) * g(3) + (-m(3) * t27 - mrSges(1,2) + t43 * (t31 * pkin(5) + t27) + ((m(3) - t43) * qJ(2) + t41) * t33 + t42 * t31) * g(2) + (-m(3) * t38 - mrSges(1,1) + t43 * (t33 * pkin(5) + t38) + t42 * t33 - t41 * t31) * g(1);
U = t1;
