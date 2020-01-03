% Calculate potential energy for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:28
% EndTime: 2019-12-31 16:32:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (84->36), mult. (74->28), div. (0->0), fcn. (49->8), ass. (0->15)
t45 = -m(3) - m(4) - m(5);
t46 = pkin(1) * t45 - mrSges(2,1);
t33 = qJ(3) + qJ(4);
t28 = sin(t33);
t29 = cos(t33);
t37 = sin(qJ(3));
t38 = cos(qJ(3));
t44 = -mrSges(3,1) - m(4) * pkin(2) - mrSges(4,1) * t38 + mrSges(4,2) * t37 - m(5) * (pkin(3) * t38 + pkin(2)) - mrSges(5,1) * t29 + mrSges(5,2) * t28;
t42 = m(4) * pkin(5) - m(5) * (-pkin(6) - pkin(5)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t35 = cos(pkin(7));
t34 = sin(pkin(7));
t32 = pkin(7) + qJ(2);
t27 = cos(t32);
t26 = sin(t32);
t1 = (-m(2) * qJ(1) - t28 * mrSges(5,1) - t38 * mrSges(4,2) - t29 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t37 + t45 * (pkin(4) + qJ(1))) * g(3) + (-t35 * mrSges(2,2) + t44 * t26 + t42 * t27 + t46 * t34 - mrSges(1,2)) * g(2) + (t34 * mrSges(2,2) - t42 * t26 + t44 * t27 + t46 * t35 - mrSges(1,1)) * g(1);
U = t1;
