% Calculate potential energy for
% S4PRRR5
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:32
% EndTime: 2019-12-31 16:33:32
% DurationCPUTime: 0.20s
% Computational Cost: add. (89->38), mult. (106->32), div. (0->0), fcn. (85->8), ass. (0->17)
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t61 = -m(5) * pkin(3) - t45 * mrSges(5,1) + t43 * mrSges(5,2) - mrSges(4,1);
t60 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t59 = m(4) + m(5);
t40 = qJ(2) + qJ(3);
t37 = sin(t40);
t38 = cos(t40);
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t57 = -m(3) * pkin(1) - t46 * mrSges(3,1) + t44 * mrSges(3,2) + t60 * t37 + t61 * t38 - mrSges(2,1);
t56 = m(3) * pkin(4) + t43 * mrSges(5,1) + t45 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t47 = -pkin(5) - pkin(4);
t42 = cos(pkin(7));
t41 = sin(pkin(7));
t36 = t46 * pkin(2) + pkin(1);
t1 = (-t44 * mrSges(3,1) - t46 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - t59 * (t44 * pkin(2) + qJ(1)) - t60 * t38 + t61 * t37 + (-m(2) - m(3)) * qJ(1)) * g(3) + (-mrSges(1,2) - t59 * (t41 * t36 + t42 * t47) + t56 * t42 + t57 * t41) * g(2) + (-mrSges(1,1) + (t59 * t47 - t56) * t41 + (-t59 * t36 + t57) * t42) * g(1);
U = t1;
