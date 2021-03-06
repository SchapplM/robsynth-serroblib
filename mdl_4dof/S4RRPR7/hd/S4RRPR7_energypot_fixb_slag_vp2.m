% Calculate potential energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:54
% EndTime: 2019-12-31 17:05:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (89->38), mult. (106->32), div. (0->0), fcn. (85->8), ass. (0->17)
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t61 = -m(5) * pkin(3) - t45 * mrSges(5,1) + t42 * mrSges(5,2) - mrSges(4,1);
t60 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t59 = m(4) + m(5);
t40 = qJ(2) + pkin(7);
t37 = sin(t40);
t38 = cos(t40);
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t57 = -m(3) * pkin(1) - t46 * mrSges(3,1) + t43 * mrSges(3,2) + t60 * t37 + t61 * t38 - mrSges(2,1);
t56 = m(3) * pkin(5) + t42 * mrSges(5,1) + t45 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t47 = cos(qJ(1));
t44 = sin(qJ(1));
t41 = -qJ(3) - pkin(5);
t36 = t46 * pkin(2) + pkin(1);
t1 = (-t43 * mrSges(3,1) - t46 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - t59 * (t43 * pkin(2) + pkin(4)) - t60 * t38 + t61 * t37 + (-m(2) - m(3)) * pkin(4)) * g(3) + (-mrSges(1,2) - t59 * (t44 * t36 + t47 * t41) + t56 * t47 + t57 * t44) * g(2) + (-mrSges(1,1) + (t59 * t41 - t56) * t44 + (-t59 * t36 + t57) * t47) * g(1);
U = t1;
