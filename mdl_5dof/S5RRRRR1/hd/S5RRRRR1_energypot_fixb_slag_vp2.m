% Calculate potential energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S5RRRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:50:56
% EndTime: 2018-11-16 14:50:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (132->43), mult. (129->35), div. (0->0), fcn. (102->10), ass. (0->20)
t64 = -m(5) - m(6);
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t63 = m(6) * pkin(4) + t49 * mrSges(6,1) - t46 * mrSges(6,2) + mrSges(5,1);
t62 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t45 = qJ(2) + qJ(3);
t43 = qJ(4) + t45;
t38 = sin(t43);
t39 = cos(t43);
t50 = cos(qJ(2));
t40 = t50 * pkin(2) + pkin(1);
t41 = sin(t45);
t42 = cos(t45);
t47 = sin(qJ(2));
t61 = -m(3) * pkin(1) - m(4) * t40 - mrSges(3,1) * t50 - mrSges(4,1) * t42 + mrSges(3,2) * t47 + mrSges(4,2) * t41 - mrSges(2,1) - t63 * t39 + t64 * (pkin(3) * t42 + t40) + t62 * t38;
t60 = t46 * mrSges(6,1) + t49 * mrSges(6,2) + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t59 = -pkin(2) * t47 + pkin(5);
t51 = cos(qJ(1));
t48 = sin(qJ(1));
t1 = (-m(4) * t59 + mrSges(3,1) * t47 + t41 * mrSges(4,1) + mrSges(3,2) * t50 + t42 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t64 * (-pkin(3) * t41 + t59) + t62 * t39 + t63 * t38 + (-m(2) - m(3)) * pkin(5)) * g(3) + (t61 * t48 - t60 * t51 - mrSges(1,2)) * g(2) + (t60 * t48 + t61 * t51 - mrSges(1,1)) * g(1);
U  = t1;
