% Calculate potential energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:41
% EndTime: 2019-12-05 17:04:41
% DurationCPUTime: 0.15s
% Computational Cost: add. (86->43), mult. (66->32), div. (0->0), fcn. (36->8), ass. (0->17)
t50 = -m(5) - m(6);
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t49 = -t41 * mrSges(6,1) + t39 * mrSges(6,2) - mrSges(5,1);
t48 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t40 = sin(qJ(2));
t35 = t40 * pkin(2);
t42 = cos(qJ(2));
t47 = t42 * pkin(2) + pkin(1);
t46 = pkin(4) + qJ(1);
t38 = qJ(2) + qJ(3);
t34 = qJ(4) + t38;
t33 = cos(t38);
t32 = sin(t38);
t31 = cos(t34);
t30 = sin(t34);
t1 = (-m(4) * t46 - t39 * mrSges(6,1) - t41 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t50 * (pkin(5) + t46) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(4) * t35 - t40 * mrSges(3,1) - t32 * mrSges(4,1) - t42 * mrSges(3,2) - t33 * mrSges(4,2) - mrSges(1,2) - mrSges(2,2) + t50 * (pkin(3) * t32 + t35) + t48 * t31 + t49 * t30) * g(2) + (-m(3) * pkin(1) - m(4) * t47 - t42 * mrSges(3,1) - t33 * mrSges(4,1) + t40 * mrSges(3,2) + t32 * mrSges(4,2) - mrSges(1,1) - mrSges(2,1) + t50 * (pkin(3) * t33 + t47) + t49 * t31 - t48 * t30) * g(1);
U = t1;
