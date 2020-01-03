% Calculate potential energy for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:35
% EndTime: 2019-12-31 16:43:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (83->35), mult. (81->29), div. (0->0), fcn. (56->6), ass. (0->15)
t48 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t47 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t46 = -m(4) - m(5);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t45 = -t47 * t36 + t48 * t38 - mrSges(3,1);
t44 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t37 = sin(qJ(1));
t32 = t37 * pkin(1);
t39 = cos(qJ(1));
t33 = t39 * pkin(1);
t34 = qJ(1) + pkin(6);
t31 = cos(t34);
t30 = sin(t34);
t1 = (-m(2) * pkin(4) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t47 * t38 + t48 * t36 + (-m(3) + t46) * (qJ(2) + pkin(4))) * g(3) + (-m(3) * t32 - t37 * mrSges(2,1) - mrSges(2,2) * t39 - mrSges(1,2) + t46 * (t30 * pkin(2) - pkin(5) * t31 + t32) - t44 * t31 + t45 * t30) * g(2) + (-m(3) * t33 - mrSges(2,1) * t39 + t37 * mrSges(2,2) - mrSges(1,1) + t46 * (t31 * pkin(2) + t30 * pkin(5) + t33) + t45 * t31 + t44 * t30) * g(1);
U = t1;
