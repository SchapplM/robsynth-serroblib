% Calculate potential energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:56
% EndTime: 2019-12-05 17:37:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (80->49), mult. (98->36), div. (0->0), fcn. (67->6), ass. (0->18)
t51 = -m(5) - m(6);
t35 = qJ(4) + qJ(5);
t28 = sin(t35);
t29 = cos(t35);
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t50 = -t28 * mrSges(6,1) - t38 * mrSges(5,2) - t29 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t36;
t49 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t48 = pkin(2) + pkin(5);
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t46 = t39 * pkin(1) + t37 * qJ(2);
t45 = pkin(3) + t48;
t33 = t37 * pkin(1);
t43 = -t39 * qJ(2) + t33;
t40 = -pkin(7) - pkin(6);
t30 = t37 * qJ(3);
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t48 - mrSges(4,1) - m(5) * t45 - t38 * mrSges(5,1) + t36 * mrSges(5,2) - m(6) * (t38 * pkin(4) + t45) - t29 * mrSges(6,1) + t28 * mrSges(6,2) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) - m(3) * t43 - m(4) * (t30 + t43) + t51 * (t30 + t33) + (-m(5) * (pkin(6) - qJ(2)) - m(6) * (-qJ(2) - t40) + t49) * t39 + t50 * t37) * g(2) + (-m(3) * t46 - mrSges(1,1) + (-m(4) + t51) * (t39 * qJ(3) + t46) + t50 * t39 + (m(5) * pkin(6) - m(6) * t40 - t49) * t37) * g(1);
U = t1;
