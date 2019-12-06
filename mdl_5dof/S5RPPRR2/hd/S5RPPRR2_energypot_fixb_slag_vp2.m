% Calculate potential energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (102->48), mult. (110->37), div. (0->0), fcn. (79->8), ass. (0->18)
t57 = m(3) + m(4) + m(5) + m(6);
t42 = pkin(8) + qJ(4);
t36 = qJ(5) + t42;
t32 = sin(t36);
t33 = cos(t36);
t34 = sin(t42);
t35 = cos(t42);
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t53 = t43 * pkin(3);
t56 = m(5) * t53 + m(6) * (pkin(4) * t34 + t53) - mrSges(2,2) + mrSges(3,3) + t32 * mrSges(6,1) + t33 * mrSges(6,2) + t34 * mrSges(5,1) + t35 * mrSges(5,2) + t43 * mrSges(4,1) + t44 * mrSges(4,2);
t45 = -pkin(6) - qJ(3);
t55 = -m(4) * qJ(3) + m(5) * t45 + m(6) * (-pkin(7) + t45) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t54 = pkin(2) + pkin(5);
t51 = t44 * pkin(3) + t54;
t47 = cos(qJ(1));
t46 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t54 - t44 * mrSges(4,1) + t43 * mrSges(4,2) - m(5) * t51 - t35 * mrSges(5,1) + t34 * mrSges(5,2) - m(6) * (pkin(4) * t35 + t51) - t33 * mrSges(6,1) + t32 * mrSges(6,2) + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + (t57 * qJ(2) + t56) * t47 + (-t57 * pkin(1) + t55) * t46) * g(2) + (-mrSges(1,1) - t57 * (t47 * pkin(1) + t46 * qJ(2)) + t55 * t47 - t56 * t46) * g(1);
U = t1;
