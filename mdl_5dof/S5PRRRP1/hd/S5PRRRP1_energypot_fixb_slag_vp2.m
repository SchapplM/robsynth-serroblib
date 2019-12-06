% Calculate potential energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:42
% EndTime: 2019-12-05 16:39:43
% DurationCPUTime: 0.27s
% Computational Cost: add. (130->45), mult. (92->34), div. (0->0), fcn. (61->8), ass. (0->20)
t60 = -m(5) - m(6);
t59 = mrSges(5,2) + mrSges(6,2);
t58 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t57 = -m(4) + t60;
t48 = sin(qJ(4));
t49 = cos(qJ(4));
t56 = t60 * pkin(3) + t59 * t48 + t58 * t49 - mrSges(4,1);
t55 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t45 = sin(pkin(8));
t41 = t45 * pkin(1);
t46 = cos(pkin(8));
t42 = t46 * pkin(1);
t54 = pkin(5) + qJ(1);
t44 = pkin(8) + qJ(2);
t40 = qJ(3) + t44;
t39 = cos(t44);
t38 = sin(t44);
t36 = cos(t40);
t35 = sin(t40);
t1 = (-m(2) * qJ(1) - m(3) * t54 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - t59 * t49 + t58 * t48 + t57 * (pkin(6) + t54)) * g(3) + (-m(3) * t41 - t45 * mrSges(2,1) - t38 * mrSges(3,1) - t46 * mrSges(2,2) - t39 * mrSges(3,2) - mrSges(1,2) + t57 * (pkin(2) * t38 + t41) + t55 * t36 + t56 * t35) * g(2) + (-m(3) * t42 - t46 * mrSges(2,1) - t39 * mrSges(3,1) + t45 * mrSges(2,2) + t38 * mrSges(3,2) - mrSges(1,1) + t57 * (pkin(2) * t39 + t42) + t56 * t36 - t55 * t35) * g(1);
U = t1;
