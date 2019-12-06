% Calculate potential energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:16
% EndTime: 2019-12-05 18:14:16
% DurationCPUTime: 0.22s
% Computational Cost: add. (134->50), mult. (81->41), div. (0->0), fcn. (50->10), ass. (0->23)
t59 = -m(5) - m(6);
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t58 = m(6) * pkin(4) + t47 * mrSges(6,1) - t45 * mrSges(6,2) + mrSges(5,1);
t57 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t46 = sin(qJ(1));
t56 = t46 * pkin(1);
t48 = cos(qJ(1));
t43 = t48 * pkin(1);
t55 = qJ(2) + pkin(5);
t44 = qJ(1) + pkin(9);
t40 = cos(t44);
t54 = pkin(2) * t40 + t43;
t53 = pkin(6) + t55;
t41 = qJ(3) + t44;
t39 = sin(t44);
t51 = -pkin(2) * t39 - t56;
t38 = qJ(4) + t41;
t37 = cos(t41);
t36 = sin(t41);
t34 = cos(t38);
t33 = sin(t38);
t1 = (-m(3) * t43 - m(4) * t54 - t48 * mrSges(2,1) - t40 * mrSges(3,1) - t37 * mrSges(4,1) + t46 * mrSges(2,2) + t39 * mrSges(3,2) + t36 * mrSges(4,2) - mrSges(1,3) + t59 * (pkin(3) * t37 + t54) - t58 * t34 + t57 * t33) * g(3) + (m(3) * t56 - m(4) * t51 + t46 * mrSges(2,1) + t39 * mrSges(3,1) + t36 * mrSges(4,1) + t48 * mrSges(2,2) + t40 * mrSges(3,2) + t37 * mrSges(4,2) - mrSges(1,2) + t59 * (-pkin(3) * t36 + t51) + t57 * t34 + t58 * t33) * g(2) + (-m(2) * pkin(5) - m(3) * t55 - m(4) * t53 - t45 * mrSges(6,1) - t47 * mrSges(6,2) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t59 * (pkin(7) + t53)) * g(1);
U = t1;
