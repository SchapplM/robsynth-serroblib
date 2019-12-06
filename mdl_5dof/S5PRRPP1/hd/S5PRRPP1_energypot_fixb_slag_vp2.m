% Calculate potential energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:08
% EndTime: 2019-12-05 16:06:08
% DurationCPUTime: 0.31s
% Computational Cost: add. (137->48), mult. (111->38), div. (0->0), fcn. (80->8), ass. (0->22)
t66 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t65 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t64 = -m(3) - m(4);
t63 = -m(5) - m(6);
t48 = qJ(3) + pkin(8);
t41 = sin(t48);
t43 = cos(t48);
t53 = sin(qJ(3));
t54 = cos(qJ(3));
t62 = -m(4) * pkin(2) - t54 * mrSges(4,1) + t53 * mrSges(4,2) - t65 * t41 + t66 * t43 - mrSges(3,1);
t61 = m(4) * pkin(6) - mrSges(3,2) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3);
t49 = sin(pkin(7));
t44 = t49 * pkin(1);
t50 = cos(pkin(7));
t45 = t50 * pkin(1);
t52 = pkin(5) + qJ(1);
t51 = -qJ(4) - pkin(6);
t47 = pkin(7) + qJ(2);
t42 = cos(t47);
t40 = sin(t47);
t39 = t54 * pkin(3) + pkin(2);
t1 = (-m(2) * qJ(1) - t53 * mrSges(4,1) - t54 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t63 * (t53 * pkin(3) + t52) + t64 * t52 + t65 * t43 + t66 * t41) * g(3) + (-t49 * mrSges(2,1) - t50 * mrSges(2,2) - mrSges(1,2) + t63 * (t40 * t39 + t42 * t51 + t44) + t64 * t44 + t61 * t42 + t62 * t40) * g(2) + (-t50 * mrSges(2,1) + t49 * mrSges(2,2) - mrSges(1,1) + t63 * (t42 * t39 - t40 * t51 + t45) + t64 * t45 + t62 * t42 - t61 * t40) * g(1);
U = t1;
