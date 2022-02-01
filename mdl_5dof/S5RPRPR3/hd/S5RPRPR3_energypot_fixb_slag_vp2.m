% Calculate potential energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:32
% EndTime: 2022-01-23 09:20:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (152->50), mult. (112->42), div. (0->0), fcn. (85->10), ass. (0->23)
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t74 = -m(6) * pkin(4) - t58 * mrSges(6,1) + t56 * mrSges(6,2) - mrSges(5,1);
t73 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t72 = m(5) + m(6);
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t71 = t73 * t54 + t74 * t55 - mrSges(4,1);
t70 = -t56 * mrSges(6,1) - t58 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t57 = sin(qJ(1));
t50 = t57 * pkin(1);
t59 = cos(qJ(1));
t51 = t59 * pkin(1);
t68 = qJ(2) + pkin(5);
t53 = qJ(1) + pkin(8);
t47 = sin(t53);
t67 = pkin(2) * t47 + t50;
t48 = cos(t53);
t66 = pkin(2) * t48 + t51;
t49 = qJ(3) + t53;
t46 = cos(t49);
t45 = sin(t49);
t1 = (-m(2) * pkin(5) - m(3) * t68 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t74 * t54 - t73 * t55 + (-m(4) - t72) * (pkin(6) + t68)) * g(3) + (-m(3) * t50 - m(4) * t67 - t57 * mrSges(2,1) - t47 * mrSges(3,1) - t59 * mrSges(2,2) - t48 * mrSges(3,2) - mrSges(1,2) - t72 * (t45 * pkin(3) + t67) + (t72 * qJ(4) - t70) * t46 + t71 * t45) * g(2) + (-m(3) * t51 - m(4) * t66 - t59 * mrSges(2,1) - t48 * mrSges(3,1) + t57 * mrSges(2,2) + t47 * mrSges(3,2) - mrSges(1,1) - t72 * (t46 * pkin(3) + t45 * qJ(4) + t66) + t71 * t46 + t70 * t45) * g(1);
U = t1;
