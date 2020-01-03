% Calculate potential energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:06
% EndTime: 2019-12-31 17:20:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (79->46), mult. (143->49), div. (0->0), fcn. (130->6), ass. (0->19)
t72 = -m(4) - m(5);
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t71 = -t55 * mrSges(3,1) + t52 * mrSges(3,2) - mrSges(2,1);
t70 = mrSges(2,2) - mrSges(3,3);
t69 = -mrSges(4,3) - mrSges(5,2);
t68 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t67 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t53 = sin(qJ(1));
t66 = t53 * t52;
t65 = t53 * t55;
t56 = cos(qJ(1));
t64 = t56 * t52;
t63 = t56 * t55;
t62 = t56 * pkin(1) + t53 * pkin(5);
t61 = t53 * pkin(1) - t56 * pkin(5);
t54 = cos(qJ(3));
t51 = sin(qJ(3));
t1 = (-mrSges(1,3) - mrSges(2,3) + t72 * (t52 * pkin(2) - t55 * pkin(6) + pkin(4)) + (-m(2) - m(3)) * pkin(4) + (-mrSges(3,2) - t69) * t55 + (t67 * t51 + t68 * t54 - mrSges(3,1)) * t52) * g(3) + (-m(3) * t61 - mrSges(1,2) + t69 * t66 + t72 * (pkin(2) * t65 + pkin(6) * t66 + t61) - t70 * t56 + t71 * t53 + t68 * (-t56 * t51 + t54 * t65) + t67 * (t51 * t65 + t56 * t54)) * g(2) + (-m(3) * t62 - mrSges(1,1) + t69 * t64 + t72 * (pkin(2) * t63 + pkin(6) * t64 + t62) + t71 * t56 + t70 * t53 + t68 * (t53 * t51 + t54 * t63) + t67 * (t51 * t63 - t53 * t54)) * g(1);
U = t1;
