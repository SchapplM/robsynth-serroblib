% Calculate potential energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:43
% DurationCPUTime: 0.27s
% Computational Cost: add. (128->48), mult. (122->39), div. (0->0), fcn. (91->6), ass. (0->20)
t57 = sin(qJ(3));
t58 = cos(qJ(3));
t76 = pkin(3) * t58 + qJ(4) * t57;
t75 = -m(6) * pkin(4) - mrSges(4,1) - mrSges(5,1) - mrSges(6,1);
t74 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t73 = -m(5) - m(6);
t53 = pkin(7) + qJ(2);
t48 = sin(t53);
t72 = t76 * t48;
t70 = t74 * t57 + t75 * t58 - mrSges(3,1);
t69 = mrSges(5,2) + mrSges(4,3) - mrSges(3,2) - mrSges(6,3);
t54 = sin(pkin(7));
t50 = t54 * pkin(1);
t55 = cos(pkin(7));
t51 = t55 * pkin(1);
t56 = pkin(5) + qJ(1);
t67 = t48 * pkin(2) + t50;
t49 = cos(t53);
t64 = -t49 * pkin(6) + t67;
t1 = (-m(2) * qJ(1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t73 * (t57 * pkin(3) - t58 * qJ(4) + t56) + (-m(3) - m(4)) * t56 - t74 * t58 + t75 * t57) * g(3) + (-mrSges(1,2) - t54 * mrSges(2,1) - t55 * mrSges(2,2) - m(3) * t50 - m(4) * t64 - m(5) * (t64 + t72) - m(6) * (t67 + t72) + (-m(6) * (-pkin(6) + qJ(5)) + t69) * t49 + t70 * t48) * g(2) + (-m(3) * t51 - t55 * mrSges(2,1) + t54 * mrSges(2,2) - mrSges(1,1) + (m(6) * qJ(5) - t69) * t48 + (-m(4) + t73) * (t49 * pkin(2) + t48 * pkin(6) + t51) + (t73 * t76 + t70) * t49) * g(1);
U = t1;
