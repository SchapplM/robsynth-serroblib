% Calculate potential energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:34
% EndTime: 2019-12-31 17:39:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (135->49), mult. (115->42), div. (0->0), fcn. (98->8), ass. (0->22)
t70 = -m(5) - m(6);
t69 = -mrSges(3,1) - mrSges(4,1);
t55 = sin(qJ(5));
t56 = cos(qJ(5));
t68 = m(6) * pkin(4) + t56 * mrSges(6,1) - t55 * mrSges(6,2) + mrSges(5,1);
t67 = mrSges(3,2) - mrSges(4,3);
t66 = m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3);
t52 = sin(pkin(8));
t49 = t52 * pkin(1);
t53 = cos(pkin(8));
t50 = t53 * pkin(1);
t65 = cos(qJ(4));
t64 = sin(qJ(4));
t54 = pkin(5) + qJ(1);
t63 = pkin(8) + qJ(2);
t48 = cos(t63);
t60 = sin(t63);
t62 = t48 * pkin(2) + t60 * qJ(3) + t50;
t59 = t60 * pkin(2) - t48 * qJ(3) + t49;
t39 = t48 * t64 - t60 * t65;
t38 = -t48 * t65 - t60 * t64;
t1 = (-m(2) * qJ(1) + t55 * mrSges(6,1) + t56 * mrSges(6,2) - mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + mrSges(5,3) + (-m(3) - m(4)) * t54 + t70 * (-pkin(6) + t54)) * g(3) + (-m(3) * t49 - m(4) * t59 - t52 * mrSges(2,1) - t53 * mrSges(2,2) - mrSges(1,2) + t69 * t60 + t70 * (t60 * pkin(3) + t59) - t67 * t48 + t68 * t39 + t66 * t38) * g(2) + (-m(3) * t50 - m(4) * t62 - t53 * mrSges(2,1) + t52 * mrSges(2,2) - mrSges(1,1) + t70 * (t48 * pkin(3) + t62) + t67 * t60 + t69 * t48 - t66 * t39 + t68 * t38) * g(1);
U = t1;
