% Calculate potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:10
% EndTime: 2019-03-09 01:37:10
% DurationCPUTime: 0.34s
% Computational Cost: add. (136->57), mult. (206->47), div. (0->0), fcn. (205->8), ass. (0->26)
t83 = m(6) + m(7);
t57 = sin(qJ(6));
t59 = cos(qJ(6));
t82 = m(7) * pkin(5) + t59 * mrSges(7,1) - t57 * mrSges(7,2) + mrSges(6,1);
t81 = m(7) * pkin(8) - mrSges(6,2) + mrSges(7,3);
t79 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t78 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t77 = t57 * mrSges(7,1) + t59 * mrSges(7,2) + t83 * pkin(7) - mrSges(5,2) + mrSges(6,3);
t58 = sin(qJ(5));
t60 = cos(qJ(5));
t76 = t81 * t58 + t82 * t60 + mrSges(5,1);
t75 = pkin(2) + pkin(6);
t74 = sin(qJ(1));
t61 = cos(qJ(1));
t73 = t61 * pkin(1) + t74 * qJ(2);
t72 = sin(pkin(9));
t71 = t61 * qJ(3) + t73;
t53 = t74 * pkin(1);
t69 = -t61 * qJ(2) + t53;
t68 = t74 * pkin(3) + t71;
t49 = t74 * qJ(3);
t65 = t49 + t53 + (-pkin(3) - qJ(2)) * t61;
t56 = cos(pkin(9));
t47 = t56 * t74 + t61 * t72;
t46 = t61 * t56 - t74 * t72;
t1 = (-m(4) * t75 - mrSges(3,1) + mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(5,3) - t82 * t58 + t81 * t60 + (-m(5) - t83) * (qJ(4) + t75) + (-m(2) - m(3)) * pkin(6)) * g(3) + (-mrSges(1,2) - m(3) * t69 - m(4) * (t49 + t69) - m(5) * t65 - t83 * (-t46 * pkin(4) + t65) + t79 * t74 - t78 * t61 + t77 * t47 + t76 * t46) * g(2) + (-m(3) * t73 - m(4) * t71 - m(5) * t68 - mrSges(1,1) - t83 * (t47 * pkin(4) + t68) + t78 * t74 + t79 * t61 - t76 * t47 + t77 * t46) * g(1);
U  = t1;
