% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:19
% EndTime: 2019-03-08 18:26:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (93->47), mult. (180->49), div. (0->0), fcn. (177->6), ass. (0->26)
t77 = m(4) + m(5);
t76 = mrSges(4,1) + mrSges(5,1) + mrSges(3,3);
t75 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2);
t74 = -m(5) * pkin(3) - t76;
t73 = -m(5) * qJ(4) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t57 = cos(pkin(4));
t72 = t57 * qJ(2) + pkin(5);
t54 = sin(pkin(6));
t58 = sin(qJ(1));
t71 = t58 * t54;
t56 = cos(pkin(6));
t69 = t58 * t56;
t59 = cos(qJ(1));
t68 = t59 * t54;
t66 = t59 * t56;
t55 = sin(pkin(4));
t64 = qJ(2) * t55;
t65 = t59 * pkin(1) + t58 * t64;
t62 = t59 * t64;
t43 = -t57 * t66 + t71;
t44 = t57 * t68 + t69;
t52 = t58 * pkin(1);
t61 = t44 * pkin(2) + t43 * qJ(3) + t52;
t46 = -t57 * t71 + t66;
t45 = t57 * t69 + t68;
t1 = (-m(2) * pkin(5) - m(3) * t72 - mrSges(1,3) - mrSges(2,3) - t77 * (t55 * t54 * pkin(2) + t72) + t74 * t57 + ((t77 * qJ(3) - t75) * t56 + t73 * t54) * t55) * g(3) + (-mrSges(1,2) - t58 * mrSges(2,1) - m(3) * (t52 - t62) - m(4) * (t61 - t62) - m(5) * t61 + t73 * t44 + t75 * t43 + (-mrSges(2,2) + (-m(5) * (-pkin(3) - qJ(2)) + t76) * t55) * t59) * g(2) + (-m(3) * t65 - t59 * mrSges(2,1) - mrSges(1,1) - t77 * (t46 * pkin(2) + t45 * qJ(3) + t65) + t73 * t46 + t75 * t45 + (t74 * t55 + mrSges(2,2)) * t58) * g(1);
U  = t1;
