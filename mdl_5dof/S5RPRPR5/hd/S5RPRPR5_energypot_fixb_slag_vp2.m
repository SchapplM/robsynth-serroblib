% Calculate potential energy for
% S5RPRPR5
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:46
% EndTime: 2019-12-05 17:55:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (152->49), mult. (187->41), div. (0->0), fcn. (168->10), ass. (0->21)
t55 = qJ(3) + pkin(9);
t49 = qJ(5) + t55;
t44 = sin(t49);
t45 = cos(t49);
t61 = cos(qJ(3));
t46 = t61 * pkin(3) + pkin(2);
t47 = sin(t55);
t48 = cos(t55);
t59 = sin(qJ(3));
t79 = -mrSges(3,1) - m(4) * pkin(2) - t61 * mrSges(4,1) + t59 * mrSges(4,2) - m(5) * t46 - t48 * mrSges(5,1) + t47 * mrSges(5,2) - m(6) * (pkin(4) * t48 + t46) - t45 * mrSges(6,1) + t44 * mrSges(6,2);
t58 = -qJ(4) - pkin(6);
t78 = mrSges(3,2) - m(4) * pkin(6) + m(5) * t58 + m(6) * (-pkin(7) + t58) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t77 = t78 * t56 + t79 * t57 - mrSges(2,1);
t74 = m(3) + m(4) + m(5) + m(6);
t71 = t59 * pkin(3);
t73 = -m(5) * t71 - t47 * mrSges(5,1) - t48 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - m(6) * (pkin(4) * t47 + t71) - t44 * mrSges(6,1) - t45 * mrSges(6,2) - t59 * mrSges(4,1) - t61 * mrSges(4,2);
t62 = cos(qJ(1));
t60 = sin(qJ(1));
t1 = (-mrSges(1,3) - t74 * (t62 * pkin(1) + t60 * qJ(2)) + t77 * t62 + t73 * t60) * g(3) + (-mrSges(1,2) + (t74 * pkin(1) - t77) * t60 + (-t74 * qJ(2) + t73) * t62) * g(2) + (-mrSges(1,1) - mrSges(2,3) + (-m(2) - t74) * pkin(5) - t78 * t57 + t79 * t56) * g(1);
U = t1;
