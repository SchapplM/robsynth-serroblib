% Calculate potential energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:46
% EndTime: 2019-12-05 18:23:47
% DurationCPUTime: 0.29s
% Computational Cost: add. (90->49), mult. (118->48), div. (0->0), fcn. (95->8), ass. (0->25)
t70 = -m(4) * pkin(1) - mrSges(3,1) - mrSges(4,1);
t69 = mrSges(3,2) + mrSges(4,2);
t43 = qJ(2) + qJ(4);
t41 = sin(t43);
t42 = cos(t43);
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t68 = -t42 * mrSges(5,1) + t41 * mrSges(5,2) + t69 * t46 + t70 * t49 - mrSges(2,1);
t67 = -m(4) * qJ(3) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t47 = sin(qJ(1));
t65 = t47 * t41;
t45 = sin(qJ(5));
t64 = t47 * t45;
t48 = cos(qJ(5));
t63 = t47 * t48;
t51 = pkin(2) + pkin(1);
t62 = t49 * t51;
t50 = cos(qJ(1));
t61 = t50 * t41;
t60 = t50 * t45;
t59 = t50 * t48;
t44 = -pkin(3) - qJ(3);
t58 = t50 * t44 + t47 * t62;
t55 = -t47 * t44 + t50 * t62;
t1 = (-mrSges(1,3) - mrSges(2,3) - t69 * t49 + (m(6) * pkin(4) - mrSges(5,2) + mrSges(6,3)) * t42 + (-t48 * mrSges(6,1) + t45 * mrSges(6,2) - mrSges(5,1)) * t41 + ((-m(5) - m(6)) * t51 + t70) * t46) * g(3) + (-mrSges(1,2) - m(5) * t58 - m(6) * (pkin(4) * t65 + t58) - (t42 * t63 - t60) * mrSges(6,1) - (-t42 * t64 - t59) * mrSges(6,2) - mrSges(6,3) * t65 - t67 * t50 + t68 * t47) * g(2) + (-mrSges(1,1) - m(5) * t55 - m(6) * (pkin(4) * t61 + t55) - (t42 * t59 + t64) * mrSges(6,1) - (-t42 * t60 + t63) * mrSges(6,2) - mrSges(6,3) * t61 + t68 * t50 + t67 * t47) * g(1);
U = t1;
