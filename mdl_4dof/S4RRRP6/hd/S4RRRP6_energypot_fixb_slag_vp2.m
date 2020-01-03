% Calculate potential energy for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:06
% EndTime: 2019-12-31 17:18:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (77->41), mult. (130->42), div. (0->0), fcn. (113->6), ass. (0->19)
t47 = cos(qJ(3));
t66 = -m(4) * pkin(2) - m(5) * (t47 * pkin(3) + pkin(2)) - mrSges(3,1);
t65 = m(4) * pkin(6) - m(5) * (-qJ(4) - pkin(6)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t64 = m(5) * pkin(3);
t63 = -mrSges(4,1) - mrSges(5,1);
t62 = mrSges(2,2) - mrSges(3,3);
t61 = mrSges(4,2) + mrSges(5,2);
t60 = -m(3) - m(4) - m(5);
t45 = sin(qJ(2));
t48 = cos(qJ(2));
t59 = -t65 * t45 + t66 * t48 - mrSges(2,1);
t44 = sin(qJ(3));
t46 = sin(qJ(1));
t58 = t46 * t44;
t57 = t46 * t48;
t49 = cos(qJ(1));
t56 = t49 * t44;
t55 = t49 * t48;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) + t60) * pkin(4) + t65 * t48 + (t61 * t44 + t63 * t47 + t66) * t45) * g(3) + (t56 * t64 - mrSges(1,2) + t60 * (t46 * pkin(1) - t49 * pkin(5)) - t62 * t49 + t63 * (t47 * t57 - t56) - t61 * (-t44 * t57 - t49 * t47) + t59 * t46) * g(2) + (-t58 * t64 - mrSges(1,1) + t60 * (t49 * pkin(1) + t46 * pkin(5)) + t62 * t46 + t63 * (t47 * t55 + t58) - t61 * (-t44 * t55 + t46 * t47) + t59 * t49) * g(1);
U = t1;
