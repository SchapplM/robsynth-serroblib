% Calculate potential energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:05
% EndTime: 2019-03-09 02:05:05
% DurationCPUTime: 0.46s
% Computational Cost: add. (181->56), mult. (307->47), div. (0->0), fcn. (340->8), ass. (0->26)
t106 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t107 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t112 = m(6) + m(7);
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t119 = t112 * pkin(4) + t106 * t82 + t107 * t84 + mrSges(5,1);
t116 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t115 = -m(5) - t112;
t114 = t106 * t84 - t107 * t82 + mrSges(4,2) - mrSges(5,3);
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t113 = mrSges(4,1) + (t112 * pkin(8) - t116) * t83 + t119 * t85;
t111 = -mrSges(2,1) - mrSges(3,1);
t110 = mrSges(2,2) - mrSges(3,3);
t103 = sin(qJ(1));
t81 = -qJ(3) + pkin(6);
t86 = cos(qJ(1));
t98 = t86 * pkin(1) + t103 * qJ(2);
t97 = cos(pkin(9));
t96 = sin(pkin(9));
t94 = t86 * pkin(2) + t98;
t93 = t103 * pkin(1) - qJ(2) * t86;
t92 = t103 * pkin(2) + t93;
t71 = -t103 * t97 + t86 * t96;
t70 = -t103 * t96 - t86 * t97;
t1 = (-mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) - t112 * (t85 * pkin(8) + t81) + (-m(4) - m(5)) * t81 + (-m(2) - m(3)) * pkin(6) + t116 * t85 + t119 * t83) * g(3) + (-m(3) * t93 - m(4) * t92 + t111 * t103 - t110 * t86 - mrSges(1,2) + t115 * (-t71 * pkin(3) - t70 * pkin(7) + t92) - t114 * t70 + t113 * t71) * g(2) + (-m(3) * t98 - m(4) * t94 + t110 * t103 + t111 * t86 - mrSges(1,1) + t115 * (-t70 * pkin(3) + pkin(7) * t71 + t94) + t114 * t71 + t113 * t70) * g(1);
U  = t1;
