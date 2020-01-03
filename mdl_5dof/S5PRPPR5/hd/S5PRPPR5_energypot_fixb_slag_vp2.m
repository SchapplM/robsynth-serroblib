% Calculate potential energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:58
% EndTime: 2019-12-31 17:37:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (120->52), mult. (222->53), div. (0->0), fcn. (216->8), ass. (0->26)
t112 = -mrSges(3,1) - mrSges(4,1);
t111 = mrSges(3,2) - mrSges(4,3);
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t103 = -m(6) * pkin(4) - t84 * mrSges(6,1) + t82 * mrSges(6,2) - mrSges(5,1);
t105 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t78 = sin(pkin(8));
t80 = cos(pkin(8));
t83 = sin(qJ(2));
t85 = cos(qJ(2));
t107 = t85 * t78 - t83 * t80;
t62 = t83 * t78 + t85 * t80;
t110 = t62 * t103 + t107 * t105 + t111 * t83 + t112 * t85 - mrSges(2,1);
t108 = -m(5) - m(6);
t109 = -t82 * mrSges(6,1) - t84 * mrSges(6,2) + t108 * qJ(4) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t79 = sin(pkin(7));
t100 = t79 * t85;
t81 = cos(pkin(7));
t99 = t81 * t85;
t96 = t81 * pkin(1) + t79 * pkin(5);
t95 = qJ(3) * t83;
t94 = t79 * pkin(1) - t81 * pkin(5);
t93 = pkin(2) * t99 + t81 * t95 + t96;
t92 = t83 * pkin(2) - t85 * qJ(3) + qJ(1);
t89 = pkin(2) * t100 + t79 * t95 + t94;
t1 = (-m(4) * t92 - mrSges(1,3) - mrSges(2,3) + t108 * (t83 * pkin(3) + t92) - t111 * t85 + t112 * t83 - t103 * t107 + t105 * t62 + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t94 - m(4) * t89 - mrSges(1,2) + t108 * (pkin(3) * t100 + t89)) * g(2) + (-m(3) * t96 - m(4) * t93 - mrSges(1,1) + t108 * (pkin(3) * t99 + t93)) * g(1) + (t110 * g(1) + t109 * g(2)) * t81 + (-t109 * g(1) + t110 * g(2)) * t79;
U = t1;
