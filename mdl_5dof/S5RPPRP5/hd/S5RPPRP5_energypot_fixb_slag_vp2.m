% Calculate potential energy for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:17
% EndTime: 2019-12-31 17:53:17
% DurationCPUTime: 0.36s
% Computational Cost: add. (112->50), mult. (201->51), div. (0->0), fcn. (188->6), ass. (0->25)
t102 = -mrSges(3,1) - mrSges(4,1);
t101 = mrSges(3,2) - mrSges(4,3);
t75 = sin(qJ(4));
t95 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t100 = t95 * t75;
t73 = sin(pkin(7));
t74 = cos(pkin(7));
t77 = cos(qJ(4));
t57 = t73 * t75 + t74 * t77;
t92 = t73 * t77;
t96 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t99 = t101 * t73 + t102 * t74 + t57 * t96 - t95 * t92 - mrSges(2,1);
t98 = -m(5) - m(6);
t94 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3) + mrSges(6,2);
t76 = sin(qJ(1));
t91 = t76 * t74;
t78 = cos(qJ(1));
t90 = t78 * t74;
t89 = t78 * pkin(1) + t76 * qJ(2);
t88 = qJ(3) * t73;
t87 = t76 * pkin(1) - t78 * qJ(2);
t86 = pkin(2) * t90 + t78 * t88 + t89;
t85 = t73 * pkin(2) - t74 * qJ(3) + pkin(5);
t81 = pkin(2) * t91 + t76 * t88 + t87;
t1 = (-m(4) * t85 - mrSges(1,3) - mrSges(2,3) + t98 * (t73 * pkin(3) + t85) - t101 * t74 + t102 * t73 + t96 * (-t74 * t75 + t92) + t95 * t57 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-m(3) * t87 - m(4) * t81 - mrSges(1,2) + t98 * (pkin(3) * t91 + t78 * pkin(6) + t81) + t91 * t100 - t94 * t78 + t99 * t76) * g(2) + (-m(3) * t89 - m(4) * t86 - mrSges(1,1) + t98 * (pkin(3) * t90 - t76 * pkin(6) + t86) + t90 * t100 + t94 * t76 + t99 * t78) * g(1);
U = t1;
