% Calculate potential energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:29
% EndTime: 2019-12-05 15:37:30
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->58), mult. (202->56), div. (0->0), fcn. (187->8), ass. (0->26)
t68 = sin(pkin(8));
t70 = cos(pkin(8));
t97 = -m(4) * pkin(2) - t70 * mrSges(4,1) + t68 * mrSges(4,2) - mrSges(3,1);
t96 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t95 = m(3) + m(4);
t94 = -m(5) - m(6);
t93 = -t68 * mrSges(4,1) - t70 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t91 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t90 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t73 = sin(qJ(2));
t74 = cos(qJ(2));
t89 = t96 * t73 + t97 * t74 - mrSges(2,1);
t88 = pkin(3) * t68;
t69 = sin(pkin(7));
t87 = t69 * t74;
t71 = cos(pkin(7));
t86 = t71 * t74;
t72 = -pkin(6) - qJ(3);
t85 = t72 * t73;
t82 = t71 * pkin(1) + t69 * pkin(5);
t67 = pkin(8) + qJ(4);
t65 = t69 * pkin(1);
t63 = cos(t67);
t62 = sin(t67);
t60 = t70 * pkin(3) + pkin(2);
t1 = (-mrSges(1,3) - mrSges(2,3) + t94 * (t73 * t60 + t74 * t72 + qJ(1)) + (-m(2) - t95) * qJ(1) - t96 * t74 + (t90 * t62 + t91 * t63 + t97) * t73) * g(3) + (-mrSges(1,2) + t94 * (-t69 * t85 + t60 * t87 + t65 + (-pkin(5) - t88) * t71) - t95 * t65 + t91 * (-t71 * t62 + t63 * t87) + t90 * (t62 * t87 + t71 * t63) + (t95 * pkin(5) - t93) * t71 + t89 * t69) * g(2) + (-mrSges(1,1) - t95 * t82 + t94 * (t60 * t86 + t69 * t88 - t71 * t85 + t82) + t91 * (t69 * t62 + t63 * t86) + t90 * (t62 * t86 - t69 * t63) + t93 * t69 + t89 * t71) * g(1);
U = t1;
