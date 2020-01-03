% Calculate potential energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:51
% EndTime: 2019-12-31 18:50:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (141->53), mult. (161->49), div. (0->0), fcn. (138->8), ass. (0->24)
t63 = cos(qJ(4));
t84 = -m(5) * pkin(3) - m(6) * (pkin(4) * t63 + pkin(3)) - mrSges(4,1);
t83 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t82 = m(6) * pkin(4);
t81 = -mrSges(5,1) - mrSges(6,1);
t80 = mrSges(5,2) + mrSges(6,2);
t79 = -m(4) - m(5) - m(6);
t78 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t56 = pkin(8) + qJ(3);
t53 = sin(t56);
t54 = cos(t56);
t57 = sin(pkin(8));
t58 = cos(pkin(8));
t77 = -m(3) * pkin(1) - t58 * mrSges(3,1) + t57 * mrSges(3,2) - t83 * t53 + t84 * t54 - mrSges(2,1);
t61 = sin(qJ(4));
t64 = cos(qJ(1));
t75 = t61 * t64;
t62 = sin(qJ(1));
t74 = t62 * t61;
t73 = t62 * t63;
t72 = t63 * t64;
t60 = -pkin(6) - qJ(2);
t50 = pkin(2) * t58 + pkin(1);
t1 = (-t57 * mrSges(3,1) - t58 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) + t79 * (t57 * pkin(2) + pkin(5)) + t83 * t54 + (t80 * t61 + t81 * t63 + t84) * t53) * g(3) + (t75 * t82 - mrSges(1,2) + t79 * (t62 * t50 + t64 * t60) + t81 * (t54 * t73 - t75) - t80 * (-t54 * t74 - t72) - t78 * t64 + t77 * t62) * g(2) + (-t74 * t82 - mrSges(1,1) + t79 * (t64 * t50 - t62 * t60) + t81 * (t54 * t72 + t74) - t80 * (-t54 * t75 + t73) + t78 * t62 + t77 * t64) * g(1);
U = t1;
