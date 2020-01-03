% Calculate potential energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:38
% EndTime: 2019-12-31 18:23:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (137->56), mult. (135->50), div. (0->0), fcn. (108->8), ass. (0->26)
t86 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t85 = -t60 * mrSges(6,1) - t63 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t82 = m(5) + m(6);
t58 = qJ(1) + pkin(8);
t53 = sin(t58);
t61 = sin(qJ(3));
t72 = qJ(4) * t61;
t64 = cos(qJ(3));
t78 = t53 * t64;
t81 = pkin(3) * t78 + t53 * t72;
t80 = t63 * mrSges(6,1) - t60 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t79 = t85 * t61 + t86 * t64 - mrSges(3,1);
t62 = sin(qJ(1));
t56 = t62 * pkin(1);
t65 = cos(qJ(1));
t57 = t65 * pkin(1);
t54 = cos(t58);
t77 = t54 * t64;
t59 = qJ(2) + pkin(5);
t73 = t53 * pkin(2) + t56;
t70 = t54 * pkin(2) + t53 * pkin(6) + t57;
t69 = -t54 * pkin(6) + t73;
t68 = pkin(3) * t77 + t54 * t72 + t70;
t1 = (-m(2) * pkin(5) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t82 * (t61 * pkin(3) + t59) + (-m(3) - m(4)) * t59 + (t82 * qJ(4) - t85) * t64 + (-m(6) * pkin(7) + t86) * t61) * g(3) + (-mrSges(1,2) - t62 * mrSges(2,1) - t65 * mrSges(2,2) - m(3) * t56 - m(4) * t69 - m(5) * (t69 + t81) - m(6) * (pkin(7) * t78 + t73 + t81) + (-m(6) * (-pkin(4) - pkin(6)) + t80) * t54 + t79 * t53) * g(2) + (-mrSges(1,1) - t65 * mrSges(2,1) + t62 * mrSges(2,2) - m(3) * t57 - m(4) * t70 - m(5) * t68 - m(6) * (pkin(7) * t77 + t68) + (-m(6) * pkin(4) - t80) * t53 + t79 * t54) * g(1);
U = t1;
