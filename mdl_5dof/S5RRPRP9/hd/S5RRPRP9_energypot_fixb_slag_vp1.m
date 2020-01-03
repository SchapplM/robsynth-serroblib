% Calculate potential energy for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:37
% EndTime: 2019-12-31 20:05:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (149->81), mult. (195->101), div. (0->0), fcn. (187->8), ass. (0->31)
t88 = rSges(6,1) + pkin(4);
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t87 = g(1) * t69 + g(2) * t67;
t86 = rSges(4,3) + qJ(3);
t85 = rSges(6,3) + qJ(5);
t66 = sin(qJ(2));
t82 = rSges(3,2) * t66;
t63 = sin(pkin(8));
t81 = t63 * t69;
t80 = t67 * t63;
t68 = cos(qJ(2));
t79 = t67 * t68;
t78 = t69 * t68;
t75 = t69 * pkin(1) + t67 * pkin(6);
t64 = cos(pkin(8));
t55 = pkin(3) * t64 + pkin(2);
t65 = -pkin(7) - qJ(3);
t73 = t66 * t55 + t68 * t65 + pkin(5);
t60 = t67 * pkin(1);
t72 = -t69 * pkin(6) + t60;
t71 = pkin(3) * t80 + t55 * t78 + t75;
t70 = -pkin(3) * t81 + t55 * t79 + t72;
t62 = pkin(8) + qJ(4);
t58 = cos(t62);
t57 = sin(t62);
t49 = t67 * t57 + t58 * t78;
t48 = t57 * t78 - t67 * t58;
t47 = -t57 * t69 + t58 * t79;
t46 = t57 * t79 + t58 * t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t69 - t67 * rSges(2,2)) + g(2) * (t67 * rSges(2,1) + rSges(2,2) * t69) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t67 * rSges(3,3) + t75) + g(2) * (rSges(3,1) * t79 - t67 * t82 + t60) + g(3) * (rSges(3,1) * t66 + rSges(3,2) * t68 + pkin(5)) + (g(1) * (rSges(3,1) * t68 - t82) + g(2) * (-rSges(3,3) - pkin(6))) * t69) - m(4) * (g(1) * (pkin(2) * t78 + (t64 * t78 + t80) * rSges(4,1) + (-t63 * t78 + t67 * t64) * rSges(4,2) + t75) + g(2) * (pkin(2) * t79 + (t64 * t79 - t81) * rSges(4,1) + (-t63 * t79 - t64 * t69) * rSges(4,2) + t72) + g(3) * (-t86 * t68 + pkin(5)) + (g(3) * (rSges(4,1) * t64 - rSges(4,2) * t63 + pkin(2)) + t87 * t86) * t66) - m(5) * (g(1) * (t49 * rSges(5,1) - t48 * rSges(5,2) + t71) + g(2) * (t47 * rSges(5,1) - t46 * rSges(5,2) + t70) + g(3) * (-rSges(5,3) * t68 + t73) + (g(3) * (rSges(5,1) * t58 - rSges(5,2) * t57) + t87 * (rSges(5,3) - t65)) * t66) - m(6) * (g(1) * (t85 * t48 + t88 * t49 + t71) + g(2) * (t85 * t46 + t88 * t47 + t70) + g(3) * (-rSges(6,2) * t68 + t73) + (g(3) * (t85 * t57 + t88 * t58) + t87 * (rSges(6,2) - t65)) * t66);
U = t1;
