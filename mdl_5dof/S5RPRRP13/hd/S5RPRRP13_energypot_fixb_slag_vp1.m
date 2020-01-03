% Calculate potential energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:27
% EndTime: 2019-12-31 18:58:28
% DurationCPUTime: 0.28s
% Computational Cost: add. (98->71), mult. (156->86), div. (0->0), fcn. (144->6), ass. (0->29)
t76 = rSges(6,1) + pkin(4);
t75 = rSges(6,3) + qJ(5);
t74 = pkin(2) + pkin(5);
t56 = sin(qJ(1));
t73 = g(1) * t56;
t59 = cos(qJ(1));
t72 = g(2) * t59;
t58 = cos(qJ(3));
t71 = rSges(4,2) * t58;
t55 = sin(qJ(3));
t70 = t55 * t56;
t54 = sin(qJ(4));
t69 = t56 * t54;
t57 = cos(qJ(4));
t68 = t56 * t57;
t67 = t59 * t54;
t66 = t59 * t57;
t65 = t59 * pkin(1) + t56 * qJ(2);
t50 = t56 * pkin(1);
t64 = t56 * pkin(6) + t50;
t63 = t59 * pkin(6) + t65;
t62 = t58 * pkin(3) + t55 * pkin(7) + t74;
t61 = pkin(3) * t70 + t63;
t60 = t64 + (-pkin(3) * t55 + t58 * pkin(7) - qJ(2)) * t59;
t43 = -t55 * t66 + t69;
t42 = t55 * t67 + t68;
t41 = t55 * t68 + t67;
t40 = t55 * t69 - t66;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t59 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t59 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t59 * rSges(3,2) + t56 * rSges(3,3) + t65) + g(2) * (-t56 * rSges(3,2) + t50 + (-rSges(3,3) - qJ(2)) * t59) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * (rSges(4,1) * t70 + t56 * t71 + t63) + g(2) * (t56 * rSges(4,3) + t64) + g(3) * (t58 * rSges(4,1) - t55 * rSges(4,2) + t74) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t55 - qJ(2) - t71)) * t59) - m(5) * (g(1) * (t41 * rSges(5,1) - t40 * rSges(5,2) + t61) + g(2) * (t43 * rSges(5,1) + t42 * rSges(5,2) + t60) + g(3) * (t55 * rSges(5,3) + t62) + (rSges(5,3) * t72 + g(3) * (rSges(5,1) * t57 - rSges(5,2) * t54) + (-rSges(5,3) - pkin(7)) * t73) * t58) - m(6) * (g(1) * (t75 * t40 + t76 * t41 + t61) + g(2) * (-t75 * t42 + t76 * t43 + t60) + g(3) * (t55 * rSges(6,2) + t62) + (rSges(6,2) * t72 + g(3) * (t75 * t54 + t76 * t57) + (-rSges(6,2) - pkin(7)) * t73) * t58);
U = t1;
