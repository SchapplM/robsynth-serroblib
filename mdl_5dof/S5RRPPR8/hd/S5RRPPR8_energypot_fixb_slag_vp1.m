% Calculate potential energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:02
% EndTime: 2019-12-31 19:38:02
% DurationCPUTime: 0.34s
% Computational Cost: add. (120->78), mult. (179->96), div. (0->0), fcn. (167->8), ass. (0->27)
t52 = sin(qJ(2));
t70 = t52 * pkin(2) + pkin(5);
t49 = sin(pkin(8));
t69 = t49 * t52;
t53 = sin(qJ(1));
t68 = t52 * t53;
t54 = cos(qJ(2));
t67 = t53 * t54;
t55 = cos(qJ(1));
t66 = t55 * pkin(1) + t53 * pkin(6);
t65 = -pkin(7) - qJ(4) - rSges(6,3);
t64 = qJ(3) * t52;
t63 = -qJ(4) - rSges(5,3);
t46 = t53 * pkin(1);
t62 = pkin(2) * t67 + t53 * t64 + t46;
t61 = t66 + (pkin(2) * t54 + t64) * t55;
t50 = cos(pkin(8));
t60 = -t49 * t54 + t50 * t52;
t59 = t50 * t54 + t69;
t58 = g(1) * t61 + g(2) * t62;
t57 = t59 * rSges(5,1) + t60 * rSges(5,2) + t54 * pkin(3);
t41 = pkin(4) * t50 + pkin(3);
t48 = pkin(8) + qJ(5);
t42 = sin(t48);
t43 = cos(t48);
t56 = t54 * t41 + pkin(4) * t69 + (t42 * t52 + t43 * t54) * rSges(6,1) + (-t42 * t54 + t43 * t52) * rSges(6,2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t55 - t53 * rSges(2,2)) + g(2) * (t53 * rSges(2,1) + rSges(2,2) * t55) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t53 * rSges(3,3) + t66) + g(2) * (rSges(3,1) * t67 - rSges(3,2) * t68 + t46) + g(3) * (rSges(3,1) * t52 + rSges(3,2) * t54 + pkin(5)) + (g(1) * (rSges(3,1) * t54 - rSges(3,2) * t52) + g(2) * (-rSges(3,3) - pkin(6))) * t55) - m(4) * (g(1) * (t53 * rSges(4,2) + t61) + g(2) * (rSges(4,1) * t67 + rSges(4,3) * t68 + t62) + g(3) * (rSges(4,1) * t52 + (-rSges(4,3) - qJ(3)) * t54 + t70) + (g(1) * (rSges(4,1) * t54 + rSges(4,3) * t52) + g(2) * (-rSges(4,2) - pkin(6))) * t55) - m(5) * (g(3) * (t60 * rSges(5,1) - t59 * rSges(5,2) + t52 * pkin(3) - t54 * qJ(3) + t70) + (g(1) * t63 + g(2) * t57) * t53 + (g(1) * t57 + g(2) * (-pkin(6) - t63)) * t55 + t58) - m(6) * (g(3) * ((t43 * rSges(6,1) - t42 * rSges(6,2) + t41) * t52 + (-t42 * rSges(6,1) - t43 * rSges(6,2) - t49 * pkin(4) - qJ(3)) * t54 + t70) + (g(1) * t65 + g(2) * t56) * t53 + (g(1) * t56 + g(2) * (-pkin(6) - t65)) * t55 + t58);
U = t1;
