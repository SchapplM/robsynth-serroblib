% Calculate potential energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:51
% EndTime: 2019-12-31 19:34:51
% DurationCPUTime: 0.32s
% Computational Cost: add. (133->74), mult. (141->91), div. (0->0), fcn. (121->8), ass. (0->28)
t74 = rSges(6,3) + pkin(7);
t73 = rSges(3,3) + pkin(6);
t54 = sin(qJ(2));
t71 = t54 * pkin(2) + pkin(5);
t51 = qJ(2) + pkin(8);
t48 = sin(t51);
t58 = cos(qJ(1));
t70 = t48 * t58;
t49 = cos(t51);
t69 = t49 * t58;
t53 = sin(qJ(5));
t68 = t53 * t58;
t55 = sin(qJ(1));
t67 = t55 * t53;
t56 = cos(qJ(5));
t66 = t55 * t56;
t65 = t56 * t58;
t57 = cos(qJ(2));
t47 = pkin(2) * t57 + pkin(1);
t52 = -qJ(3) - pkin(6);
t64 = t55 * t47 + t58 * t52;
t63 = qJ(4) * t48;
t62 = t48 * pkin(3) + t71;
t44 = t58 * t47;
t61 = pkin(3) * t69 + t58 * t63 + t44;
t60 = t64 + (pkin(3) * t49 + t63) * t55;
t59 = rSges(3,1) * t57 - rSges(3,2) * t54 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t58 - t55 * rSges(2,2)) + g(2) * (t55 * rSges(2,1) + rSges(2,2) * t58) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t54 + t57 * rSges(3,2) + pkin(5)) + (g(1) * t59 - g(2) * t73) * t58 + (g(1) * t73 + g(2) * t59) * t55) - m(4) * (g(1) * (rSges(4,1) * t69 - rSges(4,2) * t70 + t44) + g(2) * (-rSges(4,3) * t58 + t64) + g(3) * (rSges(4,1) * t48 + rSges(4,2) * t49 + t71) + (g(1) * (rSges(4,3) - t52) + g(2) * (rSges(4,1) * t49 - rSges(4,2) * t48)) * t55) - m(5) * (g(1) * (-rSges(5,2) * t69 + rSges(5,3) * t70 + t61) + g(2) * (-rSges(5,1) * t58 + t60) + g(3) * (-rSges(5,2) * t48 + (-rSges(5,3) - qJ(4)) * t49 + t62) + (g(1) * (rSges(5,1) - t52) + g(2) * (-rSges(5,2) * t49 + rSges(5,3) * t48)) * t55) - m(6) * (g(1) * ((t48 * t68 + t66) * rSges(6,1) + (t48 * t65 - t67) * rSges(6,2) + t61 + (pkin(4) - t52) * t55) + g(2) * (-t58 * pkin(4) + (t48 * t67 - t65) * rSges(6,1) + (t48 * t66 + t68) * rSges(6,2) + t60) + g(3) * (t74 * t48 + t62) + (g(3) * (-rSges(6,1) * t53 - rSges(6,2) * t56 - qJ(4)) + (g(1) * t58 + g(2) * t55) * t74) * t49);
U = t1;
