% Calculate potential energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:13
% EndTime: 2019-12-05 15:26:13
% DurationCPUTime: 0.36s
% Computational Cost: add. (133->74), mult. (141->91), div. (0->0), fcn. (121->8), ass. (0->28)
t74 = rSges(6,3) + pkin(6);
t73 = rSges(3,3) + pkin(5);
t51 = qJ(2) + pkin(8);
t48 = sin(t51);
t53 = cos(pkin(7));
t71 = t48 * t53;
t49 = cos(t51);
t70 = t49 * t53;
t52 = sin(pkin(7));
t55 = sin(qJ(5));
t69 = t52 * t55;
t57 = cos(qJ(5));
t68 = t52 * t57;
t67 = t53 * t55;
t66 = t53 * t57;
t58 = cos(qJ(2));
t47 = pkin(2) * t58 + pkin(1);
t54 = -qJ(3) - pkin(5);
t65 = t52 * t47 + t53 * t54;
t64 = qJ(4) * t48;
t56 = sin(qJ(2));
t63 = t56 * pkin(2) + qJ(1);
t44 = t53 * t47;
t62 = pkin(3) * t70 + t53 * t64 + t44;
t61 = t48 * pkin(3) + t63;
t60 = t65 + (pkin(3) * t49 + t64) * t52;
t59 = rSges(3,1) * t58 - rSges(3,2) * t56 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t53 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 + rSges(2,2) * t53) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t56 * rSges(3,1) + t58 * rSges(3,2) + qJ(1)) + (g(1) * t59 - g(2) * t73) * t53 + (g(1) * t73 + g(2) * t59) * t52) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t71 + t44) + g(2) * (-rSges(4,3) * t53 + t65) + g(3) * (rSges(4,1) * t48 + rSges(4,2) * t49 + t63) + (g(1) * (rSges(4,3) - t54) + g(2) * (rSges(4,1) * t49 - rSges(4,2) * t48)) * t52) - m(5) * (g(1) * (-rSges(5,2) * t70 + rSges(5,3) * t71 + t62) + g(2) * (-rSges(5,1) * t53 + t60) + g(3) * (-rSges(5,2) * t48 + (-rSges(5,3) - qJ(4)) * t49 + t61) + (g(1) * (rSges(5,1) - t54) + g(2) * (-rSges(5,2) * t49 + rSges(5,3) * t48)) * t52) - m(6) * (g(1) * ((t48 * t67 + t68) * rSges(6,1) + (t48 * t66 - t69) * rSges(6,2) + t62 + (pkin(4) - t54) * t52) + g(2) * (-t53 * pkin(4) + (t48 * t69 - t66) * rSges(6,1) + (t48 * t68 + t67) * rSges(6,2) + t60) + g(3) * (t74 * t48 + t61) + (g(3) * (-rSges(6,1) * t55 - rSges(6,2) * t57 - qJ(4)) + (g(1) * t53 + g(2) * t52) * t74) * t49);
U = t1;
