% Calculate potential energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:54
% EndTime: 2019-12-05 15:02:54
% DurationCPUTime: 0.36s
% Computational Cost: add. (133->74), mult. (141->91), div. (0->0), fcn. (121->8), ass. (0->28)
t74 = rSges(6,3) + pkin(6);
t51 = pkin(8) + qJ(3);
t48 = sin(t51);
t55 = cos(pkin(7));
t72 = t48 * t55;
t53 = sin(pkin(7));
t57 = sin(qJ(5));
t71 = t53 * t57;
t58 = cos(qJ(5));
t70 = t53 * t58;
t49 = cos(t51);
t69 = t55 * t49;
t68 = t55 * t57;
t67 = t55 * t58;
t54 = cos(pkin(8));
t47 = t54 * pkin(2) + pkin(1);
t56 = -pkin(5) - qJ(2);
t66 = t53 * t47 + t55 * t56;
t65 = qJ(4) * t48;
t64 = rSges(3,3) + qJ(2);
t52 = sin(pkin(8));
t63 = t52 * pkin(2) + qJ(1);
t42 = t55 * t47;
t62 = pkin(3) * t69 + t55 * t65 + t42;
t61 = t48 * pkin(3) + t63;
t60 = t66 + (pkin(3) * t49 + t65) * t53;
t59 = rSges(3,1) * t54 - rSges(3,2) * t52 + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t55 * rSges(2,1) - t53 * rSges(2,2)) + g(2) * (t53 * rSges(2,1) + t55 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t52 * rSges(3,1) + t54 * rSges(3,2) + qJ(1)) + (g(1) * t59 - g(2) * t64) * t55 + (g(1) * t64 + g(2) * t59) * t53) - m(4) * (g(1) * (rSges(4,1) * t69 - rSges(4,2) * t72 + t42) + g(2) * (-t55 * rSges(4,3) + t66) + g(3) * (t48 * rSges(4,1) + t49 * rSges(4,2) + t63) + (g(1) * (rSges(4,3) - t56) + g(2) * (rSges(4,1) * t49 - rSges(4,2) * t48)) * t53) - m(5) * (g(1) * (-rSges(5,2) * t69 + rSges(5,3) * t72 + t62) + g(2) * (-t55 * rSges(5,1) + t60) + g(3) * (-t48 * rSges(5,2) + (-rSges(5,3) - qJ(4)) * t49 + t61) + (g(1) * (rSges(5,1) - t56) + g(2) * (-rSges(5,2) * t49 + rSges(5,3) * t48)) * t53) - m(6) * (g(1) * ((t48 * t68 + t70) * rSges(6,1) + (t48 * t67 - t71) * rSges(6,2) + t62 + (pkin(4) - t56) * t53) + g(2) * (-t55 * pkin(4) + (t48 * t71 - t67) * rSges(6,1) + (t48 * t70 + t68) * rSges(6,2) + t60) + g(3) * (t74 * t48 + t61) + (g(3) * (-rSges(6,1) * t57 - rSges(6,2) * t58 - qJ(4)) + (g(1) * t55 + g(2) * t53) * t74) * t49);
U = t1;
