% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:14
% EndTime: 2019-03-08 18:26:14
% DurationCPUTime: 0.24s
% Computational Cost: add. (93->62), mult. (175->78), div. (0->0), fcn. (177->6), ass. (0->25)
t62 = rSges(5,1) + pkin(3);
t65 = rSges(5,3) + qJ(4);
t51 = sin(qJ(1));
t64 = g(1) * t51;
t52 = cos(qJ(1));
t63 = g(2) * t52;
t50 = cos(pkin(4));
t61 = t50 * qJ(2) + pkin(5);
t47 = sin(pkin(6));
t60 = t51 * t47;
t49 = cos(pkin(6));
t59 = t51 * t49;
t58 = t52 * t47;
t57 = t52 * t49;
t48 = sin(pkin(4));
t56 = t51 * t48 * qJ(2) + t52 * pkin(1);
t55 = t48 * t47 * pkin(2) + t61;
t36 = -t50 * t57 + t60;
t37 = t50 * t58 + t59;
t45 = t51 * pkin(1);
t54 = t37 * pkin(2) + t36 * qJ(3) + t45;
t38 = t50 * t59 + t58;
t39 = -t50 * t60 + t57;
t53 = t39 * pkin(2) + t38 * qJ(3) + t56;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t52 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) + t52 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t39 * rSges(3,1) - t38 * rSges(3,2) + t56) + g(2) * (t37 * rSges(3,1) - t36 * rSges(3,2) + t45) + g(3) * (t50 * rSges(3,3) + t61) + (rSges(3,3) * t64 + g(3) * (rSges(3,1) * t47 + rSges(3,2) * t49) + (-rSges(3,3) - qJ(2)) * t63) * t48) - m(4) * (g(1) * (-t39 * rSges(4,2) + t38 * rSges(4,3) + t53) + g(2) * (-t37 * rSges(4,2) + t36 * rSges(4,3) + t54) + g(3) * (t50 * rSges(4,1) + t55) + (rSges(4,1) * t64 + g(3) * (-rSges(4,2) * t47 + (-rSges(4,3) - qJ(3)) * t49) + (-rSges(4,1) - qJ(2)) * t63) * t48) - m(5) * (g(1) * (t38 * rSges(5,2) + t39 * t65 + t53) + g(2) * (t36 * rSges(5,2) + t37 * t65 + t54) + g(3) * (t62 * t50 + t55) + (g(3) * ((-rSges(5,2) - qJ(3)) * t49 + t65 * t47) + t62 * t64 + (-qJ(2) - t62) * t63) * t48);
U  = t1;
