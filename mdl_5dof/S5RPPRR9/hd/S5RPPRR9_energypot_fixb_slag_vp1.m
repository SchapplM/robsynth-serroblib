% Calculate potential energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:15
% DurationCPUTime: 0.28s
% Computational Cost: add. (114->66), mult. (177->85), div. (0->0), fcn. (189->8), ass. (0->24)
t69 = rSges(6,3) + pkin(7);
t51 = cos(qJ(4));
t68 = t51 * pkin(4);
t67 = rSges(5,3) + pkin(6);
t65 = cos(qJ(1));
t64 = sin(qJ(1));
t48 = sin(qJ(5));
t63 = t48 * t51;
t50 = cos(qJ(5));
t62 = t50 * t51;
t61 = pkin(5) - qJ(3);
t60 = t65 * pkin(1) + t64 * qJ(2);
t59 = cos(pkin(8));
t58 = sin(pkin(8));
t57 = t65 * pkin(2) + t60;
t38 = -t64 * t58 - t65 * t59;
t56 = -t38 * pkin(3) + t57;
t55 = t64 * pkin(1) - t65 * qJ(2);
t49 = sin(qJ(4));
t54 = -rSges(5,1) * t51 + rSges(5,2) * t49;
t53 = t64 * pkin(2) + t55;
t39 = t65 * t58 - t64 * t59;
t52 = -t39 * pkin(3) + t53;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t65 * rSges(2,1) - t64 * rSges(2,2)) + g(2) * (t64 * rSges(2,1) + t65 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t65 * rSges(3,1) + t64 * rSges(3,3) + t60) + g(2) * (t64 * rSges(3,1) - t65 * rSges(3,3) + t55) + g(3) * (pkin(5) + rSges(3,2))) - m(4) * (g(1) * (-t38 * rSges(4,1) - t39 * rSges(4,2) + t57) + g(2) * (-t39 * rSges(4,1) + t38 * rSges(4,2) + t53) + g(3) * (-rSges(4,3) + t61)) - m(5) * (g(1) * t56 + g(2) * t52 + g(3) * (-t49 * rSges(5,1) - t51 * rSges(5,2) + t61) + (g(1) * t67 + g(2) * t54) * t39 + (g(1) * t54 - g(2) * t67) * t38) - m(6) * (g(1) * (-t38 * t68 + t39 * pkin(6) + (-t38 * t62 + t39 * t48) * rSges(6,1) + (t38 * t63 + t39 * t50) * rSges(6,2) + t56) + g(2) * (-t39 * t68 - t38 * pkin(6) + (-t38 * t48 - t39 * t62) * rSges(6,1) + (-t38 * t50 + t39 * t63) * rSges(6,2) + t52) + g(3) * (t69 * t51 + t61) + (g(3) * (-rSges(6,1) * t50 + rSges(6,2) * t48 - pkin(4)) - (g(1) * t38 + g(2) * t39) * t69) * t49);
U = t1;
