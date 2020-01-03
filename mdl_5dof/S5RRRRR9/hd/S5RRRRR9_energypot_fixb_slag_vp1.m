% Calculate potential energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:38
% EndTime: 2019-12-31 22:27:39
% DurationCPUTime: 0.41s
% Computational Cost: add. (152->90), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->29)
t76 = rSges(4,3) + pkin(7);
t60 = -pkin(8) - pkin(7);
t75 = rSges(5,3) - t60;
t74 = rSges(6,3) + pkin(9) - t60;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t73 = g(1) * t59 + g(2) * t56;
t57 = cos(qJ(3));
t44 = t57 * pkin(3) + pkin(2);
t55 = sin(qJ(2));
t69 = rSges(3,2) * t55;
t54 = sin(qJ(3));
t68 = t56 * t54;
t58 = cos(qJ(2));
t67 = t56 * t58;
t66 = t59 * t54;
t65 = t59 * t58;
t62 = t59 * pkin(1) + t56 * pkin(6);
t53 = qJ(3) + qJ(4);
t49 = t56 * pkin(1);
t61 = -t59 * pkin(6) + t49;
t47 = qJ(5) + t53;
t46 = cos(t53);
t45 = sin(t53);
t43 = cos(t47);
t42 = sin(t47);
t41 = t54 * pkin(3) + pkin(4) * t45;
t40 = pkin(4) * t46 + t44;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t59 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t59 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t56 * rSges(3,3) + t62) + g(2) * (rSges(3,1) * t67 - t56 * t69 + t49) + g(3) * (t55 * rSges(3,1) + t58 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t58 - t69) + g(2) * (-rSges(3,3) - pkin(6))) * t59) - m(4) * (g(1) * (pkin(2) * t65 + (t57 * t65 + t68) * rSges(4,1) + (-t54 * t65 + t56 * t57) * rSges(4,2) + t62) + g(2) * (pkin(2) * t67 + (t57 * t67 - t66) * rSges(4,1) + (-t54 * t67 - t59 * t57) * rSges(4,2) + t61) + g(3) * (-t76 * t58 + pkin(5)) + (g(3) * (rSges(4,1) * t57 - rSges(4,2) * t54 + pkin(2)) + t73 * t76) * t55) - m(5) * (g(1) * (t44 * t65 + pkin(3) * t68 + (t56 * t45 + t46 * t65) * rSges(5,1) + (-t45 * t65 + t56 * t46) * rSges(5,2) + t62) + g(2) * (t44 * t67 - pkin(3) * t66 + (-t59 * t45 + t46 * t67) * rSges(5,1) + (-t45 * t67 - t59 * t46) * rSges(5,2) + t61) + g(3) * (-t75 * t58 + pkin(5)) + (g(3) * (rSges(5,1) * t46 - rSges(5,2) * t45 + t44) + t73 * t75) * t55) - m(6) * (g(1) * (t40 * t65 + t56 * t41 + (t56 * t42 + t43 * t65) * rSges(6,1) + (-t42 * t65 + t56 * t43) * rSges(6,2) + t62) + g(2) * (t40 * t67 - t59 * t41 + (-t59 * t42 + t43 * t67) * rSges(6,1) + (-t42 * t67 - t59 * t43) * rSges(6,2) + t61) + g(3) * (-t74 * t58 + pkin(5)) + (g(3) * (rSges(6,1) * t43 - rSges(6,2) * t42 + t40) + t73 * t74) * t55);
U = t1;
