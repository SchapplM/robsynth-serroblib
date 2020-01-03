% Calculate potential energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:16
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (108->78), mult. (179->99), div. (0->0), fcn. (167->6), ass. (0->26)
t70 = -rSges(5,3) - pkin(7);
t54 = sin(qJ(2));
t69 = t54 * pkin(2) + pkin(5);
t53 = sin(qJ(4));
t68 = t53 * t54;
t55 = sin(qJ(1));
t67 = t54 * t55;
t57 = cos(qJ(2));
t66 = t55 * t57;
t65 = -rSges(6,3) - qJ(5) - pkin(7);
t58 = cos(qJ(1));
t64 = t58 * pkin(1) + t55 * pkin(6);
t63 = qJ(3) * t54;
t50 = t55 * pkin(1);
t62 = pkin(2) * t66 + t55 * t63 + t50;
t61 = t64 + (pkin(2) * t57 + t63) * t58;
t56 = cos(qJ(4));
t42 = -t53 * t57 + t54 * t56;
t60 = t56 * t57 + t68;
t47 = pkin(4) * t56 + pkin(3);
t59 = pkin(4) * t68 + t47 * t57;
t40 = t60 * t58;
t39 = t42 * t58;
t38 = t60 * t55;
t37 = t42 * t55;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t58 - t55 * rSges(2,2)) + g(2) * (t55 * rSges(2,1) + rSges(2,2) * t58) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t55 * rSges(3,3) + t64) + g(2) * (rSges(3,1) * t66 - rSges(3,2) * t67 + t50) + g(3) * (rSges(3,1) * t54 + rSges(3,2) * t57 + pkin(5)) + (g(1) * (rSges(3,1) * t57 - rSges(3,2) * t54) + g(2) * (-rSges(3,3) - pkin(6))) * t58) - m(4) * (g(1) * (t55 * rSges(4,2) + t61) + g(2) * (rSges(4,1) * t66 + rSges(4,3) * t67 + t62) + g(3) * (rSges(4,1) * t54 + (-rSges(4,3) - qJ(3)) * t57 + t69) + (g(1) * (rSges(4,1) * t57 + rSges(4,3) * t54) + g(2) * (-rSges(4,2) - pkin(6))) * t58) - m(5) * (g(1) * (t40 * rSges(5,1) + t39 * rSges(5,2) + t70 * t55 + t61) + g(2) * (t38 * rSges(5,1) + t37 * rSges(5,2) + pkin(3) * t66 + t62) + g(3) * (rSges(5,1) * t42 - rSges(5,2) * t60 + pkin(3) * t54 - qJ(3) * t57 + t69) + (g(1) * pkin(3) * t57 + g(2) * (-pkin(6) - t70)) * t58) - m(6) * (g(1) * (t40 * rSges(6,1) + t39 * rSges(6,2) + t61) + g(2) * (t38 * rSges(6,1) + t37 * rSges(6,2) + t62) + g(3) * (rSges(6,1) * t42 - rSges(6,2) * t60 + t47 * t54 + (-pkin(4) * t53 - qJ(3)) * t57 + t69) + (g(1) * t65 + g(2) * t59) * t55 + (g(1) * t59 + g(2) * (-pkin(6) - t65)) * t58);
U = t1;
