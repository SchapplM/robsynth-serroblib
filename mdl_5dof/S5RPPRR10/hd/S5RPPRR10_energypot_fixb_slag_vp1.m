% Calculate potential energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:55
% EndTime: 2019-12-31 18:03:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (120->78), mult. (179->96), div. (0->0), fcn. (167->8), ass. (0->27)
t49 = sin(pkin(8));
t70 = t49 * pkin(2) + pkin(5);
t69 = -pkin(6) - rSges(5,3);
t51 = sin(qJ(4));
t68 = t49 * t51;
t52 = sin(qJ(1));
t67 = t49 * t52;
t50 = cos(pkin(8));
t66 = t52 * t50;
t54 = cos(qJ(1));
t65 = t54 * pkin(1) + t52 * qJ(2);
t64 = -pkin(7) - pkin(6) - rSges(6,3);
t63 = qJ(3) * t49;
t46 = t52 * pkin(1);
t62 = pkin(2) * t66 + t52 * t63 + t46;
t61 = t65 + (pkin(2) * t50 + t63) * t54;
t53 = cos(qJ(4));
t60 = t49 * t53 - t50 * t51;
t59 = t50 * t53 + t68;
t58 = g(1) * t61 + g(2) * t62;
t57 = t59 * rSges(5,1) + t60 * rSges(5,2) + t50 * pkin(3);
t41 = pkin(4) * t53 + pkin(3);
t48 = qJ(4) + qJ(5);
t42 = sin(t48);
t43 = cos(t48);
t56 = t50 * t41 + pkin(4) * t68 + (t42 * t49 + t43 * t50) * rSges(6,1) + (-t42 * t50 + t43 * t49) * rSges(6,2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t54 - rSges(2,2) * t52) + g(2) * (rSges(2,1) * t52 + rSges(2,2) * t54) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t52 + t65) + g(2) * (rSges(3,1) * t66 - rSges(3,2) * t67 + t46) + g(3) * (rSges(3,1) * t49 + rSges(3,2) * t50 + pkin(5)) + (g(1) * (rSges(3,1) * t50 - rSges(3,2) * t49) + g(2) * (-rSges(3,3) - qJ(2))) * t54) - m(4) * (g(1) * (rSges(4,2) * t52 + t61) + g(2) * (rSges(4,1) * t66 + rSges(4,3) * t67 + t62) + g(3) * (rSges(4,1) * t49 + (-rSges(4,3) - qJ(3)) * t50 + t70) + (g(1) * (rSges(4,1) * t50 + rSges(4,3) * t49) + g(2) * (-rSges(4,2) - qJ(2))) * t54) - m(5) * (g(3) * (t60 * rSges(5,1) - t59 * rSges(5,2) + t49 * pkin(3) - t50 * qJ(3) + t70) + (g(1) * t69 + g(2) * t57) * t52 + (g(1) * t57 + g(2) * (-qJ(2) - t69)) * t54 + t58) - m(6) * (g(3) * ((t43 * rSges(6,1) - t42 * rSges(6,2) + t41) * t49 + (-t42 * rSges(6,1) - t43 * rSges(6,2) - t51 * pkin(4) - qJ(3)) * t50 + t70) + (g(1) * t64 + g(2) * t56) * t52 + (g(1) * t56 + g(2) * (-qJ(2) - t64)) * t54 + t58);
U = t1;
