% Calculate potential energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:29
% EndTime: 2019-12-31 18:44:30
% DurationCPUTime: 0.27s
% Computational Cost: add. (155->71), mult. (154->89), div. (0->0), fcn. (142->8), ass. (0->29)
t77 = rSges(6,1) + pkin(4);
t76 = rSges(6,3) + qJ(5);
t58 = qJ(1) + pkin(8);
t53 = sin(t58);
t60 = sin(qJ(3));
t75 = t53 * t60;
t63 = cos(qJ(3));
t74 = t53 * t63;
t54 = cos(t58);
t73 = t54 * t60;
t59 = sin(qJ(4));
t72 = t59 * t63;
t62 = cos(qJ(4));
t71 = t62 * t63;
t70 = pkin(5) + qJ(2);
t61 = sin(qJ(1));
t56 = t61 * pkin(1);
t69 = t53 * pkin(2) + t56;
t68 = t60 * pkin(3) + t70;
t64 = cos(qJ(1));
t57 = t64 * pkin(1);
t67 = t54 * pkin(2) + t53 * pkin(6) + t57;
t66 = t54 * t63 * pkin(3) + pkin(7) * t73 + t67;
t65 = pkin(3) * t74 - t54 * pkin(6) + pkin(7) * t75 + t69;
t44 = t53 * t59 + t54 * t71;
t43 = -t53 * t62 + t54 * t72;
t42 = t53 * t71 - t54 * t59;
t41 = t53 * t72 + t54 * t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t64 * rSges(2,1) - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + t64 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t54 * rSges(3,1) - t53 * rSges(3,2) + t57) + g(2) * (t53 * rSges(3,1) + t54 * rSges(3,2) + t56) + g(3) * (rSges(3,3) + t70)) - m(4) * (g(1) * (t53 * rSges(4,3) + t67) + g(2) * (rSges(4,1) * t74 - rSges(4,2) * t75 + t69) + g(3) * (t60 * rSges(4,1) + t63 * rSges(4,2) + t70) + (g(1) * (rSges(4,1) * t63 - rSges(4,2) * t60) + g(2) * (-rSges(4,3) - pkin(6))) * t54) - m(5) * (g(1) * (t44 * rSges(5,1) - t43 * rSges(5,2) + rSges(5,3) * t73 + t66) + g(2) * (t42 * rSges(5,1) - t41 * rSges(5,2) + rSges(5,3) * t75 + t65) + g(3) * ((-rSges(5,3) - pkin(7)) * t63 + (rSges(5,1) * t62 - rSges(5,2) * t59) * t60 + t68)) - m(6) * (g(1) * (t76 * t43 + t77 * t44 + t66) + g(2) * (t76 * t41 + t77 * t42 + t65) + g(3) * (t68 + (-rSges(6,2) - pkin(7)) * t63) + (g(3) * (t76 * t59 + t77 * t62) + (g(1) * t54 + g(2) * t53) * rSges(6,2)) * t60);
U = t1;
