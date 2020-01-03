% Calculate potential energy for
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:09
% EndTime: 2020-01-03 11:49:09
% DurationCPUTime: 0.47s
% Computational Cost: add. (142->88), mult. (180->109), div. (0->0), fcn. (168->8), ass. (0->33)
t80 = rSges(4,3) + pkin(6);
t61 = -pkin(7) - pkin(6);
t79 = rSges(5,3) - t61;
t78 = rSges(6,3) + qJ(5) - t61;
t58 = sin(qJ(1));
t60 = cos(qJ(1));
t77 = g(2) * t58 - g(3) * t60;
t57 = sin(qJ(3));
t76 = pkin(3) * t57;
t59 = cos(qJ(3));
t48 = t59 * pkin(3) + pkin(2);
t56 = cos(pkin(8));
t72 = t56 * t58;
t71 = t56 * t60;
t70 = t57 * t58;
t69 = t57 * t60;
t68 = t58 * t59;
t67 = t59 * t60;
t64 = t58 * qJ(2);
t63 = -rSges(3,3) - qJ(2);
t55 = sin(pkin(8));
t62 = rSges(3,1) * t56 - rSges(3,2) * t55;
t54 = qJ(3) + qJ(4);
t51 = t58 * pkin(1);
t50 = cos(t54);
t49 = sin(t54);
t47 = pkin(4) * t49 + t76;
t46 = pkin(4) * t50 + t48;
t45 = -t49 * t58 - t50 * t71;
t44 = t49 * t71 - t50 * t58;
t43 = -t49 * t60 + t50 * t72;
t42 = -t49 * t72 - t50 * t60;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t60) + g(3) * (-rSges(2,1) * t60 + rSges(2,2) * t58)) - m(3) * (g(1) * (rSges(3,1) * t55 + rSges(3,2) * t56 + pkin(5)) + g(2) * t51 + (g(2) * t62 + g(3) * t63) * t58 + (g(2) * t63 + g(3) * (-pkin(1) - t62)) * t60) - m(4) * (g(1) * (-t80 * t56 + pkin(5)) + g(2) * (pkin(2) * t72 + t51 - t60 * qJ(2) + (t56 * t68 - t69) * rSges(4,1) + (-t56 * t70 - t67) * rSges(4,2)) + g(3) * (-pkin(2) * t71 - t60 * pkin(1) - t64 + (-t56 * t67 - t70) * rSges(4,1) + (t56 * t69 - t68) * rSges(4,2)) + (g(1) * (rSges(4,1) * t59 - rSges(4,2) * t57 + pkin(2)) + t77 * t80) * t55) - m(5) * (g(1) * (-t79 * t56 + pkin(5)) + g(2) * (t43 * rSges(5,1) + t42 * rSges(5,2) + t48 * t72 + t51) + g(3) * (t45 * rSges(5,1) + t44 * rSges(5,2) - pkin(3) * t70 - t64) + (g(2) * (-qJ(2) - t76) + g(3) * (-t48 * t56 - pkin(1))) * t60 + (g(1) * (rSges(5,1) * t50 - rSges(5,2) * t49 + t48) + t77 * t79) * t55) - m(6) * (g(1) * (-t78 * t56 + pkin(5)) + g(2) * (rSges(6,1) * t43 + rSges(6,2) * t42 + t46 * t72 + t51) + g(3) * (t45 * rSges(6,1) + t44 * rSges(6,2) - t58 * t47 - t64) + (g(2) * (-qJ(2) - t47) + g(3) * (-t46 * t56 - pkin(1))) * t60 + (g(1) * (rSges(6,1) * t50 - rSges(6,2) * t49 + t46) + t77 * t78) * t55);
U = t1;
