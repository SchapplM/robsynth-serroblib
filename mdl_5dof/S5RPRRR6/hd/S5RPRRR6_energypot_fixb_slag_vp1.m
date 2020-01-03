% Calculate potential energy for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:53
% EndTime: 2019-12-31 19:00:53
% DurationCPUTime: 0.29s
% Computational Cost: add. (149->70), mult. (117->86), div. (0->0), fcn. (97->10), ass. (0->28)
t69 = rSges(6,3) + pkin(8);
t68 = rSges(4,3) + pkin(6);
t50 = qJ(3) + qJ(4);
t44 = sin(t50);
t66 = rSges(5,2) * t44;
t49 = qJ(1) + pkin(9);
t43 = cos(t49);
t45 = cos(t50);
t65 = t43 * t45;
t51 = sin(qJ(5));
t64 = t45 * t51;
t54 = cos(qJ(5));
t63 = t45 * t54;
t62 = pkin(5) + qJ(2);
t55 = cos(qJ(3));
t41 = t55 * pkin(3) + pkin(2);
t56 = cos(qJ(1));
t48 = t56 * pkin(1);
t61 = t43 * t41 + t48;
t52 = sin(qJ(3));
t60 = t52 * pkin(3) + t62;
t42 = sin(t49);
t53 = sin(qJ(1));
t47 = t53 * pkin(1);
t57 = -pkin(7) - pkin(6);
t59 = t42 * t41 + t43 * t57 + t47;
t58 = rSges(4,1) * t55 - rSges(4,2) * t52 + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t56 * rSges(2,1) - t53 * rSges(2,2)) + g(2) * (t53 * rSges(2,1) + t56 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t43 * rSges(3,1) - t42 * rSges(3,2) + t48) + g(2) * (t42 * rSges(3,1) + t43 * rSges(3,2) + t47) + g(3) * (rSges(3,3) + t62)) - m(4) * (g(1) * t48 + g(2) * t47 + g(3) * (t52 * rSges(4,1) + t55 * rSges(4,2) + t62) + (g(1) * t58 - g(2) * t68) * t43 + (g(1) * t68 + g(2) * t58) * t42) - m(5) * (g(1) * (rSges(5,1) * t65 - t43 * t66 + t61) + g(2) * (-t43 * rSges(5,3) + t59) + g(3) * (t44 * rSges(5,1) + t45 * rSges(5,2) + t60) + (g(1) * (rSges(5,3) - t57) + g(2) * (rSges(5,1) * t45 - t66)) * t42) - m(6) * (g(1) * (pkin(4) * t65 - t42 * t57 + (t42 * t51 + t43 * t63) * rSges(6,1) + (t42 * t54 - t43 * t64) * rSges(6,2) + t61) + g(2) * (t42 * t45 * pkin(4) + (t42 * t63 - t43 * t51) * rSges(6,1) + (-t42 * t64 - t43 * t54) * rSges(6,2) + t59) + g(3) * (-t69 * t45 + t60) + (g(3) * (rSges(6,1) * t54 - rSges(6,2) * t51 + pkin(4)) + (g(1) * t43 + g(2) * t42) * t69) * t44);
U = t1;
