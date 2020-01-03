% Calculate potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:04
% EndTime: 2020-01-03 11:25:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (147->74), mult. (141->90), div. (0->0), fcn. (125->8), ass. (0->28)
t66 = rSges(5,3) + pkin(6);
t65 = rSges(6,3) + qJ(5) + pkin(6);
t53 = cos(qJ(1));
t64 = pkin(1) * t53;
t48 = cos(pkin(8));
t63 = pkin(3) * t48;
t50 = sin(qJ(4));
t61 = t48 * t50;
t52 = cos(qJ(4));
t60 = t48 * t52;
t59 = pkin(5) + qJ(2);
t46 = qJ(1) + pkin(7);
t43 = sin(t46);
t51 = sin(qJ(1));
t45 = t51 * pkin(1);
t58 = t43 * pkin(2) + t45;
t57 = -rSges(4,3) - qJ(3);
t56 = -pkin(4) * t50 - qJ(3);
t47 = sin(pkin(8));
t55 = rSges(4,1) * t48 - rSges(4,2) * t47;
t42 = pkin(4) * t52 + pkin(3);
t54 = t42 * t48 + t65 * t47;
t44 = cos(t46);
t40 = -t43 * t50 - t44 * t60;
t39 = -t43 * t52 + t44 * t61;
t38 = t43 * t60 - t44 * t50;
t37 = -t43 * t61 - t44 * t52;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t51 * rSges(2,1) + rSges(2,2) * t53) + g(3) * (-rSges(2,1) * t53 + t51 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t59) + g(2) * (rSges(3,1) * t43 + rSges(3,2) * t44 + t45) + g(3) * (-t44 * rSges(3,1) + t43 * rSges(3,2) - t64)) - m(4) * (g(1) * (rSges(4,1) * t47 + rSges(4,2) * t48 + t59) + g(2) * t58 - g(3) * t64 + (g(2) * t55 + g(3) * t57) * t43 + (g(2) * t57 + g(3) * (-pkin(2) - t55)) * t44) - m(5) * (g(1) * (-t66 * t48 + t59) + g(2) * (rSges(5,1) * t38 + rSges(5,2) * t37 - qJ(3) * t44 + t43 * t63 + t58) + g(3) * (t40 * rSges(5,1) + t39 * rSges(5,2) - t43 * qJ(3) - t64 + (-pkin(2) - t63) * t44) + (g(1) * (rSges(5,1) * t52 - rSges(5,2) * t50 + pkin(3)) + (g(2) * t43 - g(3) * t44) * t66) * t47) - m(6) * (g(2) * (rSges(6,1) * t38 + rSges(6,2) * t37 + t58) + g(3) * (t40 * rSges(6,1) + t39 * rSges(6,2) - t64) + (g(2) * t54 + g(3) * t56) * t43 + (g(2) * t56 + g(3) * (-pkin(2) - t54)) * t44 + (t59 + (rSges(6,1) * t52 - rSges(6,2) * t50 + t42) * t47 - t65 * t48) * g(1));
U = t1;
