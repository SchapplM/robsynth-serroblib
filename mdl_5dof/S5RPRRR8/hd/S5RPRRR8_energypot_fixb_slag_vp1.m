% Calculate potential energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (106->58), mult. (141->66), div. (0->0), fcn. (141->8), ass. (0->22)
t62 = pkin(5) - pkin(6);
t61 = rSges(5,3) + pkin(7);
t60 = cos(qJ(3));
t59 = sin(qJ(1));
t58 = sin(qJ(3));
t57 = rSges(6,3) + pkin(8) + pkin(7);
t49 = cos(qJ(1));
t56 = t49 * pkin(1) + t59 * qJ(2);
t55 = t49 * pkin(2) + t56;
t43 = t59 * pkin(1);
t54 = t59 * pkin(2) - t49 * qJ(2) + t43;
t47 = sin(qJ(4));
t48 = cos(qJ(4));
t53 = -rSges(5,1) * t48 + rSges(5,2) * t47 - pkin(3);
t46 = qJ(4) + qJ(5);
t39 = sin(t46);
t40 = cos(t46);
t52 = -rSges(6,1) * t40 + rSges(6,2) * t39 - t48 * pkin(4) - pkin(3);
t51 = g(1) * t55 + g(2) * t54;
t34 = t49 * t58 - t59 * t60;
t33 = -t49 * t60 - t59 * t58;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t49 * rSges(2,1) - t59 * rSges(2,2)) + g(2) * (t59 * rSges(2,1) + t49 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t49 * rSges(3,1) + t59 * rSges(3,3) + t56) + g(2) * (t59 * rSges(3,1) + t43 + (-rSges(3,3) - qJ(2)) * t49) + g(3) * (pkin(5) + rSges(3,2))) - m(4) * (g(1) * (-t33 * rSges(4,1) - t34 * rSges(4,2) + t55) + g(2) * (-t34 * rSges(4,1) + t33 * rSges(4,2) + t54) + g(3) * (-rSges(4,3) + t62)) - m(5) * (g(3) * (-t47 * rSges(5,1) - t48 * rSges(5,2) + t62) + (g(1) * t61 + g(2) * t53) * t34 + (g(1) * t53 - g(2) * t61) * t33 + t51) - m(6) * (g(3) * (-t39 * rSges(6,1) - t40 * rSges(6,2) - t47 * pkin(4) + t62) + (g(1) * t57 + g(2) * t52) * t34 + (g(1) * t52 - g(2) * t57) * t33 + t51);
U = t1;
