% Calculate potential energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:14
% EndTime: 2019-12-31 17:32:14
% DurationCPUTime: 0.22s
% Computational Cost: add. (106->58), mult. (141->66), div. (0->0), fcn. (141->8), ass. (0->22)
t62 = cos(qJ(3));
t61 = sin(qJ(3));
t60 = rSges(6,3) + pkin(6) + qJ(4);
t59 = qJ(1) - pkin(5);
t49 = cos(pkin(7));
t56 = sin(pkin(7));
t58 = t49 * pkin(1) + t56 * qJ(2);
t57 = rSges(5,3) + qJ(4);
t55 = t49 * pkin(2) + t58;
t43 = t56 * pkin(1);
t54 = t56 * pkin(2) - t49 * qJ(2) + t43;
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t53 = -rSges(5,1) * t48 + rSges(5,2) * t47 - pkin(3);
t46 = pkin(8) + qJ(5);
t39 = sin(t46);
t40 = cos(t46);
t52 = -rSges(6,1) * t40 + rSges(6,2) * t39 - pkin(4) * t48 - pkin(3);
t51 = g(1) * t55 + g(2) * t54;
t34 = t49 * t61 - t56 * t62;
t33 = -t49 * t62 - t56 * t61;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t49 * rSges(2,1) - t56 * rSges(2,2)) + g(2) * (t56 * rSges(2,1) + t49 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t49 * rSges(3,1) + t56 * rSges(3,3) + t58) + g(2) * (t56 * rSges(3,1) + t43 + (-rSges(3,3) - qJ(2)) * t49) + g(3) * (qJ(1) + rSges(3,2))) - m(4) * (g(1) * (-rSges(4,1) * t33 - rSges(4,2) * t34 + t55) + g(2) * (-rSges(4,1) * t34 + rSges(4,2) * t33 + t54) + g(3) * (-rSges(4,3) + t59)) - m(5) * (g(3) * (-rSges(5,1) * t47 - rSges(5,2) * t48 + t59) + (g(1) * t57 + g(2) * t53) * t34 + (g(1) * t53 - g(2) * t57) * t33 + t51) - m(6) * (g(3) * (-rSges(6,1) * t39 - rSges(6,2) * t40 - pkin(4) * t47 + t59) + (g(1) * t60 + g(2) * t52) * t34 + (g(1) * t52 - g(2) * t60) * t33 + t51);
U = t1;
