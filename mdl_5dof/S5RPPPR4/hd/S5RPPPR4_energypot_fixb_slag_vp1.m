% Calculate potential energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:04
% EndTime: 2019-12-31 17:45:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (117->57), mult. (89->61), div. (0->0), fcn. (65->8), ass. (0->22)
t53 = rSges(6,3) + pkin(6) + qJ(4);
t52 = pkin(5) + qJ(2);
t39 = qJ(1) + pkin(7);
t33 = sin(t39);
t43 = sin(qJ(1));
t36 = t43 * pkin(1);
t51 = t33 * pkin(2) + t36;
t50 = rSges(5,3) + qJ(4);
t49 = pkin(3) + t52;
t35 = cos(t39);
t44 = cos(qJ(1));
t37 = t44 * pkin(1);
t48 = t35 * pkin(2) + t33 * qJ(3) + t37;
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t47 = rSges(5,1) * t40 + rSges(5,2) * t41;
t38 = pkin(8) + qJ(5);
t32 = sin(t38);
t34 = cos(t38);
t46 = rSges(6,1) * t32 + rSges(6,2) * t34 + pkin(4) * t40;
t45 = g(1) * t48 + g(2) * t51;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t44 * rSges(2,1) - t43 * rSges(2,2)) + g(2) * (t43 * rSges(2,1) + t44 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t35 * rSges(3,1) - t33 * rSges(3,2) + t37) + g(2) * (t33 * rSges(3,1) + t35 * rSges(3,2) + t36) + g(3) * (rSges(3,3) + t52)) - m(4) * (g(1) * (-t35 * rSges(4,2) + t33 * rSges(4,3) + t48) + g(2) * (-t33 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t35 + t51) + g(3) * (rSges(4,1) + t52)) - m(5) * (g(3) * (t41 * rSges(5,1) - t40 * rSges(5,2) + t49) + (g(1) * t47 + g(2) * t50) * t33 + (g(1) * t50 + g(2) * (-qJ(3) - t47)) * t35 + t45) - m(6) * (g(3) * (t34 * rSges(6,1) - t32 * rSges(6,2) + t41 * pkin(4) + t49) + (g(1) * t46 + g(2) * t53) * t33 + (g(1) * t53 + g(2) * (-qJ(3) - t46)) * t35 + t45);
U = t1;
