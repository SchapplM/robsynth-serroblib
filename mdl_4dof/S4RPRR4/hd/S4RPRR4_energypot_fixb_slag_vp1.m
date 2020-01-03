% Calculate potential energy for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:12
% EndTime: 2019-12-31 16:50:12
% DurationCPUTime: 0.21s
% Computational Cost: add. (92->54), mult. (89->70), div. (0->0), fcn. (73->8), ass. (0->20)
t51 = rSges(5,3) + pkin(6);
t38 = sin(qJ(3));
t49 = rSges(4,2) * t38;
t36 = qJ(1) + pkin(7);
t32 = sin(t36);
t41 = cos(qJ(3));
t48 = t32 * t41;
t37 = sin(qJ(4));
t47 = t37 * t41;
t40 = cos(qJ(4));
t46 = t40 * t41;
t45 = pkin(4) + qJ(2);
t39 = sin(qJ(1));
t34 = t39 * pkin(1);
t44 = t32 * pkin(2) + t34;
t33 = cos(t36);
t42 = cos(qJ(1));
t35 = t42 * pkin(1);
t43 = t33 * pkin(2) + t32 * pkin(5) + t35;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t42 * rSges(2,1) - t39 * rSges(2,2)) + g(2) * (t39 * rSges(2,1) + t42 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t33 * rSges(3,1) - t32 * rSges(3,2) + t35) + g(2) * (t32 * rSges(3,1) + t33 * rSges(3,2) + t34) + g(3) * (rSges(3,3) + t45)) - m(4) * (g(1) * (t32 * rSges(4,3) + t43) + g(2) * (rSges(4,1) * t48 - t32 * t49 + t44) + g(3) * (t38 * rSges(4,1) + t41 * rSges(4,2) + t45) + (g(1) * (rSges(4,1) * t41 - t49) + g(2) * (-rSges(4,3) - pkin(5))) * t33) - m(5) * (g(1) * (t33 * t41 * pkin(3) + (t32 * t37 + t33 * t46) * rSges(5,1) + (t32 * t40 - t33 * t47) * rSges(5,2) + t43) + g(2) * (pkin(3) * t48 - t33 * pkin(5) + (t32 * t46 - t33 * t37) * rSges(5,1) + (-t32 * t47 - t33 * t40) * rSges(5,2) + t44) + g(3) * (-t51 * t41 + t45) + (g(3) * (rSges(5,1) * t40 - rSges(5,2) * t37 + pkin(3)) + (g(1) * t33 + g(2) * t32) * t51) * t38);
U = t1;
