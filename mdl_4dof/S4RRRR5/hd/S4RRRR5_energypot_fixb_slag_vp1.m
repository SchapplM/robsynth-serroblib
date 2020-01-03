% Calculate potential energy for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:19
% EndTime: 2019-12-31 17:27:19
% DurationCPUTime: 0.29s
% Computational Cost: add. (87->65), mult. (125->85), div. (0->0), fcn. (113->8), ass. (0->22)
t58 = rSges(4,3) + pkin(6);
t57 = rSges(5,3) + pkin(7) + pkin(6);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t56 = g(1) * t43 + g(2) * t40;
t39 = sin(qJ(2));
t52 = rSges(3,2) * t39;
t38 = sin(qJ(3));
t51 = t40 * t38;
t42 = cos(qJ(2));
t50 = t40 * t42;
t49 = t43 * t38;
t48 = t43 * t42;
t46 = t43 * pkin(1) + t40 * pkin(5);
t35 = t40 * pkin(1);
t45 = -t43 * pkin(5) + t35;
t41 = cos(qJ(3));
t37 = qJ(3) + qJ(4);
t33 = cos(t37);
t32 = sin(t37);
t31 = t41 * pkin(3) + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t43 * rSges(2,1) - t40 * rSges(2,2)) + g(2) * (t40 * rSges(2,1) + t43 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t40 * rSges(3,3) + t46) + g(2) * (rSges(3,1) * t50 - t40 * t52 + t35) + g(3) * (t39 * rSges(3,1) + t42 * rSges(3,2) + pkin(4)) + (g(1) * (rSges(3,1) * t42 - t52) + g(2) * (-rSges(3,3) - pkin(5))) * t43) - m(4) * (g(1) * (pkin(2) * t48 + (t41 * t48 + t51) * rSges(4,1) + (-t38 * t48 + t40 * t41) * rSges(4,2) + t46) + g(2) * (pkin(2) * t50 + (t41 * t50 - t49) * rSges(4,1) + (-t38 * t50 - t43 * t41) * rSges(4,2) + t45) + g(3) * (-t58 * t42 + pkin(4)) + (g(3) * (rSges(4,1) * t41 - rSges(4,2) * t38 + pkin(2)) + t56 * t58) * t39) - m(5) * (g(1) * (t31 * t48 + pkin(3) * t51 + (t40 * t32 + t33 * t48) * rSges(5,1) + (-t32 * t48 + t40 * t33) * rSges(5,2) + t46) + g(2) * (t31 * t50 - pkin(3) * t49 + (-t43 * t32 + t33 * t50) * rSges(5,1) + (-t32 * t50 - t43 * t33) * rSges(5,2) + t45) + g(3) * (-t57 * t42 + pkin(4)) + (g(3) * (rSges(5,1) * t33 - rSges(5,2) * t32 + t31) + t56 * t57) * t39);
U = t1;
