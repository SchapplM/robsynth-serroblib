% Calculate potential energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:18
% EndTime: 2019-12-31 17:09:18
% DurationCPUTime: 0.28s
% Computational Cost: add. (87->65), mult. (125->85), div. (0->0), fcn. (113->8), ass. (0->22)
t58 = rSges(5,3) + pkin(6) + qJ(3);
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t57 = g(1) * t44 + g(2) * t42;
t56 = rSges(4,3) + qJ(3);
t41 = sin(qJ(2));
t53 = rSges(3,2) * t41;
t38 = sin(pkin(7));
t52 = t42 * t38;
t43 = cos(qJ(2));
t51 = t42 * t43;
t50 = t44 * t38;
t49 = t44 * t43;
t47 = t44 * pkin(1) + t42 * pkin(5);
t35 = t42 * pkin(1);
t45 = -t44 * pkin(5) + t35;
t39 = cos(pkin(7));
t37 = pkin(7) + qJ(4);
t33 = cos(t37);
t32 = sin(t37);
t31 = t39 * pkin(3) + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t44 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t44 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t42 * rSges(3,3) + t47) + g(2) * (rSges(3,1) * t51 - t42 * t53 + t35) + g(3) * (t41 * rSges(3,1) + t43 * rSges(3,2) + pkin(4)) + (g(1) * (rSges(3,1) * t43 - t53) + g(2) * (-rSges(3,3) - pkin(5))) * t44) - m(4) * (g(1) * (pkin(2) * t49 + (t39 * t49 + t52) * rSges(4,1) + (-t38 * t49 + t42 * t39) * rSges(4,2) + t47) + g(2) * (pkin(2) * t51 + (t39 * t51 - t50) * rSges(4,1) + (-t38 * t51 - t44 * t39) * rSges(4,2) + t45) + g(3) * (-t56 * t43 + pkin(4)) + (g(3) * (rSges(4,1) * t39 - rSges(4,2) * t38 + pkin(2)) + t57 * t56) * t41) - m(5) * (g(1) * (t31 * t49 + pkin(3) * t52 + (t42 * t32 + t33 * t49) * rSges(5,1) + (-t32 * t49 + t42 * t33) * rSges(5,2) + t47) + g(2) * (t31 * t51 - pkin(3) * t50 + (-t44 * t32 + t33 * t51) * rSges(5,1) + (-t32 * t51 - t44 * t33) * rSges(5,2) + t45) + g(3) * (-t58 * t43 + pkin(4)) + (g(3) * (rSges(5,1) * t33 - rSges(5,2) * t32 + t31) + t57 * t58) * t41);
U = t1;
