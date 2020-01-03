% Calculate potential energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:06
% EndTime: 2019-12-31 17:11:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (69->60), mult. (112->77), div. (0->0), fcn. (96->6), ass. (0->20)
t55 = rSges(5,3) + pkin(6);
t38 = sin(qJ(2));
t53 = t38 * pkin(2) + pkin(4);
t39 = sin(qJ(1));
t52 = t38 * t39;
t37 = sin(qJ(4));
t51 = t39 * t37;
t40 = cos(qJ(4));
t50 = t39 * t40;
t41 = cos(qJ(2));
t49 = t39 * t41;
t42 = cos(qJ(1));
t48 = t42 * t37;
t47 = t42 * t40;
t46 = t42 * pkin(1) + t39 * pkin(5);
t45 = qJ(3) * t38;
t35 = t39 * pkin(1);
t44 = pkin(2) * t49 + t39 * t45 + t35;
t43 = t46 + (pkin(2) * t41 + t45) * t42;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t42 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 + rSges(2,2) * t42) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t39 * rSges(3,3) + t46) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t52 + t35) + g(3) * (rSges(3,1) * t38 + rSges(3,2) * t41 + pkin(4)) + (g(1) * (rSges(3,1) * t41 - rSges(3,2) * t38) + g(2) * (-rSges(3,3) - pkin(5))) * t42) - m(4) * (g(1) * (t39 * rSges(4,1) + t43) + g(2) * (-rSges(4,2) * t49 + rSges(4,3) * t52 + t44) + g(3) * (-t38 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t41 + t53) + (g(1) * (-rSges(4,2) * t41 + rSges(4,3) * t38) + g(2) * (-rSges(4,1) - pkin(5))) * t42) - m(5) * (g(1) * (t39 * pkin(3) + (t38 * t48 + t50) * rSges(5,1) + (t38 * t47 - t51) * rSges(5,2) + t43) + g(2) * ((t38 * t51 - t47) * rSges(5,1) + (t38 * t50 + t48) * rSges(5,2) + t44 + (-pkin(3) - pkin(5)) * t42) + g(3) * (t55 * t38 + t53) + (g(3) * (-rSges(5,1) * t37 - rSges(5,2) * t40 - qJ(3)) + (g(1) * t42 + g(2) * t39) * t55) * t41);
U = t1;
