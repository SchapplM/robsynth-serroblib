% Calculate potential energy for
% S4RRPR8
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:44
% EndTime: 2019-12-31 17:07:45
% DurationCPUTime: 0.20s
% Computational Cost: add. (70->56), mult. (115->72), div. (0->0), fcn. (101->6), ass. (0->19)
t37 = sin(qJ(2));
t52 = t37 * pkin(2) + pkin(4);
t51 = -pkin(6) - rSges(5,3);
t38 = sin(qJ(1));
t50 = t37 * t38;
t40 = cos(qJ(2));
t49 = t38 * t40;
t41 = cos(qJ(1));
t48 = t41 * pkin(1) + t38 * pkin(5);
t47 = qJ(3) * t37;
t34 = t38 * pkin(1);
t46 = pkin(2) * t49 + t38 * t47 + t34;
t45 = t48 + (pkin(2) * t40 + t47) * t41;
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t44 = -t36 * t40 + t37 * t39;
t43 = t36 * t37 + t39 * t40;
t42 = t43 * rSges(5,1) + t44 * rSges(5,2) + t40 * pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t41 - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) + rSges(2,2) * t41) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t38 * rSges(3,3) + t48) + g(2) * (rSges(3,1) * t49 - rSges(3,2) * t50 + t34) + g(3) * (rSges(3,1) * t37 + rSges(3,2) * t40 + pkin(4)) + (g(1) * (rSges(3,1) * t40 - rSges(3,2) * t37) + g(2) * (-rSges(3,3) - pkin(5))) * t41) - m(4) * (g(1) * (t38 * rSges(4,2) + t45) + g(2) * (rSges(4,1) * t49 + rSges(4,3) * t50 + t46) + g(3) * (rSges(4,1) * t37 + (-rSges(4,3) - qJ(3)) * t40 + t52) + (g(1) * (rSges(4,1) * t40 + rSges(4,3) * t37) + g(2) * (-rSges(4,2) - pkin(5))) * t41) - m(5) * (g(1) * t45 + g(2) * t46 + g(3) * (t44 * rSges(5,1) - t43 * rSges(5,2) + t37 * pkin(3) - t40 * qJ(3) + t52) + (g(1) * t51 + g(2) * t42) * t38 + (g(1) * t42 + g(2) * (-pkin(5) - t51)) * t41);
U = t1;
