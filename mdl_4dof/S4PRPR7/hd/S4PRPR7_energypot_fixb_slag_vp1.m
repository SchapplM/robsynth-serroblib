% Calculate potential energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:34
% EndTime: 2019-12-31 16:25:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (69->60), mult. (112->78), div. (0->0), fcn. (96->6), ass. (0->18)
t58 = rSges(5,3) + pkin(5);
t42 = sin(pkin(6));
t45 = sin(qJ(2));
t56 = t42 * t45;
t47 = cos(qJ(2));
t55 = t42 * t47;
t44 = sin(qJ(4));
t54 = t44 * t45;
t46 = cos(qJ(4));
t53 = t45 * t46;
t43 = cos(pkin(6));
t52 = t43 * pkin(1) + t42 * pkin(4);
t51 = qJ(3) * t45;
t50 = t45 * pkin(2) + qJ(1);
t39 = t42 * pkin(1);
t49 = pkin(2) * t55 + t42 * t51 + t39;
t48 = t52 + (pkin(2) * t47 + t51) * t43;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t43 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t43 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t42 * rSges(3,3) + t52) + g(2) * (rSges(3,1) * t55 - rSges(3,2) * t56 + t39) + g(3) * (t45 * rSges(3,1) + t47 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t47 - rSges(3,2) * t45) + g(2) * (-rSges(3,3) - pkin(4))) * t43) - m(4) * (g(1) * (t42 * rSges(4,1) + t48) + g(2) * (-rSges(4,2) * t55 + rSges(4,3) * t56 + t49) + g(3) * (-t45 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t47 + t50) + (g(1) * (-rSges(4,2) * t47 + rSges(4,3) * t45) + g(2) * (-rSges(4,1) - pkin(4))) * t43) - m(5) * (g(1) * (t42 * pkin(3) + (t42 * t46 + t43 * t54) * rSges(5,1) + (-t42 * t44 + t43 * t53) * rSges(5,2) + t48) + g(2) * (t49 + (t54 * rSges(5,1) + t53 * rSges(5,2)) * t42 + (-t46 * rSges(5,1) + t44 * rSges(5,2) - pkin(3) - pkin(4)) * t43) + g(3) * (t58 * t45 + t50) + (g(3) * (-rSges(5,1) * t44 - rSges(5,2) * t46 - qJ(3)) + (g(1) * t43 + g(2) * t42) * t58) * t47);
U = t1;
