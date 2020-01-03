% Calculate potential energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (95->56), mult. (132->62), div. (0->0), fcn. (132->6), ass. (0->19)
t56 = -rSges(6,3) - pkin(6);
t55 = cos(qJ(3));
t54 = sin(qJ(3));
t53 = qJ(1) - pkin(5);
t42 = cos(pkin(7));
t50 = sin(pkin(7));
t52 = t42 * pkin(1) + t50 * qJ(2);
t51 = rSges(5,3) + qJ(4);
t49 = t42 * pkin(2) + t52;
t32 = -t42 * t55 - t50 * t54;
t48 = -t32 * pkin(3) + t49;
t39 = t50 * pkin(1);
t47 = t50 * pkin(2) - t42 * qJ(2) + t39;
t33 = t42 * t54 - t50 * t55;
t46 = -t33 * pkin(3) + t47;
t43 = sin(qJ(5));
t44 = cos(qJ(5));
t45 = rSges(6,1) * t43 + rSges(6,2) * t44 + qJ(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t42 * rSges(2,1) - rSges(2,2) * t50) + g(2) * (rSges(2,1) * t50 + t42 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t42 * rSges(3,1) + rSges(3,3) * t50 + t52) + g(2) * (t50 * rSges(3,1) + t39 + (-rSges(3,3) - qJ(2)) * t42) + g(3) * (qJ(1) + rSges(3,2))) - m(4) * (g(1) * (-t32 * rSges(4,1) - t33 * rSges(4,2) + t49) + g(2) * (-t33 * rSges(4,1) + t32 * rSges(4,2) + t47) + g(3) * (-rSges(4,3) + t53)) - m(5) * (g(1) * (t32 * rSges(5,2) + t51 * t33 + t48) + g(2) * (t33 * rSges(5,2) - t51 * t32 + t46) + g(3) * (-rSges(5,1) + t53)) - m(6) * (g(1) * t48 + g(2) * t46 + g(3) * (-t44 * rSges(6,1) + t43 * rSges(6,2) - pkin(4) + t53) + (g(1) * t45 + g(2) * t56) * t33 + (g(1) * t56 - g(2) * t45) * t32);
U = t1;
