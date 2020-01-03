% Calculate potential energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:41
% EndTime: 2019-12-31 16:39:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (65->44), mult. (92->52), div. (0->0), fcn. (86->6), ass. (0->16)
t46 = rSges(5,3) + pkin(5);
t45 = sin(qJ(1));
t44 = pkin(4) - qJ(3);
t37 = cos(qJ(1));
t43 = t37 * pkin(1) + t45 * qJ(2);
t42 = cos(pkin(6));
t41 = sin(pkin(6));
t40 = t37 * pkin(2) + t43;
t32 = t45 * pkin(1);
t39 = t45 * pkin(2) - t37 * qJ(2) + t32;
t35 = sin(qJ(4));
t36 = cos(qJ(4));
t38 = -rSges(5,1) * t36 + rSges(5,2) * t35 - pkin(3);
t26 = t37 * t41 - t45 * t42;
t25 = -t37 * t42 - t45 * t41;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t37 * rSges(2,1) - t45 * rSges(2,2)) + g(2) * (t45 * rSges(2,1) + t37 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t37 * rSges(3,1) + t45 * rSges(3,3) + t43) + g(2) * (t45 * rSges(3,1) + t32 + (-rSges(3,3) - qJ(2)) * t37) + g(3) * (pkin(4) + rSges(3,2))) - m(4) * (g(1) * (-t25 * rSges(4,1) - t26 * rSges(4,2) + t40) + g(2) * (-t26 * rSges(4,1) + t25 * rSges(4,2) + t39) + g(3) * (-rSges(4,3) + t44)) - m(5) * (g(1) * t40 + g(2) * t39 + g(3) * (-t35 * rSges(5,1) - t36 * rSges(5,2) + t44) + (g(1) * t46 + g(2) * t38) * t26 + (g(1) * t38 - g(2) * t46) * t25);
U = t1;
