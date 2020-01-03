% Calculate potential energy for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:35
% EndTime: 2019-12-31 17:16:35
% DurationCPUTime: 0.15s
% Computational Cost: add. (81->47), mult. (88->56), div. (0->0), fcn. (68->6), ass. (0->19)
t48 = rSges(5,1) + pkin(3);
t33 = qJ(2) + qJ(3);
t30 = sin(t33);
t31 = cos(t33);
t47 = rSges(4,1) * t31 - rSges(4,2) * t30;
t46 = rSges(5,3) + qJ(4);
t45 = rSges(3,3) + pkin(5);
t34 = sin(qJ(2));
t44 = t34 * pkin(2) + pkin(4);
t36 = cos(qJ(2));
t28 = t36 * pkin(2) + pkin(1);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t38 = -pkin(6) - pkin(5);
t41 = t35 * t28 + t37 * t38;
t40 = rSges(3,1) * t36 - rSges(3,2) * t34 + pkin(1);
t39 = t46 * t30 + t48 * t31;
t27 = t37 * t28;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t37 * rSges(2,1) - t35 * rSges(2,2)) + g(2) * (t35 * rSges(2,1) + t37 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(3) * (t34 * rSges(3,1) + t36 * rSges(3,2) + pkin(4)) + (g(1) * t40 - g(2) * t45) * t37 + (g(1) * t45 + g(2) * t40) * t35) - m(4) * (g(1) * (t47 * t37 + t27) + g(2) * (-t37 * rSges(4,3) + t41) + g(3) * (t30 * rSges(4,1) + t31 * rSges(4,2) + t44) + (g(1) * (rSges(4,3) - t38) + g(2) * t47) * t35) - m(5) * (g(1) * t27 + g(2) * t41 + g(3) * (t48 * t30 - t46 * t31 + t44) + (-g(2) * rSges(5,2) + g(1) * t39) * t37 + (g(1) * (rSges(5,2) - t38) + g(2) * t39) * t35);
U = t1;
