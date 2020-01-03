% Calculate potential energy for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:07
% EndTime: 2019-12-31 17:15:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (77->43), mult. (81->50), div. (0->0), fcn. (61->6), ass. (0->18)
t44 = rSges(5,1) + pkin(3);
t36 = -pkin(6) - pkin(5);
t43 = rSges(3,3) + pkin(5);
t32 = sin(qJ(2));
t42 = t32 * pkin(2) + pkin(4);
t34 = cos(qJ(2));
t25 = t34 * pkin(2) + pkin(1);
t41 = rSges(4,3) - t36;
t40 = rSges(5,3) + qJ(4) - t36;
t39 = rSges(3,1) * t34 - rSges(3,2) * t32 + pkin(1);
t31 = qJ(2) + qJ(3);
t26 = sin(t31);
t27 = cos(t31);
t38 = rSges(4,1) * t27 - rSges(4,2) * t26 + t25;
t37 = -rSges(5,2) * t26 + t44 * t27 + t25;
t35 = cos(qJ(1));
t33 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t35 * rSges(2,1) - t33 * rSges(2,2)) + g(2) * (t33 * rSges(2,1) + t35 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(3) * (t32 * rSges(3,1) + t34 * rSges(3,2) + pkin(4)) + (g(1) * t39 - g(2) * t43) * t35 + (g(1) * t43 + g(2) * t39) * t33) - m(4) * (g(3) * (t26 * rSges(4,1) + t27 * rSges(4,2) + t42) + (g(1) * t38 - g(2) * t41) * t35 + (g(1) * t41 + g(2) * t38) * t33) - m(5) * (g(3) * (t27 * rSges(5,2) + t44 * t26 + t42) + (g(1) * t37 - g(2) * t40) * t35 + (g(1) * t40 + g(2) * t37) * t33);
U = t1;
