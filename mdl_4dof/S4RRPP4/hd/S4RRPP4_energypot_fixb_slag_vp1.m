% Calculate potential energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:50
% EndTime: 2019-12-31 16:58:50
% DurationCPUTime: 0.19s
% Computational Cost: add. (64->53), mult. (99->65), div. (0->0), fcn. (79->4), ass. (0->16)
t46 = rSges(5,1) + pkin(3);
t33 = sin(qJ(2));
t45 = t33 * pkin(2) + pkin(4);
t34 = sin(qJ(1));
t44 = t33 * t34;
t35 = cos(qJ(2));
t43 = t34 * t35;
t36 = cos(qJ(1));
t42 = t36 * pkin(1) + t34 * pkin(5);
t41 = qJ(3) * t33;
t40 = -rSges(5,3) - qJ(4);
t31 = t34 * pkin(1);
t39 = pkin(2) * t43 + t34 * t41 + t31;
t38 = t42 + (pkin(2) * t35 + t41) * t36;
t37 = rSges(5,2) * t33 + t46 * t35;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t36 - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) + rSges(2,2) * t36) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t34 * rSges(3,3) + t42) + g(2) * (rSges(3,1) * t43 - rSges(3,2) * t44 + t31) + g(3) * (rSges(3,1) * t33 + rSges(3,2) * t35 + pkin(4)) + (g(1) * (rSges(3,1) * t35 - rSges(3,2) * t33) + g(2) * (-rSges(3,3) - pkin(5))) * t36) - m(4) * (g(1) * (t34 * rSges(4,2) + t38) + g(2) * (rSges(4,1) * t43 + rSges(4,3) * t44 + t39) + g(3) * (rSges(4,1) * t33 + (-rSges(4,3) - qJ(3)) * t35 + t45) + (g(1) * (rSges(4,1) * t35 + rSges(4,3) * t33) + g(2) * (-rSges(4,2) - pkin(5))) * t36) - m(5) * (g(1) * t38 + g(2) * t39 + g(3) * ((-rSges(5,2) - qJ(3)) * t35 + t46 * t33 + t45) + (g(1) * t40 + g(2) * t37) * t34 + (g(1) * t37 + g(2) * (-pkin(5) - t40)) * t36);
U = t1;
