% Calculate potential energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:12
% EndTime: 2018-11-14 14:12:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (49->36), mult. (38->34), div. (0->0), fcn. (18->4), ass. (0->10)
t25 = -pkin(4) - pkin(5);
t24 = rSges(5,1) + pkin(3);
t21 = cos(qJ(2));
t23 = t21 * pkin(2) + pkin(1);
t20 = sin(qJ(2));
t22 = t20 * pkin(2) + qJ(1);
t19 = qJ(2) + qJ(3);
t16 = cos(t19);
t15 = sin(t19);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (qJ(1) + rSges(2,3)) + g(2) * rSges(2,1) + g(3) * rSges(2,2)) - m(3) * (g(1) * (t20 * rSges(3,1) + t21 * rSges(3,2) + qJ(1)) + g(2) * (t21 * rSges(3,1) - t20 * rSges(3,2) + pkin(1)) + g(3) * (-pkin(4) - rSges(3,3))) - m(4) * (g(1) * (t15 * rSges(4,1) + t16 * rSges(4,2) + t22) + g(2) * (t16 * rSges(4,1) - t15 * rSges(4,2) + t23) + g(3) * (-rSges(4,3) + t25)) - m(5) * (g(1) * (t16 * rSges(5,2) + t24 * t15 + t22) + g(2) * (-t15 * rSges(5,2) + t24 * t16 + t23) + g(3) * (-qJ(4) - rSges(5,3) + t25));
U  = t1;
