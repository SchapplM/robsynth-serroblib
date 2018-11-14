% Calculate potential energy for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PRPP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP3_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:08
% EndTime: 2018-11-14 14:01:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (41->37), mult. (42->36), div. (0->0), fcn. (22->2), ass. (0->6)
t19 = rSges(5,1) + pkin(3);
t15 = sin(qJ(2));
t18 = t15 * pkin(2) + qJ(1);
t16 = cos(qJ(2));
t17 = t16 * pkin(2) + t15 * qJ(3) + pkin(1);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * rSges(2,1) + g(2) * (qJ(1) + rSges(2,3)) - g(3) * rSges(2,2)) - m(3) * (g(1) * (t16 * rSges(3,1) - t15 * rSges(3,2) + pkin(1)) + g(2) * (t15 * rSges(3,1) + t16 * rSges(3,2) + qJ(1)) + g(3) * (pkin(4) + rSges(3,3))) - m(4) * (g(1) * (t16 * rSges(4,1) + t15 * rSges(4,3) + t17) + g(2) * (t15 * rSges(4,1) + (-rSges(4,3) - qJ(3)) * t16 + t18) + g(3) * (pkin(4) + rSges(4,2))) - m(5) * (g(1) * (t15 * rSges(5,2) + t17) + g(2) * (t19 * t15 + t18) + g(3) * (-qJ(4) + pkin(4) - rSges(5,3)) + (g(1) * t19 + g(2) * (-rSges(5,2) - qJ(3))) * t16);
U  = t1;
