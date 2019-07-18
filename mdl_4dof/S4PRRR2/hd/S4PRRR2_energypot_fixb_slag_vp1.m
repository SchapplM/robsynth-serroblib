% Calculate potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->32), mult. (38->36), div. (0->0), fcn. (18->6), ass. (0->11)
t20 = sin(qJ(2));
t22 = t20 * pkin(1);
t19 = qJ(2) + qJ(3);
t21 = cos(qJ(2));
t18 = t21 * pkin(1);
t17 = qJ(4) + t19;
t16 = cos(t19);
t15 = sin(t19);
t14 = cos(t17);
t13 = sin(t17);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (-g(1) * rSges(2,2) + g(2) * (-qJ(1) - rSges(2,3)) + g(3) * rSges(2,1)) - m(3) * (g(1) * (-t20 * rSges(3,1) - t21 * rSges(3,2)) + g(2) * (-qJ(1) - rSges(3,3)) + g(3) * (t21 * rSges(3,1) - t20 * rSges(3,2))) - m(4) * (g(1) * (-t15 * rSges(4,1) - t16 * rSges(4,2) - t22) + g(2) * (-qJ(1) - rSges(4,3)) + g(3) * (t16 * rSges(4,1) - t15 * rSges(4,2) + t18)) - m(5) * (g(1) * (-t13 * rSges(5,1) - t14 * rSges(5,2) - pkin(2) * t15 - t22) + g(2) * (-qJ(1) - rSges(5,3)) + g(3) * (t14 * rSges(5,1) - t13 * rSges(5,2) + pkin(2) * t16 + t18));
U  = t1;
