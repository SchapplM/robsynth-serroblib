% Calculate Gravitation load on the joints for
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:43
% EndTime: 2019-03-08 18:19:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (27->23), mult. (42->30), div. (0->0), fcn. (26->2), ass. (0->7)
t9 = rSges(5,1) + pkin(3);
t8 = -m(4) - m(5);
t5 = sin(qJ(2));
t6 = cos(qJ(2));
t7 = t6 * pkin(2) + t5 * qJ(3);
t3 = t6 * qJ(3);
t1 = [(-m(2) - m(3) + t8) * g(2), -m(3) * (g(1) * (-t5 * rSges(3,1) - rSges(3,2) * t6) + g(2) * (rSges(3,1) * t6 - t5 * rSges(3,2))) - m(4) * (g(1) * (rSges(4,3) * t6 + t3 + (-rSges(4,1) - pkin(2)) * t5) + g(2) * (rSges(4,1) * t6 + t5 * rSges(4,3) + t7)) - m(5) * (g(1) * (rSges(5,2) * t6 + t3) + g(2) * (t6 * t9 + t7) + (g(1) * (-pkin(2) - t9) + g(2) * rSges(5,2)) * t5) t8 * (g(1) * t5 - g(2) * t6) m(5) * g(3)];
taug  = t1(:);
