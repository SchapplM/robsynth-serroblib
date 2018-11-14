% Calculate Gravitation load on the joints for
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:18
% EndTime: 2018-11-14 14:12:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (49->23), mult. (48->31), div. (0->0), fcn. (28->4), ass. (0->13)
t15 = rSges(5,1) + pkin(3);
t8 = sin(qJ(2));
t14 = pkin(2) * t8;
t7 = qJ(2) + qJ(3);
t4 = sin(t7);
t5 = cos(t7);
t13 = t5 * rSges(4,1) - rSges(4,2) * t4;
t12 = -rSges(5,2) * t4 + t15 * t5;
t11 = -rSges(4,1) * t4 - rSges(4,2) * t5;
t10 = -rSges(5,2) * t5 - t15 * t4;
t9 = cos(qJ(2));
t6 = t9 * pkin(2);
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(1), -m(3) * (g(1) * (rSges(3,1) * t9 - t8 * rSges(3,2)) + g(2) * (-t8 * rSges(3,1) - t9 * rSges(3,2))) - m(4) * (g(1) * (t13 + t6) + g(2) * (t11 - t14)) - m(5) * (g(1) * (t12 + t6) + g(2) * (t10 - t14)) -m(4) * (g(1) * t13 + g(2) * t11) - m(5) * (g(1) * t12 + g(2) * t10) m(5) * g(3)];
taug  = t1(:);
