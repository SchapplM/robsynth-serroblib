% Calculate Gravitation load on the joints for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:22
% EndTime: 2018-11-14 13:42:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (77->29), mult. (64->37), div. (0->0), fcn. (54->4), ass. (0->13)
t18 = -m(4) - m(5);
t12 = pkin(6) + qJ(2);
t10 = sin(t12);
t11 = cos(t12);
t17 = t11 * pkin(2) + t10 * qJ(3);
t16 = cos(qJ(4));
t15 = sin(qJ(4));
t1 = -t10 * t15 - t11 * t16;
t2 = -t10 * t16 + t11 * t15;
t14 = -rSges(5,1) * t2 + rSges(5,2) * t1;
t13 = rSges(5,1) * t1 + rSges(5,2) * t2;
t8 = t11 * qJ(3);
t3 = [(-m(2) - m(3) + t18) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t10 - rSges(3,2) * t11) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10)) - m(4) * (g(1) * (rSges(4,3) * t11 + t8 + (-rSges(4,1) - pkin(2)) * t10) + g(2) * (rSges(4,1) * t11 + rSges(4,3) * t10 + t17)) - m(5) * (g(1) * (t8 + (-pkin(2) - pkin(3)) * t10 - t14) + g(2) * (pkin(3) * t11 - t13 + t17)) t18 * (g(1) * t10 - g(2) * t11) -m(5) * (g(1) * t14 + g(2) * t13)];
taug  = t3(:);
