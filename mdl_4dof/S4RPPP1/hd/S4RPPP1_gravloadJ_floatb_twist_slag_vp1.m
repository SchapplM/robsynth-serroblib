% Calculate Gravitation load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RPPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (182->54), mult. (220->72), div. (0->0), fcn. (199->10), ass. (0->29)
t38 = rSges(5,2) + qJ(3);
t37 = rSges(5,3) + qJ(4);
t36 = -m(4) - m(5);
t21 = sin(pkin(4));
t23 = sin(qJ(1));
t34 = t21 * t23;
t24 = cos(qJ(1));
t33 = t21 * t24;
t31 = qJ(2) * t21;
t32 = t24 * pkin(1) + t23 * t31;
t30 = rSges(4,3) + qJ(3);
t19 = pkin(4) - pkin(6);
t14 = cos(t19) / 0.2e1;
t18 = pkin(4) + pkin(6);
t16 = cos(t18);
t29 = t14 + t16 / 0.2e1;
t13 = -sin(t19) / 0.2e1;
t15 = sin(t18);
t28 = t15 / 0.2e1 + t13;
t22 = cos(pkin(6));
t7 = t22 * t24 - t23 * t28;
t27 = t7 * pkin(2) + t32;
t26 = -t23 * pkin(1) + t24 * t31;
t5 = t23 * t22 + t24 * t28;
t25 = -t5 * pkin(2) + t26;
t20 = sin(pkin(6));
t6 = t24 * t20 + t23 * t29;
t4 = t20 * t23 - t24 * t29;
t1 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - rSges(2,2) * t24) + g(2) * (rSges(2,1) * t24 - t23 * rSges(2,2))) - m(3) * (g(1) * (-t5 * rSges(3,1) + t4 * rSges(3,2) + rSges(3,3) * t33 + t26) + g(2) * (rSges(3,1) * t7 - rSges(3,2) * t6 + rSges(3,3) * t34 + t32)) - m(4) * (g(1) * (rSges(4,1) * t33 + t5 * rSges(4,2) - t30 * t4 + t25) + g(2) * (rSges(4,1) * t34 - rSges(4,2) * t7 + t30 * t6 + t27)) - m(5) * (g(1) * (-t37 * t5 - t38 * t4 + t25) + g(2) * (t37 * t7 + t38 * t6 + t27) + (g(1) * t24 + g(2) * t23) * t21 * (rSges(5,1) + pkin(3))) (-m(3) + t36) * (g(3) * cos(pkin(4)) + (g(1) * t23 - g(2) * t24) * t21) t36 * (g(1) * t6 + g(2) * t4 + g(3) * (-t15 / 0.2e1 + t13)) -m(5) * (g(1) * t7 + g(2) * t5 + g(3) * (t14 - t16 / 0.2e1))];
taug  = t1(:);
