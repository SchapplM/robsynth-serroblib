% Calculate Gravitation load on the joints for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:04
% EndTime: 2018-11-14 13:39:04
% DurationCPUTime: 0.08s
% Computational Cost: add. (35->16), mult. (62->25), div. (0->0), fcn. (62->4), ass. (0->10)
t13 = rSges(5,1) + pkin(3);
t12 = cos(qJ(3));
t11 = sin(qJ(3));
t10 = -rSges(5,3) - qJ(4);
t9 = sin(pkin(5));
t8 = -m(3) - m(4) - m(5);
t7 = cos(pkin(5));
t2 = t7 * t11 - t9 * t12;
t1 = -t9 * t11 - t7 * t12;
t3 = [(-m(2) + t8) * g(3), t8 * (g(1) * t9 - g(2) * t7) -m(4) * (g(1) * (-t2 * rSges(4,1) + t1 * rSges(4,2)) + g(2) * (t1 * rSges(4,1) + t2 * rSges(4,2))) - m(5) * ((-g(1) * t13 + g(2) * t10) * t2 + (g(1) * t10 + g(2) * t13) * t1) -m(5) * (g(1) * t2 - g(2) * t1)];
taug  = t3(:);
