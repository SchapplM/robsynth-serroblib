% Calculate Gravitation load on the joints for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [3x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S3RRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:53
% EndTime: 2018-11-14 10:15:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (66->24), mult. (50->34), div. (0->0), fcn. (32->6), ass. (0->17)
t11 = sin(qJ(1));
t19 = pkin(1) * t11;
t10 = qJ(1) + qJ(2);
t6 = sin(t10);
t7 = cos(t10);
t18 = t7 * rSges(3,1) - rSges(3,2) * t6;
t8 = qJ(3) + t10;
t4 = sin(t8);
t5 = cos(t8);
t17 = t5 * rSges(4,1) - rSges(4,2) * t4;
t16 = pkin(2) * t7 + t17;
t15 = -rSges(3,1) * t6 - rSges(3,2) * t7;
t14 = -rSges(4,1) * t4 - rSges(4,2) * t5;
t13 = -pkin(2) * t6 + t14;
t12 = cos(qJ(1));
t9 = t12 * pkin(1);
t1 = [-m(2) * (g(1) * (-t11 * rSges(2,1) - rSges(2,2) * t12) + g(2) * (rSges(2,1) * t12 - t11 * rSges(2,2))) - m(3) * (g(1) * (t15 - t19) + g(2) * (t18 + t9)) - m(4) * (g(1) * (t13 - t19) + g(2) * (t16 + t9)) -m(3) * (g(1) * t15 + g(2) * t18) - m(4) * (g(1) * t13 + g(2) * t16) -m(4) * (g(1) * t14 + g(2) * t17)];
taug  = t1(:);
