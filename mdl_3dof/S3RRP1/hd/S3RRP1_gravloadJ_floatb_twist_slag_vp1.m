% Calculate Gravitation load on the joints for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:51
% EndTime: 2019-03-08 18:06:51
% DurationCPUTime: 0.11s
% Computational Cost: add. (58->23), mult. (50->32), div. (0->0), fcn. (34->4), ass. (0->14)
t20 = rSges(4,1) + pkin(2);
t19 = qJ(3) + rSges(4,3);
t12 = sin(qJ(1));
t18 = t12 * pkin(1);
t11 = qJ(1) + qJ(2);
t8 = sin(t11);
t9 = cos(t11);
t17 = t19 * t8 + t20 * t9;
t16 = t9 * rSges(3,1) - t8 * rSges(3,2);
t15 = -t8 * rSges(3,1) - t9 * rSges(3,2);
t14 = t19 * t9 - t20 * t8;
t13 = cos(qJ(1));
t10 = t13 * pkin(1);
t1 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - t13 * rSges(2,2)) + g(2) * (t13 * rSges(2,1) - t12 * rSges(2,2))) - m(3) * (g(1) * (t15 - t18) + g(2) * (t10 + t16)) - m(4) * (g(1) * (t14 - t18) + g(2) * (t10 + t17)) -m(3) * (g(1) * t15 + g(2) * t16) - m(4) * (g(1) * t14 + g(2) * t17) -m(4) * (g(1) * t8 - g(2) * t9)];
taug  = t1(:);
