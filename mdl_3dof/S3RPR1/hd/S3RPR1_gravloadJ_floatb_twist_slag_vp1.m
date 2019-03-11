% Calculate Gravitation load on the joints for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S3RPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (36->26), mult. (60->36), div. (0->0), fcn. (54->4), ass. (0->11)
t10 = sin(qJ(1));
t11 = cos(qJ(1));
t16 = t11 * pkin(1) + t10 * qJ(2);
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t1 = -t10 * t14 - t11 * t15;
t2 = -t10 * t15 + t11 * t14;
t13 = -rSges(4,1) * t2 + rSges(4,2) * t1;
t12 = t1 * rSges(4,1) + t2 * rSges(4,2);
t8 = t11 * qJ(2);
t3 = [-m(2) * (g(1) * (-t10 * rSges(2,1) - rSges(2,2) * t11) + g(2) * (rSges(2,1) * t11 - t10 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t11 + t8 + (-rSges(3,1) - pkin(1)) * t10) + g(2) * (rSges(3,1) * t11 + t10 * rSges(3,3) + t16)) - m(4) * (g(1) * (t8 + (-pkin(1) - pkin(2)) * t10 - t13) + g(2) * (pkin(2) * t11 - t12 + t16)) (-m(3) - m(4)) * (g(1) * t10 - g(2) * t11) -m(4) * (g(1) * t13 + g(2) * t12)];
taug  = t3(:);
