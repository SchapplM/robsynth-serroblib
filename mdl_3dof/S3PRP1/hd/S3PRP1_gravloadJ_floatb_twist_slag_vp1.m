% Calculate Gravitation load on the joints for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S3PRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:14
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (15->13), mult. (24->18), div. (0->0), fcn. (14->2), ass. (0->5)
t4 = rSges(4,1) + pkin(2);
t3 = rSges(4,3) + qJ(3);
t2 = cos(qJ(2));
t1 = sin(qJ(2));
t5 = [(-m(2) - m(3) - m(4)) * g(2), -m(3) * (g(1) * (-t1 * rSges(3,1) - t2 * rSges(3,2)) + g(2) * (t2 * rSges(3,1) - t1 * rSges(3,2))) - m(4) * ((g(1) * t3 + g(2) * t4) * t2 + (-g(1) * t4 + g(2) * t3) * t1) -m(4) * (g(1) * t1 - g(2) * t2)];
taug  = t5(:);
