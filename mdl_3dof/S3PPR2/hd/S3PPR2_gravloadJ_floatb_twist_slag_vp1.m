% Calculate Gravitation load on the joints for
% S3PPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S3PPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:00
% EndTime: 2018-11-14 10:11:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (7->6), mult. (12->9), div. (0->0), fcn. (4->2), ass. (0->4)
t3 = m(3) + m(4);
t2 = cos(qJ(3));
t1 = sin(qJ(3));
t4 = [(-m(2) - t3) * g(1), t3 * g(2), -m(4) * (g(1) * (-rSges(4,1) * t1 - rSges(4,2) * t2) + g(2) * (-rSges(4,1) * t2 + rSges(4,2) * t1))];
taug  = t4(:);
