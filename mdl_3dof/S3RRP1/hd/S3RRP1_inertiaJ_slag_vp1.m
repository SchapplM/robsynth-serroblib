% Calculate joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->32), mult. (138->44), div. (0->0), fcn. (84->4), ass. (0->21)
t27 = rSges(4,1) + pkin(2);
t26 = qJ(3) + rSges(4,3);
t22 = sin(qJ(1));
t25 = pkin(1) * t22;
t24 = Icges(3,3) + Icges(4,2);
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t19 = cos(t21);
t8 = t19 * rSges(3,1) - rSges(3,2) * t18;
t4 = t26 * t18 + t27 * t19;
t7 = -rSges(3,1) * t18 - rSges(3,2) * t19;
t3 = -t27 * t18 + t26 * t19;
t23 = cos(qJ(1));
t20 = t23 * pkin(1);
t10 = rSges(2,1) * t23 - t22 * rSges(2,2);
t9 = -t22 * rSges(2,1) - rSges(2,2) * t23;
t6 = t20 + t8;
t5 = t7 - t25;
t2 = t20 + t4;
t1 = t3 - t25;
t11 = [Icges(2,3) + m(2) * (t10 ^ 2 + t9 ^ 2) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2) + t24; m(3) * (t7 * t5 + t8 * t6) + m(4) * (t1 * t3 + t4 * t2) + t24; m(3) * (t7 ^ 2 + t8 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + t24; m(4) * (t1 * t18 - t19 * t2); m(4) * (t18 * t3 - t19 * t4); m(4) * (t18 ^ 2 + t19 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t11(1) t11(2) t11(4); t11(2) t11(3) t11(5); t11(4) t11(5) t11(6);];
Mq  = res;
