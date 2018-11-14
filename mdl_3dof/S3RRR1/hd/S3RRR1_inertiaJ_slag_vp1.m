% Calculate joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:53
% EndTime: 2018-11-14 10:15:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (174->36), mult. (138->46), div. (0->0), fcn. (78->6), ass. (0->24)
t23 = sin(qJ(1));
t26 = t23 * pkin(1);
t25 = Icges(3,3) + Icges(4,3);
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t19 = cos(t22);
t10 = t19 * rSges(3,1) - t18 * rSges(3,2);
t20 = qJ(3) + t22;
t16 = sin(t20);
t17 = cos(t20);
t8 = t17 * rSges(4,1) - t16 * rSges(4,2);
t4 = pkin(2) * t19 + t8;
t9 = -t18 * rSges(3,1) - t19 * rSges(3,2);
t7 = -t16 * rSges(4,1) - t17 * rSges(4,2);
t3 = -pkin(2) * t18 + t7;
t24 = cos(qJ(1));
t21 = t24 * pkin(1);
t12 = t24 * rSges(2,1) - t23 * rSges(2,2);
t11 = -t23 * rSges(2,1) - t24 * rSges(2,2);
t6 = t10 + t21;
t5 = t9 - t26;
t2 = t21 + t4;
t1 = t3 - t26;
t13 = [Icges(2,3) + m(2) * (t11 ^ 2 + t12 ^ 2) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2) + t25; m(3) * (t10 * t6 + t5 * t9) + m(4) * (t1 * t3 + t2 * t4) + t25; m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + t25; m(4) * (t1 * t7 + t2 * t8) + Icges(4,3); m(4) * (t3 * t7 + t4 * t8) + Icges(4,3); m(4) * (t7 ^ 2 + t8 ^ 2) + Icges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t13(1) t13(2) t13(4); t13(2) t13(3) t13(5); t13(4) t13(5) t13(6);];
Mq  = res;
