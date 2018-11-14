% Calculate joint inertia matrix for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:23
% EndTime: 2018-11-14 10:14:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->35), mult. (170->51), div. (0->0), fcn. (144->4), ass. (0->17)
t21 = cos(qJ(3));
t20 = sin(qJ(3));
t17 = sin(qJ(1));
t18 = cos(qJ(1));
t19 = t18 * pkin(1) + t17 * qJ(2);
t7 = -t17 * t20 - t18 * t21;
t8 = t17 * t21 - t18 * t20;
t3 = t8 * rSges(4,1) + t7 * rSges(4,2);
t4 = t7 * rSges(4,1) - t8 * rSges(4,2);
t15 = t18 * qJ(2);
t10 = t18 * rSges(2,1) - t17 * rSges(2,2);
t9 = -t17 * rSges(2,1) - t18 * rSges(2,2);
t6 = t18 * rSges(3,1) + t17 * rSges(3,3) + t19;
t5 = t18 * rSges(3,3) + t15 + (-rSges(3,1) - pkin(1)) * t17;
t2 = t18 * pkin(2) + t19 - t4;
t1 = t15 + (-pkin(1) - pkin(2)) * t17 - t3;
t11 = [Icges(2,3) + Icges(3,2) + Icges(4,3) + m(2) * (t10 ^ 2 + t9 ^ 2) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2); m(3) * (t17 * t5 - t18 * t6) + m(4) * (t17 * t1 - t18 * t2); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * (t17 ^ 2 + t18 ^ 2); m(4) * (t3 * t1 + t2 * t4) - Icges(4,3); m(4) * (t3 * t17 - t4 * t18); m(4) * (t3 ^ 2 + t4 ^ 2) + Icges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t11(1) t11(2) t11(4); t11(2) t11(3) t11(5); t11(4) t11(5) t11(6);];
Mq  = res;
