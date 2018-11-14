% Calculate joint inertia matrix for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3PRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_inertiaJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR2_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR2_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR2_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:45
% EndTime: 2018-11-14 10:12:45
% DurationCPUTime: 0.07s
% Computational Cost: add. (50->19), mult. (68->25), div. (0->0), fcn. (31->4), ass. (0->12)
t10 = qJ(2) + qJ(3);
t8 = sin(t10);
t9 = cos(t10);
t4 = t9 * rSges(4,1) - t8 * rSges(4,2);
t3 = -t8 * rSges(4,1) - t9 * rSges(4,2);
t12 = cos(qJ(2));
t11 = sin(qJ(2));
t6 = t12 * rSges(3,1) - t11 * rSges(3,2);
t5 = -t11 * rSges(3,1) - t12 * rSges(3,2);
t2 = t12 * pkin(2) + t4;
t1 = -t11 * pkin(2) + t3;
t7 = [m(2) + m(3) + m(4); m(3) * t6 + m(4) * t2; Icges(3,3) + Icges(4,3) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2); m(4) * t4; m(4) * (t3 * t1 + t4 * t2) + Icges(4,3); m(4) * (t3 ^ 2 + t4 ^ 2) + Icges(4,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t7(1) t7(2) t7(4); t7(2) t7(3) t7(5); t7(4) t7(5) t7(6);];
Mq  = res;
