% Calculate joint inertia matrix for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [2x2]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_inertiaJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_inertiaJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_inertiaJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (94->27), mult. (226->52), div. (0->0), fcn. (194->4), ass. (0->21)
t21 = cos(qJ(1));
t34 = t21 ^ 2;
t19 = sin(qJ(1));
t35 = t19 ^ 2;
t40 = t34 + t35;
t18 = sin(qJ(2));
t20 = cos(qJ(2));
t25 = -rSges(3,1) * t20 + rSges(3,2) * t18;
t36 = t21 * t19;
t33 = -rSges(3,3) - pkin(1);
t28 = t25 * t21;
t22 = Icges(3,5) * t20 - Icges(3,6) * t18;
t16 = t21 * rSges(2,1) - t19 * rSges(2,2);
t15 = -t19 * rSges(2,1) - t21 * rSges(2,2);
t14 = -t18 * rSges(3,1) - t20 * rSges(3,2);
t5 = -Icges(3,3) * t19 + t22 * t21;
t4 = Icges(3,3) * t21 + t22 * t19;
t3 = t25 * t19 + t33 * t21;
t2 = t33 * t19 - t28;
t1 = t21 * t28 + t25 * t35;
t6 = [t20 ^ 2 * Icges(3,2) + Icges(2,3) + m(2) * (t15 ^ 2 + t16 ^ 2) + m(3) * (t2 ^ 2 + t3 ^ 2) + (Icges(3,1) * t18 + 0.2e1 * Icges(3,4) * t20) * t18; m(3) * (t19 * t2 + t21 * t3) * t14 - t40 * (Icges(3,5) * t18 + Icges(3,6) * t20); m(3) * (t40 * t14 ^ 2 + t1 ^ 2) + t21 * (t34 * t4 - t5 * t36) - t19 * (t35 * t5 - t4 * t36);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t6(1), t6(2); t6(2), t6(3);];
Mq = res;
