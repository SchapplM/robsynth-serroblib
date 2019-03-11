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
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
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
% StartTime: 2019-03-08 18:00:00
% EndTime: 2019-03-08 18:00:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (94->29), mult. (226->54), div. (0->0), fcn. (194->4), ass. (0->21)
t22 = cos(qJ(1));
t34 = t22 ^ 2;
t20 = sin(qJ(1));
t35 = t20 ^ 2;
t42 = t34 + t35;
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t38 = -rSges(3,1) * t21 + rSges(3,2) * t19;
t37 = t20 * t22;
t36 = t38 * t20;
t24 = -Icges(3,5) * t21 + Icges(3,6) * t19;
t23 = t20 * rSges(3,3) + t38 * t22;
t16 = -rSges(2,1) * t22 + t20 * rSges(2,2);
t15 = t20 * rSges(2,1) + rSges(2,2) * t22;
t14 = -rSges(3,1) * t19 - rSges(3,2) * t21;
t5 = Icges(3,3) * t20 + t24 * t22;
t4 = -Icges(3,3) * t22 + t24 * t20;
t3 = (rSges(3,3) + pkin(1)) * t22 - t36;
t2 = t20 * pkin(1) + t23;
t1 = t22 * t23 + (-t22 * rSges(3,3) + t36) * t20;
t6 = [t21 ^ 2 * Icges(3,2) + Icges(2,3) + m(2) * (t15 ^ 2 + t16 ^ 2) + m(3) * (t2 ^ 2 + t3 ^ 2) + (Icges(3,1) * t19 + 0.2e1 * Icges(3,4) * t21) * t19; m(3) * (-t2 * t20 - t22 * t3) * t14 - t42 * (Icges(3,5) * t19 + Icges(3,6) * t21); m(3) * (t42 * t14 ^ 2 + t1 ^ 2) - t22 * (t4 * t34 - t5 * t37) + t20 * (t35 * t5 - t4 * t37);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t6(1) t6(2); t6(2) t6(3);];
Mq  = res;
