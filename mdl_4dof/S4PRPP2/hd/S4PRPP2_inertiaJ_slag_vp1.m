% Calculate joint inertia matrix for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:49
% EndTime: 2019-03-08 18:18:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->29), mult. (87->33), div. (0->0), fcn. (45->4), ass. (0->17)
t16 = m(4) + m(5);
t11 = sin(qJ(2));
t15 = pkin(2) * t11;
t14 = rSges(5,1) + pkin(3);
t13 = rSges(5,3) + qJ(4);
t12 = cos(qJ(2));
t10 = qJ(2) + pkin(5);
t9 = t12 * pkin(2);
t8 = cos(t10);
t7 = sin(t10);
t6 = rSges(3,1) * t12 - t11 * rSges(3,2);
t5 = -t11 * rSges(3,1) - rSges(3,2) * t12;
t4 = rSges(4,1) * t8 - rSges(4,2) * t7 + t9;
t3 = -rSges(4,1) * t7 - rSges(4,2) * t8 - t15;
t2 = t13 * t7 + t14 * t8 + t9;
t1 = t13 * t8 - t14 * t7 - t15;
t17 = [m(2) + m(3) + t16; m(3) * t6 + m(4) * t4 + m(5) * t2; Icges(3,3) + Icges(4,3) + Icges(5,2) + m(3) * (t5 ^ 2 + t6 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2); 0; 0; t16; -m(5) * t8; m(5) * (t1 * t7 - t2 * t8); 0; m(5) * (t7 ^ 2 + t8 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1) t17(2) t17(4) t17(7); t17(2) t17(3) t17(5) t17(8); t17(4) t17(5) t17(6) t17(9); t17(7) t17(8) t17(9) t17(10);];
Mq  = res;
