% Calculate joint inertia matrix for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_inertiaJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (177->39), mult. (138->46), div. (0->0), fcn. (78->6), ass. (0->24)
t21 = sin(qJ(2));
t25 = t21 * pkin(1);
t22 = cos(qJ(2));
t24 = t22 * pkin(1);
t23 = Icges(4,3) + Icges(5,3);
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t10 = -t18 * rSges(4,1) + t17 * rSges(4,2);
t19 = qJ(4) + t20;
t15 = sin(t19);
t16 = cos(t19);
t8 = -t16 * rSges(5,1) + t15 * rSges(5,2);
t9 = -t17 * rSges(4,1) - t18 * rSges(4,2);
t7 = -t15 * rSges(5,1) - t16 * rSges(5,2);
t4 = -pkin(2) * t18 + t8;
t3 = -pkin(2) * t17 + t7;
t12 = -t22 * rSges(3,1) + t21 * rSges(3,2);
t11 = -t21 * rSges(3,1) - t22 * rSges(3,2);
t6 = t10 - t24;
t5 = t9 - t25;
t2 = t4 - t24;
t1 = t3 - t25;
t13 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t23; 0; m(4) * (t10 * t6 + t9 * t5) + m(5) * (t3 * t1 + t4 * t2) + t23; m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t23; 0; Icges(5,3) + m(5) * (t7 * t1 + t8 * t2); m(5) * (t7 * t3 + t8 * t4) + Icges(5,3); m(5) * (t7 ^ 2 + t8 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1), t13(2), t13(4), t13(7); t13(2), t13(3), t13(5), t13(8); t13(4), t13(5), t13(6), t13(9); t13(7), t13(8), t13(9), t13(10);];
Mq  = res;
