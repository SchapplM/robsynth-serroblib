% Calculate joint inertia matrix for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:14
% EndTime: 2019-03-08 18:25:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (255->40), mult. (138->46), div. (0->0), fcn. (78->6), ass. (0->25)
t25 = pkin(7) + qJ(2);
t22 = sin(t25);
t27 = pkin(2) * t22;
t26 = Icges(4,3) + Icges(5,3);
t24 = qJ(3) + t25;
t19 = sin(t24);
t20 = cos(t24);
t10 = t20 * rSges(4,1) - t19 * rSges(4,2);
t21 = qJ(4) + t24;
t16 = sin(t21);
t17 = cos(t21);
t8 = t17 * rSges(5,1) - t16 * rSges(5,2);
t4 = pkin(3) * t20 + t8;
t9 = -t19 * rSges(4,1) - t20 * rSges(4,2);
t7 = -t16 * rSges(5,1) - t17 * rSges(5,2);
t3 = -pkin(3) * t19 + t7;
t23 = cos(t25);
t18 = pkin(2) * t23;
t12 = t23 * rSges(3,1) - t22 * rSges(3,2);
t11 = -t22 * rSges(3,1) - t23 * rSges(3,2);
t6 = t10 + t18;
t5 = t9 - t27;
t2 = t18 + t4;
t1 = t3 - t27;
t13 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t26; 0; m(4) * (t10 * t6 + t9 * t5) + m(5) * (t3 * t1 + t4 * t2) + t26; m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t26; 0; m(5) * (t7 * t1 + t8 * t2) + Icges(5,3); m(5) * (t7 * t3 + t8 * t4) + Icges(5,3); m(5) * (t7 ^ 2 + t8 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t13(1) t13(2) t13(4) t13(7); t13(2) t13(3) t13(5) t13(8); t13(4) t13(5) t13(6) t13(9); t13(7) t13(8) t13(9) t13(10);];
Mq  = res;
