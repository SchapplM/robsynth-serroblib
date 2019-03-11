% Calculate joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:03
% EndTime: 2019-03-08 18:36:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (366->56), mult. (246->68), div. (0->0), fcn. (150->6), ass. (0->34)
t37 = rSges(5,1) + pkin(3);
t30 = qJ(1) + qJ(2);
t26 = sin(t30);
t36 = pkin(2) * t26;
t31 = sin(qJ(1));
t35 = t31 * pkin(1);
t34 = Icges(4,3) + Icges(5,3);
t33 = Icges(3,3) + t34;
t27 = cos(t30);
t16 = t27 * rSges(3,1) - t26 * rSges(3,2);
t28 = qJ(3) + t30;
t24 = sin(t28);
t25 = cos(t28);
t14 = t25 * rSges(4,1) - t24 * rSges(4,2);
t23 = pkin(2) * t27;
t10 = t14 + t23;
t8 = -t24 * rSges(5,2) + t37 * t25;
t15 = -t26 * rSges(3,1) - t27 * rSges(3,2);
t13 = -t24 * rSges(4,1) - t25 * rSges(4,2);
t4 = t23 + t8;
t7 = -t25 * rSges(5,2) - t37 * t24;
t9 = t13 - t36;
t3 = t7 - t36;
t32 = cos(qJ(1));
t29 = t32 * pkin(1);
t18 = t32 * rSges(2,1) - t31 * rSges(2,2);
t17 = -t31 * rSges(2,1) - t32 * rSges(2,2);
t12 = t16 + t29;
t11 = t15 - t35;
t6 = t10 + t29;
t5 = t9 - t35;
t2 = t29 + t4;
t1 = t3 - t35;
t19 = [Icges(2,3) + m(2) * (t17 ^ 2 + t18 ^ 2) + m(3) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t33; m(3) * (t15 * t11 + t16 * t12) + m(4) * (t10 * t6 + t9 * t5) + m(5) * (t3 * t1 + t4 * t2) + t33; m(3) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t33; m(4) * (t13 * t5 + t14 * t6) + m(5) * (t7 * t1 + t8 * t2) + t34; m(4) * (t14 * t10 + t13 * t9) + m(5) * (t7 * t3 + t8 * t4) + t34; m(4) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t7 ^ 2 + t8 ^ 2) + t34; 0; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1) t19(2) t19(4) t19(7); t19(2) t19(3) t19(5) t19(8); t19(4) t19(5) t19(6) t19(9); t19(7) t19(8) t19(9) t19(10);];
Mq  = res;
