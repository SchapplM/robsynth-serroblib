% Calculate joint inertia matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:54
% EndTime: 2019-03-08 18:34:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (316->52), mult. (196->61), div. (0->0), fcn. (118->8), ass. (0->33)
t31 = qJ(1) + qJ(2);
t28 = sin(t31);
t36 = pkin(2) * t28;
t32 = sin(qJ(1));
t35 = t32 * pkin(1);
t34 = Icges(3,3) + Icges(4,3) + Icges(5,3);
t29 = cos(t31);
t14 = t29 * rSges(3,1) - t28 * rSges(3,2);
t27 = pkin(7) + t31;
t26 = qJ(4) + t27;
t20 = sin(t26);
t21 = cos(t26);
t10 = t21 * rSges(5,1) - t20 * rSges(5,2);
t23 = sin(t27);
t24 = cos(t27);
t25 = pkin(2) * t29;
t8 = t24 * rSges(4,1) - t23 * rSges(4,2) + t25;
t13 = -t28 * rSges(3,1) - t29 * rSges(3,2);
t9 = -t20 * rSges(5,1) - t21 * rSges(5,2);
t4 = pkin(3) * t24 + t10 + t25;
t7 = -t23 * rSges(4,1) - t24 * rSges(4,2) - t36;
t3 = -pkin(3) * t23 - t36 + t9;
t33 = cos(qJ(1));
t30 = t33 * pkin(1);
t16 = t33 * rSges(2,1) - t32 * rSges(2,2);
t15 = -t32 * rSges(2,1) - t33 * rSges(2,2);
t12 = t14 + t30;
t11 = t13 - t35;
t6 = t30 + t8;
t5 = t7 - t35;
t2 = t30 + t4;
t1 = t3 - t35;
t17 = [Icges(2,3) + m(2) * (t15 ^ 2 + t16 ^ 2) + m(3) * (t11 ^ 2 + t12 ^ 2) + m(4) * (t5 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t34; m(3) * (t13 * t11 + t14 * t12) + m(4) * (t7 * t5 + t8 * t6) + m(5) * (t3 * t1 + t4 * t2) + t34; m(3) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t34; 0; 0; m(4) + m(5); Icges(5,3) + m(5) * (t9 * t1 + t10 * t2); m(5) * (t10 * t4 + t9 * t3) + Icges(5,3); 0; m(5) * (t10 ^ 2 + t9 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1) t17(2) t17(4) t17(7); t17(2) t17(3) t17(5) t17(8); t17(4) t17(5) t17(6) t17(9); t17(7) t17(8) t17(9) t17(10);];
Mq  = res;
