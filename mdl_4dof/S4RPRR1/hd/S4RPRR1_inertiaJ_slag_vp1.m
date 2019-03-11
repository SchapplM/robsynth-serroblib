% Calculate joint inertia matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:39
% EndTime: 2019-03-08 18:31:39
% DurationCPUTime: 0.10s
% Computational Cost: add. (271->48), mult. (160->55), div. (0->0), fcn. (94->8), ass. (0->31)
t29 = sin(qJ(1));
t34 = t29 * pkin(1);
t28 = qJ(1) + pkin(7);
t25 = cos(t28);
t30 = cos(qJ(1));
t27 = t30 * pkin(1);
t33 = pkin(2) * t25 + t27;
t32 = Icges(4,3) + Icges(5,3);
t26 = qJ(3) + t28;
t21 = sin(t26);
t22 = cos(t26);
t12 = t22 * rSges(4,1) - t21 * rSges(4,2);
t23 = qJ(4) + t26;
t18 = sin(t23);
t19 = cos(t23);
t8 = t19 * rSges(5,1) - t18 * rSges(5,2);
t6 = pkin(3) * t22 + t8;
t24 = sin(t28);
t31 = -pkin(2) * t24 - t34;
t11 = -t21 * rSges(4,1) - t22 * rSges(4,2);
t7 = -t18 * rSges(5,1) - t19 * rSges(5,2);
t5 = -pkin(3) * t21 + t7;
t14 = t30 * rSges(2,1) - t29 * rSges(2,2);
t13 = -t29 * rSges(2,1) - t30 * rSges(2,2);
t10 = t25 * rSges(3,1) - t24 * rSges(3,2) + t27;
t9 = -t24 * rSges(3,1) - t25 * rSges(3,2) - t34;
t4 = t12 + t33;
t3 = t11 + t31;
t2 = t6 + t33;
t1 = t31 + t5;
t15 = [Icges(2,3) + Icges(3,3) + m(2) * (t13 ^ 2 + t14 ^ 2) + m(3) * (t10 ^ 2 + t9 ^ 2) + m(4) * (t3 ^ 2 + t4 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t32; 0; m(3) + m(4) + m(5); m(4) * (t11 * t3 + t12 * t4) + m(5) * (t5 * t1 + t6 * t2) + t32; 0; m(4) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + t32; Icges(5,3) + m(5) * (t7 * t1 + t8 * t2); 0; m(5) * (t7 * t5 + t8 * t6) + Icges(5,3); m(5) * (t7 ^ 2 + t8 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;
