% Calculate joint inertia matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:34:03
% EndTime: 2019-05-28 15:34:03
% DurationCPUTime: 0.17s
% Computational Cost: add. (449->63), mult. (398->87), div. (0->0), fcn. (340->6), ass. (0->31)
t36 = sin(qJ(1));
t42 = t36 * pkin(1);
t41 = cos(qJ(4));
t40 = sin(qJ(4));
t35 = qJ(1) + qJ(2);
t32 = sin(t35);
t33 = cos(t35);
t15 = -t32 * t40 - t33 * t41;
t16 = t32 * t41 - t33 * t40;
t6 = t15 * rSges(5,1) - t16 * rSges(5,2);
t39 = t33 * pkin(2) + t32 * qJ(3);
t38 = Icges(3,3) + Icges(4,2) + Icges(5,3);
t18 = t33 * rSges(3,1) - t32 * rSges(3,2);
t10 = t33 * rSges(4,1) + t32 * rSges(4,3) + t39;
t4 = t33 * pkin(3) + t39 - t6;
t17 = -t32 * rSges(3,1) - t33 * rSges(3,2);
t5 = t16 * rSges(5,1) + t15 * rSges(5,2);
t25 = t33 * qJ(3);
t9 = t25 + t33 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t32;
t3 = t25 + (-pkin(2) - pkin(3)) * t32 - t5;
t37 = cos(qJ(1));
t34 = t37 * pkin(1);
t21 = t37 * rSges(2,1) - t36 * rSges(2,2);
t20 = -t36 * rSges(2,1) - t37 * rSges(2,2);
t14 = t18 + t34;
t13 = t17 - t42;
t8 = t34 + t10;
t7 = t9 - t42;
t2 = t34 + t4;
t1 = t3 - t42;
t11 = [Icges(2,3) + m(2) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + t38; m(3) * (t17 * t13 + t18 * t14) + m(4) * (t10 * t8 + t9 * t7) + m(5) * (t3 * t1 + t4 * t2) + t38; m(3) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + t38; m(4) * (t32 * t7 - t33 * t8) + m(5) * (t32 * t1 - t33 * t2); m(4) * (-t33 * t10 + t32 * t9) + m(5) * (t32 * t3 - t33 * t4); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * (t32 ^ 2 + t33 ^ 2); -Icges(5,3) + m(5) * (t5 * t1 + t6 * t2); m(5) * (t5 * t3 + t6 * t4) - Icges(5,3); m(5) * (t5 * t32 - t6 * t33); m(5) * (t5 ^ 2 + t6 ^ 2) + Icges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1), t11(2), t11(4), t11(7); t11(2), t11(3), t11(5), t11(8); t11(4), t11(5), t11(6), t11(9); t11(7), t11(8), t11(9), t11(10);];
Mq  = res;
