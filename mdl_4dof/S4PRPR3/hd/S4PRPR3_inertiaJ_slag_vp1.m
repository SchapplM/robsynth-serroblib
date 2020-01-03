% Calculate joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:47
% EndTime: 2019-12-31 16:20:49
% DurationCPUTime: 0.28s
% Computational Cost: add. (406->54), mult. (330->84), div. (0->0), fcn. (272->6), ass. (0->32)
t30 = pkin(6) + qJ(2);
t26 = sin(t30);
t23 = t26 ^ 2;
t28 = cos(t30);
t24 = t28 ^ 2;
t44 = t23 + t24;
t48 = t26 * t28;
t29 = pkin(7) + qJ(4);
t27 = cos(t29);
t47 = rSges(5,1) * t27;
t25 = sin(t29);
t46 = rSges(5,2) * t25;
t45 = t26 * rSges(5,3) + t28 * t47;
t41 = rSges(4,3) + qJ(3);
t38 = -t46 + t47;
t35 = Icges(5,5) * t27 - Icges(5,6) * t25;
t31 = sin(pkin(7));
t32 = cos(pkin(7));
t34 = rSges(4,1) * t32 - rSges(4,2) * t31 + pkin(2);
t33 = -pkin(5) - qJ(3);
t22 = t32 * pkin(3) + pkin(2);
t18 = t28 * rSges(3,1) - t26 * rSges(3,2);
t17 = -t26 * rSges(3,1) - t28 * rSges(3,2);
t16 = t25 * rSges(5,1) + t27 * rSges(5,2);
t7 = Icges(5,3) * t26 + t28 * t35;
t6 = -Icges(5,3) * t28 + t26 * t35;
t5 = t26 * t41 + t34 * t28;
t4 = -t34 * t26 + t28 * t41;
t3 = -t26 * t33 + (t22 - t46) * t28 + t45;
t2 = (rSges(5,3) - t33) * t28 + (-t22 - t38) * t26;
t1 = t28 * (-t28 * t46 + t45) + (-t28 * rSges(5,3) + t26 * t38) * t26;
t8 = [m(2) + m(3) + m(4) + m(5); 0; Icges(4,2) * t32 ^ 2 + Icges(5,2) * t27 ^ 2 + Icges(3,3) + (Icges(4,1) * t31 + 0.2e1 * Icges(4,4) * t32) * t31 + m(3) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (Icges(5,1) * t25 + 0.2e1 * Icges(5,4) * t27) * t25; 0; m(4) * (t26 * t4 - t28 * t5) + m(5) * (t26 * t2 - t28 * t3); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t44; m(5) * t1; m(5) * (-t2 * t28 - t26 * t3) * t16 + t44 * (Icges(5,5) * t25 + Icges(5,6) * t27); 0; m(5) * (t44 * t16 ^ 2 + t1 ^ 2) + t26 * (t23 * t7 - t6 * t48) - t28 * (t24 * t6 - t7 * t48);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
