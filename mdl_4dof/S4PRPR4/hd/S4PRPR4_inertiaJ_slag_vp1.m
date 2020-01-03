% Calculate joint inertia matrix for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:56
% DurationCPUTime: 0.25s
% Computational Cost: add. (306->48), mult. (324->74), div. (0->0), fcn. (260->4), ass. (0->27)
t30 = pkin(6) + qJ(2);
t28 = sin(t30);
t26 = t28 ^ 2;
t29 = cos(t30);
t27 = t29 ^ 2;
t43 = t26 + t27;
t31 = sin(qJ(4));
t32 = cos(qJ(4));
t49 = rSges(5,1) * t31 + rSges(5,2) * t32;
t48 = t28 * t29;
t47 = t49 * t29;
t44 = t29 * pkin(2) + t28 * qJ(3);
t40 = t29 * rSges(5,3) + t49 * t28;
t34 = Icges(5,5) * t31 + Icges(5,6) * t32;
t23 = t29 * qJ(3);
t2 = t23 + t47 + (-rSges(5,3) - pkin(2) - pkin(5)) * t28;
t3 = t29 * pkin(5) + t40 + t44;
t33 = m(5) * (t28 * t2 - t29 * t3);
t19 = t32 * rSges(5,1) - t31 * rSges(5,2);
t13 = t29 * rSges(3,1) - t28 * rSges(3,2);
t12 = -t28 * rSges(3,1) - t29 * rSges(3,2);
t7 = Icges(5,3) * t28 - t34 * t29;
t6 = Icges(5,3) * t29 + t34 * t28;
t5 = -t29 * rSges(4,2) + t28 * rSges(4,3) + t44;
t4 = t29 * rSges(4,3) + t23 + (rSges(4,2) - pkin(2)) * t28;
t1 = -t28 * t40 + (t28 * rSges(5,3) - t47) * t29;
t8 = [m(2) + m(3) + m(4) + m(5); 0; Icges(5,1) * t32 ^ 2 + Icges(4,1) + Icges(3,3) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Icges(5,4) * t32 + Icges(5,2) * t31) * t31; 0; m(4) * (t28 * t4 - t29 * t5) + t33; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t43; m(5) * t1; t19 * t33 + t43 * (Icges(5,5) * t32 - Icges(5,6) * t31); m(5) * t43 * t19; m(5) * (t43 * t19 ^ 2 + t1 ^ 2) + t29 * (t27 * t6 + t7 * t48) + t28 * (t26 * t7 + t6 * t48);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
