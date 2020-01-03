% Calculate joint inertia matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:35
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.31s
% Computational Cost: add. (306->66), mult. (436->98), div. (0->0), fcn. (346->6), ass. (0->36)
t38 = sin(qJ(1));
t33 = t38 ^ 2;
t39 = cos(qJ(1));
t34 = t39 ^ 2;
t25 = t33 + t34;
t32 = pkin(6) + qJ(4);
t26 = sin(t32);
t27 = cos(t32);
t58 = rSges(5,1) * t26 + rSges(5,2) * t27;
t57 = t38 * t39;
t19 = m(5) * t25;
t35 = sin(pkin(6));
t56 = pkin(3) * t35;
t53 = m(4) * t25 + t19;
t52 = t39 * pkin(1) + t38 * qJ(2);
t49 = rSges(4,3) + qJ(3);
t48 = t39 * rSges(5,3) + t58 * t38;
t36 = cos(pkin(6));
t45 = rSges(4,1) * t35 + rSges(4,2) * t36;
t41 = Icges(5,5) * t26 + Icges(5,6) * t27;
t29 = t39 * qJ(2);
t37 = -pkin(5) - qJ(3);
t2 = t29 + (t58 + t56) * t39 + (-rSges(5,3) - pkin(1) + t37) * t38;
t3 = -t39 * t37 + t38 * t56 + t48 + t52;
t40 = m(5) * (t38 * t2 - t39 * t3);
t22 = t39 * rSges(2,1) - t38 * rSges(2,2);
t21 = -t38 * rSges(2,1) - t39 * rSges(2,2);
t18 = t27 * rSges(5,1) - t26 * rSges(5,2);
t13 = -t39 * rSges(3,2) + t38 * rSges(3,3) + t52;
t12 = t39 * rSges(3,3) + t29 + (rSges(3,2) - pkin(1)) * t38;
t7 = Icges(5,3) * t38 - t39 * t41;
t6 = Icges(5,3) * t39 + t38 * t41;
t5 = t45 * t38 + t49 * t39 + t52;
t4 = t29 + t45 * t39 + (-pkin(1) - t49) * t38;
t1 = -t38 * t48 + (t38 * rSges(5,3) - t39 * t58) * t39;
t8 = [Icges(4,1) * t36 ^ 2 + Icges(5,1) * t27 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t36 + Icges(4,2) * t35) * t35 + m(2) * (t21 ^ 2 + t22 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t2 ^ 2 + t3 ^ 2) + (-0.2e1 * Icges(5,4) * t27 + Icges(5,2) * t26) * t26; m(3) * (t38 * t12 - t39 * t13) + m(4) * (t38 * t4 - t39 * t5) + t40; m(3) * t25 + t53; m(4) * (t38 * t5 + t39 * t4) + m(5) * (t39 * t2 + t38 * t3); 0; t53; t18 * t40 + t25 * (Icges(5,5) * t27 - Icges(5,6) * t26); t18 * t19; 0; m(5) * (t25 * t18 ^ 2 + t1 ^ 2) + t39 * (t34 * t6 + t7 * t57) + t38 * (t33 * t7 + t6 * t57);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
