% Calculate joint inertia matrix for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:35
% DurationCPUTime: 0.29s
% Computational Cost: add. (608->59), mult. (436->86), div. (0->0), fcn. (360->6), ass. (0->39)
t43 = pkin(7) + qJ(2);
t42 = qJ(3) + t43;
t38 = sin(t42);
t65 = t38 ^ 2;
t39 = cos(t42);
t64 = t39 ^ 2;
t44 = sin(qJ(4));
t68 = Icges(5,5) * t44;
t45 = cos(qJ(4));
t67 = Icges(5,6) * t45;
t27 = t67 + t68;
t66 = t38 * t39;
t30 = t44 * rSges(5,1) + t45 * rSges(5,2);
t61 = m(5) * t30;
t40 = sin(t43);
t60 = pkin(2) * t40;
t59 = rSges(5,1) * t45;
t58 = rSges(5,2) * t44;
t57 = t39 * rSges(5,3) + t38 * t58;
t54 = Icges(5,2) * t45 ^ 2 + Icges(4,3) + (Icges(5,1) * t44 + 0.2e1 * Icges(5,4) * t45) * t44;
t53 = t27 * t64 + (t68 / 0.2e1 + t67 / 0.2e1 + t27 / 0.2e1) * t65;
t19 = t39 * rSges(4,1) - t38 * rSges(4,2);
t18 = -t38 * rSges(4,1) - t39 * rSges(4,2);
t47 = Icges(5,5) * t45 - Icges(5,6) * t44;
t46 = t38 * rSges(5,3) + (-t58 + t59) * t39;
t9 = t39 * pkin(3) + t38 * pkin(6) + t46;
t8 = t39 * pkin(6) + (-pkin(3) - t59) * t38 + t57;
t41 = cos(t43);
t37 = pkin(2) * t41;
t21 = t41 * rSges(3,1) - t40 * rSges(3,2);
t20 = -t40 * rSges(3,1) - t41 * rSges(3,2);
t17 = t19 + t37;
t16 = t18 - t60;
t11 = Icges(5,3) * t38 + t47 * t39;
t10 = -Icges(5,3) * t39 + t47 * t38;
t7 = t37 + t9;
t6 = t8 - t60;
t3 = t38 * (t38 * t59 - t57) + t39 * t46;
t1 = [m(2) + m(3) + m(4) + m(5); 0; Icges(3,3) + m(3) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + t54; 0; m(4) * (t18 * t16 + t19 * t17) + m(5) * (t8 * t6 + t9 * t7) + t54; m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t8 ^ 2 + t9 ^ 2) + t54; m(5) * t3; (-t38 * t7 - t39 * t6) * t61 + t53; (-t38 * t9 - t39 * t8) * t61 + t53; m(5) * (t3 ^ 2 + (t64 + t65) * t30 ^ 2) + t38 * (-t10 * t66 + t65 * t11) - t39 * (t64 * t10 - t11 * t66);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
