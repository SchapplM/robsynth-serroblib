% Calculate joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:22
% DurationCPUTime: 0.39s
% Computational Cost: add. (559->94), mult. (898->138), div. (0->0), fcn. (976->8), ass. (0->49)
t43 = pkin(8) + qJ(5);
t37 = sin(t43);
t74 = Icges(6,5) * t37;
t38 = cos(t43);
t73 = Icges(6,6) * t38;
t72 = -t74 / 0.2e1 - t73 / 0.2e1;
t58 = m(5) / 0.2e1 + m(6) / 0.2e1;
t71 = 0.2e1 * t58;
t47 = sin(qJ(1));
t48 = cos(qJ(1));
t59 = sin(pkin(7));
t60 = cos(pkin(7));
t27 = -t47 * t59 - t48 * t60;
t28 = -t47 * t60 + t48 * t59;
t70 = t28 * t27;
t26 = t27 ^ 2;
t25 = t28 ^ 2;
t24 = -t37 * rSges(6,1) - t38 * rSges(6,2);
t69 = m(6) * t24;
t68 = rSges(6,1) * t38;
t67 = rSges(6,2) * t37;
t66 = t28 * rSges(6,3) - t27 * t68;
t65 = t25 + t26;
t64 = t48 * pkin(1) + t47 * qJ(2);
t61 = rSges(5,3) + qJ(4);
t57 = t48 * pkin(2) + t64;
t40 = t48 * qJ(2);
t56 = t40 + (-pkin(1) - pkin(2)) * t47;
t53 = t67 - t68;
t50 = -Icges(6,5) * t38 + Icges(6,6) * t37;
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t49 = rSges(5,1) * t45 - rSges(5,2) * t44 + pkin(3);
t46 = -pkin(6) - qJ(4);
t36 = t45 * pkin(4) + pkin(3);
t31 = t48 * rSges(2,1) - t47 * rSges(2,2);
t30 = -t47 * rSges(2,1) - t48 * rSges(2,2);
t18 = t48 * rSges(3,1) + t47 * rSges(3,3) + t64;
t17 = t48 * rSges(3,3) + t40 + (-rSges(3,1) - pkin(1)) * t47;
t13 = -t27 * rSges(4,1) - t28 * rSges(4,2) + t57;
t12 = t28 * rSges(4,1) - t27 * rSges(4,2) + t56;
t7 = Icges(6,3) * t28 + t50 * t27;
t6 = -Icges(6,3) * t27 + t50 * t28;
t5 = -t49 * t27 + t61 * t28 + t57;
t4 = t61 * t27 + t49 * t28 + t56;
t3 = -t28 * t46 + (-t36 + t67) * t27 + t57 + t66;
t2 = (rSges(6,3) - t46) * t27 + (t36 - t53) * t28 + t56;
t1 = t27 * (t27 * t67 + t66) + (-t27 * rSges(6,3) + t53 * t28) * t28;
t8 = [t45 ^ 2 * Icges(5,2) + t38 ^ 2 * Icges(6,2) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(5,1) * t44 + 0.2e1 * Icges(5,4) * t45) * t44 + m(2) * (t30 ^ 2 + t31 ^ 2) + m(3) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t4 ^ 2 + t5 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2) + (Icges(6,1) * t37 + 0.2e1 * Icges(6,4) * t38) * t37; m(3) * (t47 * t17 - t48 * t18) + m(4) * (t47 * t12 - t48 * t13) + m(5) * (t47 * t4 - t48 * t5) + m(6) * (t47 * t2 - t48 * t3); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t58) * (t47 ^ 2 + t48 ^ 2); 0; 0; m(4) + m(5) + m(6); m(5) * (-t27 * t5 + t28 * t4) + m(6) * (t28 * t2 - t27 * t3); (t27 * t48 + t28 * t47) * t71; 0; t65 * t71; (t25 / 0.2e1 + t26 / 0.2e1) * (-t73 - t74) + (t72 * t28 - t3 * t69) * t28 + (-t2 * t69 + t72 * t27) * t27; (-t27 * t47 + t28 * t48) * t69; -m(6) * t1; 0; m(6) * (t65 * t24 ^ 2 + t1 ^ 2) + t28 * (t25 * t7 - t6 * t70) - t27 * (t26 * t6 - t7 * t70);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
