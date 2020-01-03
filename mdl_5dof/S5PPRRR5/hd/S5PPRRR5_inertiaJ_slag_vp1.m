% Calculate joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:40
% DurationCPUTime: 0.41s
% Computational Cost: add. (975->75), mult. (1070->119), div. (0->0), fcn. (1188->8), ass. (0->45)
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t66 = qJ(3) + qJ(4);
t62 = sin(t66);
t63 = cos(t66);
t30 = -t46 * t63 + t47 * t62;
t74 = t30 ^ 2;
t29 = -t46 * t62 - t47 * t63;
t75 = t29 ^ 2;
t81 = t74 + t75;
t51 = cos(qJ(3));
t77 = t51 * pkin(3);
t76 = t30 * t29;
t48 = sin(qJ(5));
t50 = cos(qJ(5));
t38 = -t48 * rSges(6,1) - t50 * rSges(6,2);
t73 = m(6) * t38;
t72 = rSges(6,1) * t50;
t71 = rSges(6,2) * t48;
t49 = sin(qJ(3));
t70 = t46 * t49;
t69 = t47 * t49;
t18 = -t30 * rSges(5,1) + t29 * rSges(5,2);
t65 = t81 * (Icges(6,5) * t48 + Icges(6,6) * t50);
t64 = -t30 * rSges(6,3) - t29 * t71;
t19 = t29 * rSges(5,1) + t30 * rSges(5,2);
t56 = -Icges(6,5) * t50 + Icges(6,6) * t48;
t55 = -t29 * rSges(6,3) + (t71 - t72) * t30;
t54 = t50 ^ 2 * Icges(6,2) + Icges(5,3) + (Icges(6,1) * t48 + 0.2e1 * Icges(6,4) * t50) * t48;
t53 = -pkin(3) * t69 + t77 * t46;
t52 = -pkin(3) * t70 - t77 * t47;
t8 = -t30 * pkin(4) - t29 * pkin(7) + t55;
t9 = -t30 * pkin(7) + (pkin(4) + t72) * t29 + t64;
t33 = t46 * t51 - t69;
t32 = -t47 * t51 - t70;
t21 = t32 * rSges(4,1) - t33 * rSges(4,2);
t20 = t33 * rSges(4,1) + t32 * rSges(4,2);
t17 = t19 + t52;
t16 = t53 + t18;
t11 = Icges(6,3) * t30 + t56 * t29;
t10 = -Icges(6,3) * t29 + t56 * t30;
t5 = t52 + t9;
t4 = t53 + t8;
t1 = t30 * t55 + t29 * (-t29 * t72 - t64);
t2 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t46 ^ 2 + t47 ^ 2); 0; m(4) * (t20 * t46 - t21 * t47) + m(5) * (t16 * t46 - t17 * t47) + m(6) * (t4 * t46 - t5 * t47); Icges(4,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + t54; 0; m(5) * (t18 * t46 - t19 * t47) + m(6) * (t8 * t46 - t9 * t47); m(6) * (t8 * t4 + t9 * t5) + m(5) * (t18 * t16 + t19 * t17) + t54; m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t54; m(6) * t1; (-t29 * t46 + t30 * t47) * t73; (-t29 * t4 - t30 * t5) * t73 + t65; (-t29 * t8 - t30 * t9) * t73 + t65; m(6) * (t81 * t38 ^ 2 + t1 ^ 2) + t30 * (-t10 * t76 + t74 * t11) - t29 * (t75 * t10 - t11 * t76);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
