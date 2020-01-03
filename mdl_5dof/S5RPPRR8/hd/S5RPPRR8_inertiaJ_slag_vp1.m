% Calculate joint inertia matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:05
% DurationCPUTime: 0.45s
% Computational Cost: add. (1021->111), mult. (1134->163), div. (0->0), fcn. (1236->8), ass. (0->58)
t53 = sin(qJ(5));
t90 = t53 / 0.2e1;
t55 = cos(qJ(5));
t89 = t55 / 0.2e1;
t88 = Icges(6,5) * t90 + Icges(6,6) * t89;
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t73 = pkin(8) + qJ(4);
t67 = sin(t73);
t68 = cos(t73);
t28 = -t54 * t67 - t56 * t68;
t29 = -t54 * t68 + t56 * t67;
t87 = t29 * t28;
t81 = rSges(6,2) * t53;
t70 = -pkin(4) + t81;
t82 = rSges(6,1) * t55;
t77 = t29 * rSges(6,3) - t28 * t82;
t7 = -t29 * pkin(7) - t70 * t28 - t77;
t78 = t28 * rSges(6,3) + t29 * t82;
t6 = -t28 * pkin(7) + t70 * t29 - t78;
t74 = Icges(6,4) * t55;
t60 = Icges(6,2) * t53 - t74;
t75 = Icges(6,4) * t53;
t61 = -Icges(6,1) * t55 + t75;
t86 = -(t28 * t88 - t55 * (-Icges(6,6) * t28 + t60 * t29) / 0.2e1 - t53 * (-Icges(6,5) * t28 + t61 * t29) / 0.2e1) * t28 - (t29 * t88 + (Icges(6,6) * t29 + t60 * t28) * t89 + (Icges(6,5) * t29 + t61 * t28) * t90) * t29;
t85 = t28 ^ 2;
t84 = t29 ^ 2;
t37 = -t53 * rSges(6,1) - t55 * rSges(6,2);
t83 = m(6) * t37;
t51 = sin(pkin(8));
t80 = t54 * t51;
t79 = t56 * t51;
t76 = t56 * pkin(1) + t54 * qJ(2);
t52 = cos(pkin(8));
t47 = t52 * pkin(3) + pkin(2);
t69 = pkin(3) * t80 + t56 * t47 + t76;
t20 = -t29 * rSges(5,1) + t28 * rSges(5,2);
t21 = t28 * rSges(5,1) + t29 * rSges(5,2);
t59 = -Icges(6,5) * t55 + Icges(6,6) * t53;
t49 = t56 * qJ(2);
t58 = pkin(3) * t79 + t49 + (-pkin(1) - t47) * t54;
t57 = t55 * (-Icges(6,2) * t55 - t75) + t53 * (-Icges(6,1) * t53 - t74) - Icges(5,3);
t39 = t56 * rSges(2,1) - t54 * rSges(2,2);
t38 = -t54 * rSges(2,1) - t56 * rSges(2,2);
t32 = t54 * t52 - t79;
t31 = -t56 * t52 - t80;
t27 = t56 * rSges(3,1) + t54 * rSges(3,3) + t76;
t26 = t56 * rSges(3,3) + t49 + (-rSges(3,1) - pkin(1)) * t54;
t19 = -t31 * rSges(4,1) + t32 * rSges(4,2) + t56 * pkin(2) + t76;
t18 = -t32 * rSges(4,1) - t31 * rSges(4,2) + t49 + (-pkin(1) - pkin(2)) * t54;
t13 = Icges(6,3) * t29 + t59 * t28;
t12 = -Icges(6,3) * t28 + t59 * t29;
t11 = -t21 + t69;
t10 = -t20 + t58;
t5 = t69 - t7;
t4 = t58 - t6;
t1 = t29 * (t29 * t81 - t78) + t28 * (t28 * t81 + t77);
t2 = [Icges(3,2) + Icges(2,3) + Icges(4,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2) + m(2) * (t38 ^ 2 + t39 ^ 2) - t57; m(6) * (t54 * t4 - t56 * t5) + m(5) * (t54 * t10 - t56 * t11) + m(3) * (t54 * t26 - t56 * t27) + m(4) * (t54 * t18 - t56 * t19); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t54 ^ 2 + t56 ^ 2); 0; 0; m(4) + m(5) + m(6); m(6) * (t6 * t4 + t7 * t5) + m(5) * (t20 * t10 + t21 * t11) + t57; m(5) * (t20 * t54 - t21 * t56) + m(6) * (t6 * t54 - t7 * t56); 0; m(5) * (t20 ^ 2 + t21 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) - t57; (-t28 * t4 - t29 * t5) * t83 + t86; (-t28 * t54 + t29 * t56) * t83; -m(6) * t1; (-t28 * t6 - t29 * t7) * t83 - t86; m(6) * (t1 ^ 2 + (t84 + t85) * t37 ^ 2) + t29 * (-t12 * t87 + t84 * t13) - t28 * (t85 * t12 - t13 * t87);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
