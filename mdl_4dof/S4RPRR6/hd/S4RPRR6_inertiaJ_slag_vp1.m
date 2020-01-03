% Calculate joint inertia matrix for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:30
% DurationCPUTime: 0.52s
% Computational Cost: add. (1015->125), mult. (958->196), div. (0->0), fcn. (840->8), ass. (0->66)
t65 = sin(qJ(1));
t66 = cos(qJ(1));
t98 = t65 * t66;
t59 = pkin(7) + qJ(3);
t54 = qJ(4) + t59;
t49 = sin(t54);
t50 = cos(t54);
t80 = rSges(5,1) * t50 - rSges(5,2) * t49;
t60 = t65 ^ 2;
t61 = t66 ^ 2;
t97 = t65 / 0.2e1;
t96 = -t66 / 0.2e1;
t63 = cos(pkin(7));
t51 = t63 * pkin(2) + pkin(1);
t53 = cos(t59);
t95 = rSges(4,1) * t53;
t52 = sin(t59);
t93 = rSges(4,2) * t52;
t64 = -pkin(5) - qJ(2);
t68 = t65 * rSges(5,3) + t80 * t66;
t8 = t65 * (-t66 * rSges(5,3) + t80 * t65) + t66 * t68;
t91 = t65 * rSges(4,3) + t66 * t95;
t90 = t60 + t61;
t89 = Icges(4,4) * t52;
t88 = Icges(4,4) * t53;
t87 = Icges(5,4) * t49;
t86 = Icges(5,4) * t50;
t85 = rSges(3,3) + qJ(2);
t33 = Icges(5,5) * t49 + Icges(5,6) * t50;
t71 = -Icges(5,2) * t49 + t86;
t73 = Icges(5,1) * t50 - t87;
t34 = Icges(5,2) * t50 + t87;
t35 = Icges(5,1) * t49 + t86;
t75 = -t34 * t49 + t35 * t50;
t84 = (t50 * (Icges(5,6) * t65 + t71 * t66) + t49 * (Icges(5,5) * t65 + t73 * t66) + t65 * t33 + t75 * t66) * t97 + (t50 * (-Icges(5,6) * t66 + t71 * t65) + t49 * (-Icges(5,5) * t66 + t73 * t65) - t66 * t33 + t75 * t65) * t96;
t69 = Icges(5,5) * t50 - Icges(5,6) * t49;
t20 = -Icges(5,3) * t66 + t69 * t65;
t21 = Icges(5,3) * t65 + t69 * t66;
t83 = -t66 * (t61 * t20 - t21 * t98) + t65 * (-t20 * t98 + t60 * t21);
t36 = t49 * rSges(5,1) + t50 * rSges(5,2);
t82 = -pkin(3) * t52 - t36;
t81 = -t93 + t95;
t74 = Icges(4,1) * t53 - t89;
t72 = -Icges(4,2) * t52 + t88;
t70 = Icges(4,5) * t53 - Icges(4,6) * t52;
t62 = sin(pkin(7));
t67 = rSges(3,1) * t63 - rSges(3,2) * t62 + pkin(1);
t58 = -pkin(6) + t64;
t46 = t66 * rSges(2,1) - t65 * rSges(2,2);
t45 = -t65 * rSges(2,1) - t66 * rSges(2,2);
t43 = pkin(3) * t53 + t51;
t42 = t52 * rSges(4,1) + t53 * rSges(4,2);
t37 = t66 * t43;
t27 = Icges(4,3) * t65 + t70 * t66;
t26 = -Icges(4,3) * t66 + t70 * t65;
t19 = t85 * t65 + t67 * t66;
t18 = -t67 * t65 + t85 * t66;
t15 = t82 * t66;
t14 = t82 * t65;
t13 = -t65 * t64 + (t51 - t93) * t66 + t91;
t12 = (rSges(4,3) - t64) * t66 + (-t51 - t81) * t65;
t11 = -t65 * t58 + t37 + t68;
t10 = (rSges(5,3) - t58) * t66 + (-t43 - t80) * t65;
t9 = t66 * (-t66 * t93 + t91) + (-t66 * rSges(4,3) + t81 * t65) * t65;
t3 = t66 * (-t66 * t51 + t37) + (t43 - t51) * t60 + t8;
t1 = [Icges(3,2) * t63 ^ 2 + t50 * t34 + t49 * t35 + t53 * (Icges(4,2) * t53 + t89) + t52 * (Icges(4,1) * t52 + t88) + Icges(2,3) + (Icges(3,1) * t62 + 0.2e1 * Icges(3,4) * t63) * t62 + m(5) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(3) * (t18 ^ 2 + t19 ^ 2) + m(2) * (t45 ^ 2 + t46 ^ 2); m(5) * (t65 * t10 - t66 * t11) + m(4) * (t65 * t12 - t66 * t13) + m(3) * (t65 * t18 - t66 * t19); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t90; (t53 * (-Icges(4,6) * t66 + t72 * t65) + t52 * (-Icges(4,5) * t66 + t74 * t65)) * t96 + (t53 * (Icges(4,6) * t65 + t72 * t66) + t52 * (Icges(4,5) * t65 + t74 * t66)) * t97 + m(5) * (t15 * t10 + t14 * t11) + m(4) * (-t12 * t66 - t13 * t65) * t42 + (t60 / 0.2e1 + t61 / 0.2e1) * (Icges(4,5) * t52 + Icges(4,6) * t53) + t84; m(5) * (-t14 * t66 + t15 * t65); m(4) * (t90 * t42 ^ 2 + t9 ^ 2) + t65 * (-t26 * t98 + t60 * t27) - t66 * (t61 * t26 - t27 * t98) + m(5) * (t14 ^ 2 + t15 ^ 2 + t3 ^ 2) + t83; m(5) * (-t10 * t66 - t11 * t65) * t36 + t84; 0; m(5) * (t8 * t3 + (-t14 * t65 - t15 * t66) * t36) + t83; m(5) * (t90 * t36 ^ 2 + t8 ^ 2) + t83;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
