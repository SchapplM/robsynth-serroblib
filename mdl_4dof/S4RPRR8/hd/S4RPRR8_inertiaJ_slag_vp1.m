% Calculate joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:07
% DurationCPUTime: 0.50s
% Computational Cost: add. (653->119), mult. (960->184), div. (0->0), fcn. (830->6), ass. (0->65)
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t100 = t63 * t61;
t59 = qJ(3) + qJ(4);
t50 = sin(t59);
t51 = cos(t59);
t101 = rSges(5,1) * t50 + rSges(5,2) * t51;
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t93 = rSges(4,2) * t62;
t99 = (rSges(4,1) * t60 + t93) * t63;
t57 = t61 ^ 2;
t58 = t63 ^ 2;
t98 = t61 / 0.2e1;
t97 = t63 / 0.2e1;
t67 = Icges(5,5) * t50 + Icges(5,6) * t51;
t17 = Icges(5,3) * t63 + t67 * t61;
t18 = Icges(5,3) * t61 - t67 * t63;
t96 = t61 * (t17 * t100 + t57 * t18) + t63 * (t18 * t100 + t58 * t17);
t95 = pkin(3) * t60;
t91 = t60 * t61;
t23 = t63 * rSges(5,3) + t101 * t61;
t90 = -pkin(3) * t91 - t23;
t89 = t63 * pkin(1) + t61 * qJ(2);
t88 = t57 + t58;
t87 = Icges(4,4) * t60;
t86 = Icges(4,4) * t62;
t85 = Icges(5,4) * t50;
t84 = Icges(5,4) * t51;
t83 = rSges(4,1) * t91 + t63 * rSges(4,3) + t61 * t93;
t33 = Icges(5,5) * t51 - Icges(5,6) * t50;
t69 = Icges(5,2) * t51 + t85;
t71 = Icges(5,1) * t50 + t84;
t34 = -Icges(5,2) * t50 + t84;
t35 = Icges(5,1) * t51 - t85;
t73 = t34 * t51 + t35 * t50;
t82 = (-t50 * (Icges(5,6) * t61 - t69 * t63) + t51 * (Icges(5,5) * t61 - t71 * t63) + t61 * t33 - t73 * t63) * t98 + (-t50 * (Icges(5,6) * t63 + t69 * t61) + t51 * (Icges(5,5) * t63 + t71 * t61) + t63 * t33 + t73 * t61) * t97;
t36 = t51 * rSges(5,1) - t50 * rSges(5,2);
t81 = pkin(3) * t62 + t36;
t15 = t81 * t61;
t16 = t81 * t63;
t78 = t15 * t61 + t16 * t63;
t72 = Icges(4,1) * t60 + t86;
t70 = Icges(4,2) * t62 + t87;
t68 = Icges(4,5) * t60 + Icges(4,6) * t62;
t53 = t63 * qJ(2);
t12 = t53 + t99 + (-rSges(4,3) - pkin(1) - pkin(5)) * t61;
t13 = t63 * pkin(5) + t83 + t89;
t66 = m(4) * (t61 * t12 - t63 * t13);
t64 = -pkin(6) - pkin(5);
t10 = t53 + (t101 + t95) * t63 + (-rSges(5,3) - pkin(1) + t64) * t61;
t11 = -t63 * t64 + t89 - t90;
t65 = m(5) * (t61 * t10 - t63 * t11);
t43 = t63 * rSges(2,1) - t61 * rSges(2,2);
t42 = t62 * rSges(4,1) - t60 * rSges(4,2);
t41 = -t61 * rSges(2,1) - t63 * rSges(2,2);
t31 = -t63 * rSges(3,2) + t61 * rSges(3,3) + t89;
t30 = t63 * rSges(3,3) + t53 + (rSges(3,2) - pkin(1)) * t61;
t25 = Icges(4,3) * t61 - t68 * t63;
t24 = Icges(4,3) * t63 + t68 * t61;
t14 = t63 * (t61 * rSges(5,3) - t101 * t63);
t9 = -t61 * t83 + (t61 * rSges(4,3) - t99) * t63;
t8 = -t61 * t23 + t14;
t3 = -t58 * t95 + t90 * t61 + t14;
t1 = [-t50 * t34 + t51 * t35 - t60 * (-Icges(4,2) * t60 + t86) + t62 * (Icges(4,1) * t62 - t87) + Icges(3,1) + Icges(2,3) + m(2) * (t41 ^ 2 + t43 ^ 2) + m(3) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2); m(3) * (t61 * t30 - t63 * t31) + t66 + t65; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t88; (-t60 * (Icges(4,6) * t63 + t70 * t61) + t62 * (Icges(4,5) * t63 + t72 * t61)) * t97 + (-t60 * (Icges(4,6) * t61 - t70 * t63) + t62 * (Icges(4,5) * t61 - t72 * t63)) * t98 + m(5) * (t15 * t10 - t16 * t11) + t42 * t66 + (t58 / 0.2e1 + t57 / 0.2e1) * (Icges(4,5) * t62 - Icges(4,6) * t60) + t82; m(4) * t88 * t42 + m(5) * t78; m(4) * (t88 * t42 ^ 2 + t9 ^ 2) + t63 * (t25 * t100 + t58 * t24) + t61 * (t24 * t100 + t57 * t25) + m(5) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + t96; t36 * t65 + t82; m(5) * t88 * t36; m(5) * (t8 * t3 + t78 * t36) + t96; m(5) * (t88 * t36 ^ 2 + t8 ^ 2) + t96;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
