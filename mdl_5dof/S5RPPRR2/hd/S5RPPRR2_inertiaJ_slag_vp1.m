% Calculate joint inertia matrix for
% S5RPPRR2
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:37
% DurationCPUTime: 0.61s
% Computational Cost: add. (1144->145), mult. (1158->214), div. (0->0), fcn. (986->8), ass. (0->78)
t75 = sin(qJ(1));
t76 = cos(qJ(1));
t116 = t76 * t75;
t69 = pkin(8) + qJ(4);
t62 = qJ(5) + t69;
t57 = sin(t62);
t58 = cos(t62);
t118 = rSges(6,1) * t57 + rSges(6,2) * t58;
t60 = sin(t69);
t61 = cos(t69);
t117 = rSges(5,1) * t60 + rSges(5,2) * t61;
t115 = t117 * t76;
t25 = t76 * rSges(6,3) + t118 * t75;
t72 = sin(pkin(8));
t110 = t72 * pkin(3);
t48 = pkin(4) * t60 + t110;
t114 = -t75 * t48 - t25;
t70 = t75 ^ 2;
t71 = t76 ^ 2;
t113 = t75 / 0.2e1;
t112 = t76 / 0.2e1;
t80 = Icges(6,5) * t57 + Icges(6,6) * t58;
t19 = Icges(6,3) * t76 + t80 * t75;
t20 = Icges(6,3) * t75 - t80 * t76;
t111 = t75 * (t19 * t116 + t70 * t20) + t76 * (t20 * t116 + t71 * t19);
t55 = t70 + t71;
t46 = m(5) * t55;
t45 = m(6) * t55;
t74 = -pkin(6) - qJ(3);
t105 = t76 * t110 + t75 * t74;
t104 = t76 * pkin(1) + t75 * qJ(2);
t103 = Icges(5,4) * t60;
t102 = Icges(5,4) * t61;
t101 = Icges(6,4) * t57;
t100 = Icges(6,4) * t58;
t99 = rSges(4,3) + qJ(3);
t98 = m(4) * t55 + t45 + t46;
t97 = t76 * rSges(5,3) + t117 * t75;
t35 = Icges(6,5) * t58 - Icges(6,6) * t57;
t82 = Icges(6,2) * t58 + t101;
t84 = Icges(6,1) * t57 + t100;
t36 = -Icges(6,2) * t57 + t100;
t37 = Icges(6,1) * t58 - t101;
t86 = t36 * t58 + t37 * t57;
t96 = (-t57 * (Icges(6,6) * t75 - t82 * t76) + t58 * (Icges(6,5) * t75 - t84 * t76) + t75 * t35 - t86 * t76) * t113 + (-t57 * (Icges(6,6) * t76 + t82 * t75) + t58 * (Icges(6,5) * t76 + t84 * t75) + t76 * t35 + t86 * t75) * t112;
t38 = t58 * rSges(6,1) - t57 * rSges(6,2);
t95 = pkin(4) * t61 + t38;
t73 = cos(pkin(8));
t94 = rSges(4,1) * t72 + rSges(4,2) * t73;
t16 = t95 * t75;
t17 = t95 * t76;
t91 = t16 * t75 + t17 * t76;
t85 = Icges(5,1) * t60 + t102;
t83 = Icges(5,2) * t61 + t103;
t81 = Icges(5,5) * t60 + Icges(5,6) * t61;
t79 = t75 * t110 - t76 * t74;
t68 = -pkin(7) + t74;
t10 = -t76 * t68 + t104 - t114;
t64 = t76 * qJ(2);
t9 = t64 + (t48 + t118) * t76 + (-rSges(6,3) - pkin(1) + t68) * t75;
t78 = m(6) * (-t76 * t10 + t75 * t9);
t12 = t64 + t115 + (-rSges(5,3) - pkin(1)) * t75 + t105;
t13 = t79 + t97 + t104;
t77 = m(5) * (t75 * t12 - t76 * t13);
t52 = t76 * rSges(2,1) - t75 * rSges(2,2);
t51 = -t75 * rSges(2,1) - t76 * rSges(2,2);
t44 = t61 * rSges(5,1) - t60 * rSges(5,2);
t33 = -t76 * rSges(3,2) + t75 * rSges(3,3) + t104;
t32 = t76 * rSges(3,3) + t64 + (rSges(3,2) - pkin(1)) * t75;
t27 = Icges(5,3) * t75 - t81 * t76;
t26 = Icges(5,3) * t76 + t81 * t75;
t18 = t76 * (t75 * rSges(6,3) - t118 * t76);
t15 = t94 * t75 + t99 * t76 + t104;
t14 = t64 + t94 * t76 + (-pkin(1) - t99) * t75;
t11 = -t75 * t97 + (t75 * rSges(5,3) - t115) * t76;
t8 = -t75 * t25 + t18;
t3 = t76 * (-t76 * t48 + t105) + t18 + (t79 + t114) * t75;
t1 = [Icges(4,1) * t73 ^ 2 - t57 * t36 + t58 * t37 - t60 * (-Icges(5,2) * t60 + t102) + t61 * (Icges(5,1) * t61 - t103) + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t73 + Icges(4,2) * t72) * t72 + m(6) * (t10 ^ 2 + t9 ^ 2) + m(2) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2); t78 + m(3) * (t75 * t32 - t76 * t33) + m(4) * (t75 * t14 - t76 * t15) + t77; m(3) * t55 + t98; m(6) * (t75 * t10 + t76 * t9) + m(4) * (t76 * t14 + t75 * t15) + m(5) * (t76 * t12 + t75 * t13); 0; t98; (-t60 * (Icges(5,6) * t76 + t83 * t75) + t61 * (Icges(5,5) * t76 + t85 * t75)) * t112 + (-t60 * (Icges(5,6) * t75 - t83 * t76) + t61 * (Icges(5,5) * t75 - t85 * t76)) * t113 + m(6) * (-t17 * t10 + t16 * t9) + t44 * t77 + (t71 / 0.2e1 + t70 / 0.2e1) * (Icges(5,5) * t61 - Icges(5,6) * t60) + t96; m(6) * t91 + t44 * t46; m(6) * (t16 * t76 - t17 * t75); m(5) * (t55 * t44 ^ 2 + t11 ^ 2) + t76 * (t27 * t116 + t71 * t26) + t75 * (t26 * t116 + t70 * t27) + m(6) * (t16 ^ 2 + t17 ^ 2 + t3 ^ 2) + t111; t38 * t78 + t96; t38 * t45; 0; m(6) * (t8 * t3 + t91 * t38) + t111; m(6) * (t55 * t38 ^ 2 + t8 ^ 2) + t111;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
