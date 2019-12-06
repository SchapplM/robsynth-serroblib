% Calculate joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:07
% DurationCPUTime: 0.91s
% Computational Cost: add. (1268->184), mult. (1414->261), div. (0->0), fcn. (1224->8), ass. (0->93)
t88 = qJ(3) + pkin(8);
t78 = sin(t88);
t79 = cos(t88);
t92 = sin(qJ(3));
t94 = cos(qJ(3));
t151 = Icges(4,5) * t92 + Icges(5,5) * t78 + Icges(4,6) * t94 + Icges(5,6) * t79;
t150 = Icges(4,3) + Icges(5,3);
t93 = sin(qJ(1));
t95 = cos(qJ(1));
t149 = t93 * t95;
t80 = qJ(5) + t88;
t75 = sin(t80);
t76 = cos(t80);
t148 = rSges(6,1) * t75 + rSges(6,2) * t76;
t116 = rSges(5,1) * t78 + rSges(5,2) * t79;
t26 = t95 * rSges(6,3) + t148 * t93;
t139 = t92 * pkin(3);
t57 = pkin(4) * t78 + t139;
t147 = -t93 * t57 - t26;
t146 = t150 * t93 - t151 * t95;
t145 = t150 * t95 + t151 * t93;
t144 = (rSges(4,1) * t92 + rSges(4,2) * t94) * t95;
t89 = t93 ^ 2;
t90 = t95 ^ 2;
t143 = t93 / 0.2e1;
t142 = t95 / 0.2e1;
t98 = Icges(6,5) * t75 + Icges(6,6) * t76;
t20 = Icges(6,3) * t95 + t98 * t93;
t21 = Icges(6,3) * t93 - t98 * t95;
t141 = t93 * (t20 * t149 + t89 * t21) + t95 * (t21 * t149 + t90 * t20);
t69 = t89 + t90;
t55 = m(6) * t69;
t140 = pkin(3) * t94;
t134 = t92 * t93;
t133 = t93 * t94;
t84 = t95 * rSges(5,3);
t91 = -qJ(4) - pkin(6);
t132 = m(5) * t69 + t55;
t131 = t95 * t139 + t93 * t91;
t130 = t95 * pkin(1) + t93 * qJ(2);
t129 = Icges(4,4) * t92;
t128 = Icges(4,4) * t94;
t127 = Icges(5,4) * t78;
t126 = Icges(5,4) * t79;
t125 = Icges(6,4) * t75;
t124 = Icges(6,4) * t76;
t123 = -t116 * t93 - t84;
t122 = rSges(4,1) * t134 + rSges(4,2) * t133 + t95 * rSges(4,3);
t101 = Icges(6,2) * t76 + t125;
t104 = Icges(6,1) * t75 + t124;
t47 = -Icges(6,2) * t75 + t124;
t48 = Icges(6,1) * t76 - t125;
t107 = t47 * t76 + t48 * t75;
t46 = Icges(6,5) * t76 - Icges(6,6) * t75;
t121 = (-t107 * t95 - t75 * (Icges(6,6) * t93 - t101 * t95) + t76 * (Icges(6,5) * t93 - t104 * t95) + t93 * t46) * t143 + (t107 * t93 - t75 * (Icges(6,6) * t95 + t101 * t93) + t76 * (Icges(6,5) * t95 + t104 * t93) + t95 * t46) * t142;
t49 = t76 * rSges(6,1) - t75 * rSges(6,2);
t120 = pkin(4) * t79 + t49;
t72 = pkin(3) * t134;
t119 = -t95 * t91 + t72;
t73 = pkin(3) * t133;
t15 = t120 * t93 + t73;
t16 = (-t120 - t140) * t95;
t114 = t15 * t93 - t16 * t95;
t106 = Icges(4,1) * t92 + t128;
t105 = Icges(5,1) * t78 + t126;
t103 = Icges(4,2) * t94 + t129;
t102 = Icges(5,2) * t79 + t127;
t82 = t95 * qJ(2);
t17 = t82 + t144 + (-rSges(4,3) - pkin(1) - pkin(6)) * t93;
t18 = t95 * pkin(6) + t122 + t130;
t97 = m(4) * (t93 * t17 - t95 * t18);
t87 = -pkin(7) + t91;
t10 = t82 + (t148 + t57) * t95 + (-rSges(6,3) - pkin(1) + t87) * t93;
t11 = -t95 * t87 + t130 - t147;
t96 = m(6) * (t93 * t10 - t95 * t11);
t66 = t95 * rSges(2,1) - t93 * rSges(2,2);
t65 = t94 * rSges(4,1) - t92 * rSges(4,2);
t64 = -t93 * rSges(2,1) - t95 * rSges(2,2);
t54 = t79 * rSges(5,1) - t78 * rSges(5,2);
t44 = t72 + (-pkin(6) - t91) * t95;
t43 = -t95 * rSges(3,2) + t93 * rSges(3,3) + t130;
t42 = t95 * rSges(3,3) + t82 + (rSges(3,2) - pkin(1)) * t93;
t35 = t95 * (-t93 * pkin(6) - t131);
t28 = (-t54 - t140) * t95;
t27 = t93 * t54 + t73;
t19 = t95 * (t93 * rSges(6,3) - t148 * t95);
t14 = t119 - t123 + t130;
t13 = t82 + t116 * t95 + (-rSges(5,3) - pkin(1)) * t93 + t131;
t12 = -t93 * t122 + (t93 * rSges(4,3) - t144) * t95;
t9 = -t93 * t26 + t19;
t4 = t35 - t116 * t90 + (t123 - t44 + t84) * t93;
t3 = t35 + t95 * (-t95 * t57 + t131) + t19 + (t119 - t44 + t147) * t93;
t1 = [-t75 * t47 + t76 * t48 - t78 * (-Icges(5,2) * t78 + t126) + t79 * (Icges(5,1) * t79 - t127) - t92 * (-Icges(4,2) * t92 + t128) + t94 * (Icges(4,1) * t94 - t129) + Icges(3,1) + Icges(2,3) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(2) * (t64 ^ 2 + t66 ^ 2); m(5) * (t93 * t13 - t95 * t14) + t96 + t97 + m(3) * (t93 * t42 - t95 * t43); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t69 + t132; m(5) * (t27 * t13 + t28 * t14) + m(6) * (t15 * t10 + t16 * t11) + t65 * t97 + t121 + (-t78 * (Icges(5,6) * t93 - t102 * t95) + t79 * (Icges(5,5) * t93 - t105 * t95) - t92 * (Icges(4,6) * t93 - t103 * t95) + t94 * (Icges(4,5) * t93 - t106 * t95)) * t143 + (-t78 * (Icges(5,6) * t95 + t102 * t93) + t79 * (Icges(5,5) * t95 + t105 * t93) - t92 * (Icges(4,6) * t95 + t103 * t93) + t94 * (Icges(4,5) * t95 + t106 * t93)) * t142 + (Icges(4,5) * t94 + Icges(5,5) * t79 - Icges(4,6) * t92 - Icges(5,6) * t78) * (t90 / 0.2e1 + t89 / 0.2e1); m(5) * (t27 * t93 - t28 * t95) + m(6) * t114 + m(4) * t69 * t65; m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2 + t4 ^ 2) + m(4) * (t69 * t65 ^ 2 + t12 ^ 2) + t141 + t145 * t95 * t90 + (t146 * t89 + (t145 * t93 + t146 * t95) * t95) * t93; m(5) * (t95 * t13 + t93 * t14) + m(6) * (t95 * t10 + t93 * t11); 0; m(6) * (t95 * t15 + t93 * t16) + m(5) * (t95 * t27 + t93 * t28); t132; t49 * t96 + t121; t49 * t55; m(6) * (t114 * t49 + t9 * t3) + t141; 0; m(6) * (t69 * t49 ^ 2 + t9 ^ 2) + t141;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
