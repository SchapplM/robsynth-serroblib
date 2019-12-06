% Calculate joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:15
% EndTime: 2019-12-05 17:59:20
% DurationCPUTime: 1.37s
% Computational Cost: add. (1352->177), mult. (1784->244), div. (0->0), fcn. (1572->6), ass. (0->93)
t177 = Icges(5,4) + Icges(6,4);
t176 = Icges(5,1) + Icges(6,1);
t175 = Icges(5,5) + Icges(6,5);
t174 = Icges(5,2) + Icges(6,2);
t173 = Icges(5,6) + Icges(6,6);
t101 = qJ(3) + qJ(4);
t91 = cos(t101);
t172 = t177 * t91;
t90 = sin(t101);
t171 = t177 * t90;
t168 = t176 * t91 - t171;
t169 = -t174 * t90 + t172;
t170 = (t169 + t172) * t91 + (t168 - t171 + (-t174 + t176) * t91) * t90;
t167 = Icges(5,3) + Icges(6,3);
t166 = t173 * t91 + t175 * t90;
t163 = -0.2e1 * t173 * t90 + 0.2e1 * t175 * t91;
t105 = cos(qJ(1));
t102 = sin(qJ(3));
t148 = t102 * pkin(3);
t70 = pkin(4) * t90 + t148;
t162 = (-rSges(6,1) * t90 - rSges(6,2) * t91 - t70) * t105;
t106 = -pkin(7) - pkin(6);
t98 = -qJ(5) + t106;
t161 = rSges(6,3) - t98;
t103 = sin(qJ(1));
t160 = t166 * t103 + t167 * t105;
t159 = t167 * t103 - t166 * t105;
t138 = t103 * t91;
t139 = t103 * t90;
t155 = -rSges(6,1) * t139 - rSges(6,2) * t138 - t105 * rSges(6,3) - t103 * t70;
t154 = t103 * t105;
t104 = cos(qJ(3));
t153 = (rSges(4,1) * t102 + rSges(4,2) * t104) * t105;
t152 = (rSges(5,1) * t90 + rSges(5,2) * t91) * t105;
t99 = t103 ^ 2;
t100 = t105 ^ 2;
t151 = t103 / 0.2e1;
t150 = t105 / 0.2e1;
t149 = pkin(3) * t104;
t145 = t103 * t106 + t105 * t148;
t147 = (t161 * t103 + t145 + t162) * t105;
t135 = t102 * t103;
t86 = pkin(3) * t135;
t146 = t86 - (t106 - t98) * t105 + t155;
t67 = t91 * rSges(6,1) - t90 * rSges(6,2);
t30 = pkin(4) * t138 + t103 * t67;
t144 = t105 * pkin(1) + t103 * qJ(2);
t83 = t99 + t100;
t137 = Icges(4,4) * t102;
t136 = Icges(4,4) * t104;
t134 = t103 * t104;
t47 = rSges(5,1) * t139 + rSges(5,2) * t138 + t105 * rSges(5,3);
t133 = rSges(4,1) * t135 + rSges(4,2) * t134 + t105 * rSges(4,3);
t132 = t159 * t99 * t103 + (t160 * t100 + (t160 * t103 + t159 * t105) * t103) * t105;
t131 = -pkin(4) * t91 - t67;
t68 = t91 * rSges(5,1) - t90 * rSges(5,2);
t87 = pkin(3) * t134;
t32 = t103 * t68 + t87;
t33 = (-t68 - t149) * t105;
t113 = t32 * t103 - t33 * t105;
t112 = Icges(4,1) * t102 + t136;
t111 = Icges(4,2) * t104 + t137;
t110 = Icges(4,5) * t102 + Icges(4,6) * t104;
t93 = t105 * qJ(2);
t26 = t93 + t153 + (-rSges(4,3) - pkin(1) - pkin(6)) * t103;
t27 = t105 * pkin(6) + t133 + t144;
t109 = m(4) * (t103 * t26 - t105 * t27);
t20 = t93 + t152 + (-rSges(5,3) - pkin(1)) * t103 + t145;
t21 = -t105 * t106 + t144 + t47 + t86;
t108 = m(5) * (t103 * t20 - t105 * t21);
t107 = (t163 * t150 - t170 * t151) * t105 + (t170 * t150 + t163 * t151) * t103;
t77 = t105 * rSges(2,1) - t103 * rSges(2,2);
t76 = t104 * rSges(4,1) - t102 * rSges(4,2);
t75 = -t103 * rSges(2,1) - t105 * rSges(2,2);
t69 = m(6) * t83;
t57 = t86 + (-pkin(6) - t106) * t105;
t56 = -t105 * rSges(3,2) + t103 * rSges(3,3) + t144;
t55 = t105 * rSges(3,3) + t93 + (rSges(3,2) - pkin(1)) * t103;
t50 = Icges(4,3) * t103 - t110 * t105;
t49 = Icges(4,3) * t105 + t110 * t103;
t48 = t105 * (-t103 * pkin(6) - t145);
t31 = t131 * t105;
t29 = t105 * (t103 * rSges(5,3) - t152);
t25 = (t131 - t149) * t105;
t24 = t87 + t30;
t19 = -t103 * t133 + (t103 * rSges(4,3) - t153) * t105;
t18 = -t105 * t98 + t144 - t155;
t17 = t93 - t162 + (-pkin(1) - t161) * t103;
t16 = -t103 * t47 + t29;
t7 = t29 + t48 + (-t47 - t57) * t103;
t6 = t146 * t103 + t147;
t5 = t48 + (-t57 + t146) * t103 + t147;
t1 = [-t102 * (-Icges(4,2) * t102 + t136) + t104 * (Icges(4,1) * t104 - t137) + Icges(3,1) + Icges(2,3) + t168 * t91 - t169 * t90 + m(6) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t75 ^ 2 + t77 ^ 2); m(6) * (t103 * t17 - t105 * t18) + t108 + t109 + m(3) * (t103 * t55 - t105 * t56); t69 + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t83; (-t102 * (Icges(4,6) * t105 + t111 * t103) + t104 * (Icges(4,5) * t105 + t112 * t103)) * t150 + (-t102 * (Icges(4,6) * t103 - t111 * t105) + t104 * (Icges(4,5) * t103 - t112 * t105)) * t151 + m(6) * (t24 * t17 + t25 * t18) + m(5) * (t32 * t20 + t33 * t21) + t76 * t109 + (t99 / 0.2e1 + t100 / 0.2e1) * (Icges(4,5) * t104 - Icges(4,6) * t102) + t107; m(5) * t113 + m(6) * (t24 * t103 - t25 * t105) + m(4) * t83 * t76; m(6) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t7 ^ 2) + m(4) * (t83 * t76 ^ 2 + t19 ^ 2) + t105 * (t100 * t49 + t50 * t154) + t103 * (t49 * t154 + t99 * t50) + t132; m(6) * (t30 * t17 + t31 * t18) + t68 * t108 + t107; m(6) * (t30 * t103 - t31 * t105) + m(5) * t83 * t68; m(6) * (t30 * t24 + t31 * t25 + t6 * t5) + m(5) * (t113 * t68 + t16 * t7) + t132; m(5) * (t83 * t68 ^ 2 + t16 ^ 2) + m(6) * (t30 ^ 2 + t31 ^ 2 + t6 ^ 2) + t132; m(6) * (t103 * t18 + t105 * t17); 0; m(6) * (t103 * t25 + t105 * t24); m(6) * (t103 * t31 + t105 * t30); t69;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
