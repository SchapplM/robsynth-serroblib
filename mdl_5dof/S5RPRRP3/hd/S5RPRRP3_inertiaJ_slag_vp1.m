% Calculate joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:55
% EndTime: 2020-01-03 11:47:01
% DurationCPUTime: 1.68s
% Computational Cost: add. (2263->177), mult. (1730->238), div. (0->0), fcn. (1546->8), ass. (0->100)
t180 = Icges(5,4) + Icges(6,4);
t179 = Icges(5,1) + Icges(6,1);
t178 = Icges(5,2) + Icges(6,2);
t100 = qJ(3) + qJ(4);
t94 = cos(t100);
t177 = t180 * t94;
t93 = sin(t100);
t176 = t180 * t93;
t175 = Icges(5,5) + Icges(6,5);
t174 = -Icges(5,6) - Icges(6,6);
t173 = -t178 * t93 + t177;
t172 = t179 * t94 - t176;
t99 = qJ(1) + pkin(8);
t91 = sin(t99);
t92 = cos(t99);
t171 = t173 * t91 + t174 * t92;
t170 = -t173 * t92 + t174 * t91;
t169 = t172 * t91 - t175 * t92;
t168 = -t172 * t92 - t175 * t91;
t167 = Icges(5,3) + Icges(6,3);
t166 = t174 * t93 + t175 * t94;
t105 = -pkin(7) - pkin(6);
t159 = rSges(6,3) + qJ(5) - t105;
t165 = t178 * t94 + t176;
t164 = t179 * t93 + t177;
t163 = -t166 * t91 + t167 * t92;
t162 = t166 * t92 + t167 * t91;
t161 = t168 * t94 - t170 * t93;
t160 = -t169 * t94 + t171 * t93;
t158 = t174 * t94 - t175 * t93;
t89 = t91 ^ 2;
t90 = t92 ^ 2;
t139 = t90 + t89;
t157 = t94 * rSges(6,2) + (rSges(6,1) + pkin(4)) * t93;
t145 = t91 * t94;
t146 = t91 * t93;
t103 = cos(qJ(3));
t88 = t103 * pkin(3) + pkin(2);
t69 = pkin(4) * t94 + t88;
t156 = rSges(6,1) * t145 - rSges(6,2) * t146 - t159 * t92 + t91 * t69;
t155 = -t164 * t94 + t165 * t93;
t154 = (rSges(6,1) * t94 - rSges(6,2) * t93 + t69) * t92;
t153 = t163 * t90 + (t161 * t91 + (-t160 + t162) * t92) * t91;
t152 = t162 * t89 + (t160 * t92 + (-t161 + t163) * t91) * t92;
t151 = (rSges(5,1) * t94 - rSges(5,2) * t93) * t92;
t150 = -t91 / 0.2e1;
t149 = -t92 / 0.2e1;
t101 = sin(qJ(3));
t147 = pkin(3) * t101;
t141 = t92 * t105 + t91 * t88;
t144 = (-t141 + t156) * t91;
t71 = t92 * t88;
t143 = -t71 + t154 + (t105 + t159) * t91;
t33 = t157 * t92;
t140 = t92 * pkin(2) + t91 * pkin(6);
t138 = rSges(4,1) * t103;
t137 = rSges(4,2) * t101;
t132 = Icges(4,4) * t101;
t131 = Icges(4,4) * t103;
t121 = -t137 + t138;
t112 = Icges(4,1) * t103 - t132;
t111 = -Icges(4,2) * t101 + t131;
t110 = Icges(4,5) * t103 - Icges(4,6) * t101;
t109 = rSges(5,1) * t145 - rSges(5,2) * t146 - t92 * rSges(5,3);
t107 = (t155 * t92 + t158 * t91 + t168 * t93 + t170 * t94) * t150 + (-t155 * t91 + t158 * t92 + t169 * t93 + t171 * t94) * t149;
t106 = t152 * t91 + t153 * t92;
t104 = cos(qJ(1));
t102 = sin(qJ(1));
t97 = t104 * pkin(1);
t95 = t102 * pkin(1);
t84 = t92 * t147;
t82 = t91 * t138;
t81 = t104 * rSges(2,1) - t102 * rSges(2,2);
t80 = t102 * rSges(2,1) + t104 * rSges(2,2);
t79 = t101 * rSges(4,1) + t103 * rSges(4,2);
t68 = t93 * rSges(5,1) + t94 * rSges(5,2);
t58 = t92 * rSges(3,1) - t91 * rSges(3,2) + t97;
t57 = t91 * rSges(3,1) + t92 * rSges(3,2) + t95;
t50 = -Icges(4,3) * t92 + t110 * t91;
t49 = t92 * t68 + t84;
t48 = (-t68 - t147) * t91;
t47 = -t91 * rSges(5,3) - t151;
t32 = t157 * t91;
t31 = t91 * t105 + t140 - t71;
t30 = t91 * t109;
t28 = t91 * (-t91 * pkin(2) + t92 * pkin(6) + t141);
t27 = t84 + t33;
t26 = (-t157 - t147) * t91;
t25 = t91 * rSges(4,3) + t121 * t92 + t140 + t97;
t24 = t82 + t95 + (-rSges(4,3) - pkin(6)) * t92 + (pkin(2) - t137) * t91;
t21 = t71 + t97 + t151 + (rSges(5,3) - t105) * t91;
t20 = t109 + t95 + t141;
t19 = t159 * t91 + t154 + t97;
t18 = t95 + t156;
t17 = t91 * (-t91 * t137 + t82) + t121 * t90;
t16 = -t92 * t47 + t30;
t7 = t28 + t30 + (-t31 - t47) * t92;
t6 = t143 * t92 + t144;
t1 = t28 + (-t31 + t143) * t92 + t144;
t2 = [t101 * (Icges(4,1) * t101 + t131) + t103 * (Icges(4,2) * t103 + t132) + Icges(2,3) + Icges(3,3) + t165 * t94 + t164 * t93 + m(2) * (t80 ^ 2 + t81 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(6) * (t18 ^ 2 + t19 ^ 2); 0; m(3) + m(4) + m(5) + m(6); (t101 * (-Icges(4,5) * t92 + t112 * t91) + t103 * (-Icges(4,6) * t92 + t111 * t91)) * t149 + (t101 * (-Icges(4,5) * t91 - t112 * t92) + t103 * (-Icges(4,6) * t91 - t111 * t92)) * t150 + m(5) * (t49 * t20 + t48 * t21) + m(6) * (t27 * t18 + t26 * t19) + m(4) * (t24 * t92 - t25 * t91) * t79 + (t90 / 0.2e1 + t89 / 0.2e1) * (Icges(4,5) * t101 + Icges(4,6) * t103) + t107; m(4) * t17 + m(5) * t7 + m(6) * t1; m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2 + t7 ^ 2) + m(4) * (t139 * t79 ^ 2 + t17 ^ 2) + (-t90 * t50 + t153) * t92 + (-t91 * t50 * t92 - t139 * (-Icges(4,3) * t91 - t110 * t92) + t152) * t91; m(6) * (t33 * t18 - t32 * t19) + m(5) * (t20 * t92 - t21 * t91) * t68 + t107; m(5) * t16 + m(6) * t6; m(6) * (t6 * t1 - t32 * t26 + t33 * t27) + m(5) * (t16 * t7 + (-t48 * t91 + t49 * t92) * t68) + t106; m(5) * (t139 * t68 ^ 2 + t16 ^ 2) + m(6) * (t32 ^ 2 + t33 ^ 2 + t6 ^ 2) + t106; m(6) * (-t91 * t18 - t92 * t19); 0; m(6) * (-t92 * t26 - t91 * t27); m(6) * (t92 * t32 - t91 * t33); m(6) * t139;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
