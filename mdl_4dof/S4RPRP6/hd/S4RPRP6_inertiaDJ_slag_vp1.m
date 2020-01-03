% Calculate time derivative of joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:06
% DurationCPUTime: 4.78s
% Computational Cost: add. (1027->236), mult. (2628->337), div. (0->0), fcn. (1972->4), ass. (0->124)
t99 = cos(qJ(3));
t165 = Icges(5,4) * t99;
t167 = Icges(4,4) * t99;
t97 = sin(qJ(3));
t206 = t165 + t167 + (-Icges(4,2) - Icges(5,2)) * t97;
t166 = Icges(5,4) * t97;
t168 = Icges(4,4) * t97;
t205 = -t166 - t168 + (Icges(4,1) + Icges(5,1)) * t99;
t100 = cos(qJ(1));
t150 = qJD(1) * t100;
t141 = t99 * t150;
t152 = qJD(3) * t99;
t98 = sin(qJ(1));
t198 = t97 * t150 + t98 * t152;
t145 = rSges(4,1) * t198 + rSges(4,2) * t141;
t148 = -rSges(4,3) - pkin(1) - pkin(5);
t154 = qJD(3) * t97;
t171 = qJ(2) * t150 + qJD(2) * t98;
t11 = (-rSges(4,2) * t154 + qJD(1) * t148) * t98 + t145 + t171;
t135 = rSges(4,1) * t97 + rSges(4,2) * t99;
t149 = qJD(3) * t100;
t180 = rSges(4,2) * t97;
t70 = rSges(4,1) * t99 - t180;
t88 = qJD(2) * t100;
t12 = t88 + t70 * t149 + ((-qJ(2) - t135) * t98 + t148 * t100) * qJD(1);
t108 = t100 * t135;
t90 = t100 * qJ(2);
t29 = t148 * t98 + t108 + t90;
t170 = t100 * pkin(1) + t98 * qJ(2);
t174 = t98 * t99;
t177 = t97 * t98;
t92 = t100 * rSges(4,3);
t48 = rSges(4,1) * t177 + rSges(4,2) * t174 + t92;
t30 = pkin(5) * t100 + t170 + t48;
t204 = -t100 * t11 + t98 * t12 + (t100 * t29 + t30 * t98) * qJD(1);
t201 = t206 * qJD(3);
t200 = t205 * qJD(3);
t183 = rSges(5,1) + pkin(3);
t199 = t183 * t99;
t197 = -rSges(5,2) * t174 - t100 * rSges(5,3) - t177 * t183;
t120 = Icges(4,2) * t99 + t168;
t41 = Icges(4,6) * t100 + t120 * t98;
t124 = Icges(4,1) * t97 + t167;
t45 = Icges(4,5) * t100 + t124 * t98;
t127 = t41 * t99 + t45 * t97;
t106 = t127 * t100;
t118 = Icges(5,2) * t99 + t166;
t39 = Icges(5,6) * t100 + t118 * t98;
t122 = Icges(5,1) * t97 + t165;
t43 = Icges(5,5) * t100 + t122 * t98;
t129 = t39 * t99 + t43 * t97;
t107 = t129 * t100;
t155 = qJD(1) * t98;
t96 = -qJ(4) - pkin(5);
t196 = rSges(5,2) * t141 + qJD(4) * t100 + t96 * t155 + t183 * t198;
t178 = rSges(5,2) * t99;
t112 = t183 * t97 + t178;
t114 = Icges(5,5) * t97 + Icges(5,6) * t99;
t194 = -Icges(5,3) * t98 + t100 * t114;
t116 = Icges(4,5) * t97 + Icges(4,6) * t99;
t193 = -Icges(4,3) * t98 + t100 * t116;
t192 = -Icges(5,6) * t98 + t100 * t118;
t191 = -Icges(4,6) * t98 + t100 * t120;
t190 = -Icges(5,5) * t98 + t100 * t122;
t189 = -Icges(4,5) * t98 + t100 * t124;
t188 = 2 * m(4);
t187 = 2 * m(5);
t94 = t98 ^ 2;
t95 = t100 ^ 2;
t182 = rSges(3,2) - pkin(1);
t181 = -pkin(5) - t96;
t179 = rSges(5,2) * t97;
t176 = t98 * rSges(4,3);
t173 = -t100 * t181 + t197;
t134 = -rSges(5,1) * t97 - t178;
t172 = (rSges(5,3) + t181) * t98 + (-pkin(3) * t97 + t134) * t100;
t169 = t96 - rSges(5,3);
t35 = Icges(5,3) * t100 + t114 * t98;
t157 = qJD(1) * t35;
t37 = Icges(4,3) * t100 + t116 * t98;
t156 = qJD(1) * t37;
t153 = qJD(3) * t98;
t151 = qJD(4) * t98;
t147 = -pkin(1) + t169;
t144 = t97 * t153;
t140 = t97 * t149;
t69 = rSges(5,1) * t99 - t179;
t138 = pkin(3) * t99 + t69;
t62 = t135 * qJD(3);
t137 = t62 * (t94 + t95);
t128 = -t190 * t97 - t192 * t99;
t126 = -t189 * t97 - t191 * t99;
t117 = Icges(4,5) * t99 - Icges(4,6) * t97;
t115 = Icges(5,5) * t99 - Icges(5,6) * t97;
t111 = t128 * t98;
t110 = t126 * t98;
t109 = rSges(3,3) * t100 + t182 * t98;
t61 = t134 * qJD(3);
t52 = -rSges(3,2) * t100 + t98 * rSges(3,3) + t170;
t51 = t90 + t109;
t50 = t176 - t108;
t34 = t138 * t100;
t33 = t138 * t98;
t32 = t88 + ((-rSges(3,3) - qJ(2)) * t98 + t182 * t100) * qJD(1);
t31 = qJD(1) * t109 + t171;
t28 = -t100 * t96 + t170 - t197;
t27 = t112 * t100 + t147 * t98 + t90;
t18 = qJD(1) * t193 + t117 * t153;
t17 = -t117 * t149 + t156;
t16 = qJD(1) * t194 + t115 * t153;
t15 = -t115 * t149 + t157;
t14 = t69 * t150 + t98 * t61 + (t141 - t144) * pkin(3);
t13 = t69 * t155 - t100 * t61 + (t155 * t99 + t140) * pkin(3);
t10 = -t100 * t126 - t193 * t98;
t9 = t98 * t37 - t106;
t8 = -t100 * t128 - t194 * t98;
t7 = t98 * t35 - t107;
t6 = -t100 * t193 + t110;
t5 = t100 * t37 + t127 * t98;
t4 = -t100 * t194 + t111;
t3 = t100 * t35 + t129 * t98;
t2 = -t151 + t88 + (-t179 + t199) * t149 + (t147 * t100 + (-qJ(2) - t112) * t98) * qJD(1);
t1 = (-rSges(5,2) * t154 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t98 + t171 + t196;
t19 = [0.2e1 * m(3) * (t31 * t52 + t32 * t51) + (t11 * t30 + t12 * t29) * t188 + (t1 * t28 + t2 * t27) * t187 - t201 * t99 - t200 * t97 + (t120 + t118) * t154 + (-t124 - t122) * t152; m(3) * (-t100 * t31 + t98 * t32 + (t100 * t51 + t52 * t98) * qJD(1)) + m(4) * t204 + m(5) * (-t1 * t100 + t98 * t2 + (t100 * t27 + t28 * t98) * qJD(1)); 0; m(5) * (-t1 * t34 + t13 * t28 + t14 * t27 + t2 * t33) + (-t106 / 0.2e1 - t107 / 0.2e1 + (-t126 / 0.2e1 - t128 / 0.2e1) * t98 + (-t116 - t114) * (t95 / 0.2e1 + t94 / 0.2e1)) * qJD(3) + ((-t205 * t99 + t206 * t97) * t149 + (t200 * t99 - t201 * t97) * t100) * t98 / 0.2e1 + (-(-t100 * t30 + t29 * t98) * t62 + t204 * t70) * m(4); m(5) * (-t13 * t100 + t14 * t98 + (t100 * t33 - t34 * t98) * qJD(1)) - m(4) * t137; ((t100 * t50 - t98 * t48) * (-t98 * t145 + (t180 * t94 - t70 * t95) * qJD(3) + ((-t48 + t92) * t100 + (t108 - t50 + t176) * t98) * qJD(1)) - t70 * t137) * t188 + t100 * ((t100 * t18 + (t6 + t106) * qJD(1)) * t100 + (-t5 * qJD(1) + (-t152 * t189 + t154 * t191) * t98 + (t17 + (-t41 * t97 + t45 * t99) * qJD(3) + (t126 - t37) * qJD(1)) * t100) * t98) + t98 * ((t98 * t17 + (-t9 + t110) * qJD(1)) * t98 + (t10 * qJD(1) + (-t152 * t45 + t154 * t41 + t156) * t100 + (t18 + (t189 * t99 - t191 * t97) * qJD(3) + t127 * qJD(1)) * t98) * t100) + (-t34 * t13 + t33 * t14 + ((rSges(5,2) * t144 - t196) * t98 + (rSges(5,2) * t140 - t149 * t199 + t151) * t100 + (((-pkin(5) + rSges(5,3)) * t98 - t172) * t98 + ((-pkin(5) - t169) * t100 + t112 * t98 + t173) * t100) * qJD(1)) * (t100 * t172 + t173 * t98)) * t187 + t100 * ((t100 * t16 + (t4 + t107) * qJD(1)) * t100 + (-t3 * qJD(1) + (-t152 * t190 + t154 * t192) * t98 + (t15 + (-t39 * t97 + t43 * t99) * qJD(3) + (t128 - t35) * qJD(1)) * t100) * t98) + t98 * ((t98 * t15 + (-t7 + t111) * qJD(1)) * t98 + (t8 * qJD(1) + (-t152 * t43 + t154 * t39 + t157) * t100 + (t16 + (t190 * t99 - t192 * t97) * qJD(3) + t129 * qJD(1)) * t98) * t100) + ((-t4 - t6) * t98 + (-t3 - t5) * t100) * t155 + ((t10 + t8) * t98 + (t7 + t9) * t100) * t150; m(5) * (t98 * t1 + t100 * t2 + (t100 * t28 - t27 * t98) * qJD(1)); 0; m(5) * (t100 * t14 + t98 * t13 + (-t100 * t34 - t33 * t98) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
