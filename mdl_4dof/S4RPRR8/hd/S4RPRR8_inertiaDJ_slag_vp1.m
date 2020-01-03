% Calculate time derivative of joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR8_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:10
% DurationCPUTime: 2.45s
% Computational Cost: add. (2794->286), mult. (4244->426), div. (0->0), fcn. (3260->6), ass. (0->168)
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t117 = qJ(3) + qJ(4);
t107 = sin(t117);
t108 = cos(t117);
t151 = rSges(5,1) * t107 + rSges(5,2) * t108;
t60 = t121 * rSges(5,3) + t119 * t151;
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t176 = qJD(3) * t120;
t179 = qJD(1) * t121;
t235 = t118 * t179 + t119 * t176;
t201 = Icges(4,4) * t118;
t139 = Icges(4,2) * t120 + t201;
t68 = Icges(4,6) * t121 + t119 * t139;
t200 = Icges(4,4) * t120;
t142 = Icges(4,1) * t118 + t200;
t70 = Icges(4,5) * t121 + t119 * t142;
t145 = t118 * t70 + t120 * t68;
t234 = t121 * t145;
t199 = Icges(5,4) * t107;
t138 = Icges(5,2) * t108 + t199;
t56 = Icges(5,6) * t121 + t119 * t138;
t198 = Icges(5,4) * t108;
t141 = Icges(5,1) * t107 + t198;
t58 = Icges(5,5) * t121 + t119 * t141;
t150 = t107 * t58 + t108 * t56;
t233 = t121 * t150;
t210 = rSges(4,2) * t120;
t152 = rSges(4,1) * t118 + t210;
t131 = t121 * t152;
t184 = t118 * t119;
t103 = pkin(3) * t184;
t122 = -pkin(6) - pkin(5);
t215 = -pkin(5) - t122;
t232 = t215 * t121 + t103;
t166 = t120 * t179;
t171 = t235 * rSges(4,1) + rSges(4,2) * t166;
t174 = -rSges(4,3) - pkin(1) - pkin(5);
t178 = qJD(3) * t118;
t183 = qJ(2) * t179 + qJD(2) * t119;
t26 = (-rSges(4,2) * t178 + qJD(1) * t174) * t119 + t171 + t183;
t106 = qJD(2) * t121;
t175 = qJD(3) * t121;
t211 = rSges(4,2) * t118;
t92 = rSges(4,1) * t120 - t211;
t27 = t106 + t92 * t175 + (t174 * t121 + (-qJ(2) - t152) * t119) * qJD(1);
t231 = t119 * t27 - t121 * t26;
t135 = Icges(5,5) * t107 + Icges(5,6) * t108;
t230 = -Icges(5,3) * t119 + t121 * t135;
t136 = Icges(4,5) * t118 + Icges(4,6) * t120;
t229 = -Icges(4,3) * t119 + t121 * t136;
t228 = -Icges(5,6) * t119 + t121 * t138;
t227 = -Icges(4,6) * t119 + t121 * t139;
t226 = -Icges(5,5) * t119 + t121 * t141;
t225 = -Icges(4,5) * t119 + t121 * t142;
t114 = qJD(3) + qJD(4);
t63 = t138 * t114;
t64 = t141 * t114;
t78 = Icges(5,5) * t108 - Icges(5,6) * t107;
t79 = -Icges(5,2) * t107 + t198;
t80 = Icges(5,1) * t108 - t199;
t224 = t107 * (t114 * t79 + t64) - t108 * (t114 * t80 - t63) + qJD(1) * t78;
t223 = 2 * m(4);
t222 = 2 * m(5);
t115 = t119 ^ 2;
t116 = t121 ^ 2;
t221 = m(4) * t92;
t220 = t119 / 0.2e1;
t219 = t121 / 0.2e1;
t218 = rSges(3,2) - pkin(1);
t217 = -rSges(5,3) - pkin(1);
t216 = pkin(3) * t118;
t214 = -t60 - t232;
t212 = rSges(5,1) * t108;
t208 = t118 * t68;
t207 = t118 * t227;
t206 = t119 * rSges(4,3);
t204 = t120 * t70;
t203 = t120 * t225;
t112 = t121 * rSges(4,3);
t54 = Icges(5,3) * t121 + t119 * t135;
t191 = qJD(1) * t54;
t66 = Icges(4,3) * t121 + t119 * t136;
t190 = qJD(1) * t66;
t188 = t107 * t114;
t187 = t108 * t114;
t186 = t114 * t119;
t185 = t114 * t121;
t182 = t121 * pkin(1) + t119 * qJ(2);
t181 = t115 + t116;
t180 = qJD(1) * t119;
t177 = qJD(3) * t119;
t173 = t122 + t217;
t172 = t151 * t179 + t186 * t212;
t170 = -t235 * pkin(3) - t122 * t180;
t169 = rSges(5,2) * t188;
t168 = pkin(3) * t176;
t72 = rSges(4,1) * t184 + t119 * t210 + t112;
t81 = -rSges(5,2) * t107 + t212;
t164 = pkin(3) * t120 + t81;
t34 = qJD(1) * t56 - t185 * t79;
t160 = t114 * t226 - t34;
t35 = t228 * qJD(1) + t79 * t186;
t159 = t114 * t58 + t35;
t36 = qJD(1) * t58 - t185 * t80;
t158 = -t114 * t228 - t36;
t37 = t226 * qJD(1) + t80 * t186;
t157 = -t114 * t56 + t37;
t65 = t151 * t114;
t154 = t181 * t65;
t86 = t152 * qJD(3);
t153 = t181 * t86;
t149 = -t107 * t226 - t108 * t228;
t148 = t107 * t80 + t108 * t79;
t144 = -t118 * t225 - t120 * t227;
t143 = Icges(4,1) * t120 - t201;
t140 = -Icges(4,2) * t118 + t200;
t137 = Icges(4,5) * t120 - Icges(4,6) * t118;
t130 = t149 * t119;
t14 = t119 * t150 + t121 * t54;
t15 = -t121 * t230 + t130;
t16 = t119 * t54 - t233;
t17 = -t119 * t230 - t121 * t149;
t32 = -t185 * t78 + t191;
t33 = t230 * qJD(1) + t78 * t186;
t134 = -t180 * (t119 * t15 + t121 * t14) + t119 * ((t119 * t32 + (-t16 + t130) * qJD(1)) * t119 + (t17 * qJD(1) + (-t107 * t37 - t108 * t35 - t187 * t58 + t188 * t56 + t191) * t121 + (t150 * qJD(1) + t158 * t107 + t160 * t108 + t33) * t119) * t121) + t121 * ((t121 * t33 + (t15 + t233) * qJD(1)) * t121 + (-t14 * qJD(1) + (t107 * t36 + t108 * t34 - t187 * t226 + t188 * t228) * t119 + (t32 + t159 * t108 + t157 * t107 + (t149 - t54) * qJD(1)) * t121) * t119) + (t119 * t17 + t121 * t16) * t179;
t133 = rSges(3,3) * t121 + t119 * t218;
t125 = qJD(1) * t148 - t135 * t114;
t132 = (t107 * t160 - t108 * t158 + t125 * t119 + t224 * t121) * t220 + (-t107 * t159 + t108 * t157 - t224 * t119 + t125 * t121) * t219 - (-t107 * t56 + t108 * t58 + t119 * t148 + t121 * t78) * t180 / 0.2e1 + (t107 * t228 - t108 * t226 + t119 * t78 - t121 * t148) * t179 / 0.2e1;
t129 = t144 * t119;
t128 = qJD(3) * t143;
t127 = qJD(3) * t140;
t126 = t151 + t216;
t12 = (qJD(1) * t217 - t169) * t119 - t170 + t172 + t183;
t13 = t106 + (t114 * t81 + t168) * t121 + (t173 * t121 + (-qJ(2) - t126) * t119) * qJD(1);
t110 = t121 * qJ(2);
t39 = t119 * t173 + t121 * t126 + t110;
t40 = -t121 * t122 + t103 + t182 + t60;
t124 = t119 * t13 - t12 * t121 + (t119 * t40 + t121 * t39) * qJD(1);
t28 = t81 * t180 + t121 * t65 + (t118 * t175 + t120 * t180) * pkin(3);
t29 = t81 * t179 - t119 * t65 + (-t118 * t177 + t166) * pkin(3);
t52 = t164 * t119;
t53 = t164 * t121;
t123 = t29 * t119 - t28 * t121 + (-t119 * t53 + t121 * t52) * qJD(1);
t76 = t119 * t215 - t121 * t216;
t75 = -rSges(3,2) * t121 + rSges(3,3) * t119 + t182;
t74 = t110 + t133;
t73 = t206 - t131;
t61 = t119 * rSges(5,3) - t121 * t151;
t51 = t121 * t61;
t50 = t106 + (t218 * t121 + (-rSges(3,3) - qJ(2)) * t119) * qJD(1);
t49 = qJD(1) * t133 + t183;
t48 = t121 * pkin(5) + t182 + t72;
t47 = t119 * t174 + t110 + t131;
t42 = t229 * qJD(1) + t137 * t177;
t41 = -t137 * t175 + t190;
t38 = (-rSges(5,3) * qJD(1) - t169) * t119 + t172;
t31 = -t119 * t60 + t51;
t30 = t121 * (qJD(1) * t60 - t81 * t185);
t21 = -t119 * t229 - t121 * t144;
t20 = t119 * t66 - t234;
t19 = -t121 * t229 + t129;
t18 = t145 * t119 + t121 * t66;
t11 = t119 * t214 + t121 * t76 + t51;
t10 = -t119 * t38 + t30 + (-t119 * t61 - t121 * t60) * qJD(1);
t3 = -t116 * t168 + t30 + (-t38 + t170) * t119 + ((-t119 * pkin(5) - t61 - t76) * t119 + (t214 + t232) * t121) * qJD(1);
t1 = [0.2e1 * m(3) * (t49 * t75 + t50 * t74) - t118 * t128 - t142 * t176 - t120 * t127 + t139 * t178 + (t26 * t48 + t27 * t47) * t223 - t80 * t188 - t108 * t64 - t79 * t187 + t107 * t63 + (t12 * t40 + t13 * t39) * t222; m(3) * (t119 * t50 - t121 * t49 + (t119 * t75 + t121 * t74) * qJD(1)) + m(4) * ((t119 * t48 + t121 * t47) * qJD(1) + t231) + m(5) * t124; 0; (-qJD(3) * t145 - t118 * (t227 * qJD(1) + t119 * t127) + t120 * (t225 * qJD(1) + t119 * t128)) * t219 + (-qJD(3) * t144 - t118 * (qJD(1) * t68 - t140 * t175) + t120 * (qJD(1) * t70 - t143 * t175)) * t220 + m(4) * (t231 * t92 - (t119 * t47 - t121 * t48) * t86) + m(5) * (-t12 * t53 + t13 * t52 + t28 * t40 + t29 * t39) - (t116 / 0.2e1 + t115 / 0.2e1) * t136 * qJD(3) + ((t208 / 0.2e1 - t204 / 0.2e1 + t48 * t221) * t119 + (t207 / 0.2e1 - t203 / 0.2e1 + t47 * t221) * t121) * qJD(1) + t132; -m(4) * t153 + m(5) * t123; ((-t119 * t72 + t121 * t73) * (-t119 * t171 + (t115 * t211 - t116 * t92) * qJD(3) + ((-t72 + t112) * t121 + (t131 - t73 + t206) * t119) * qJD(1)) - t92 * t153) * t223 - (t119 * t19 + t18 * t121) * t180 + t121 * ((t121 * t42 + (t19 + t234) * qJD(1)) * t121 + (-t18 * qJD(1) + (-t176 * t225 + t178 * t227) * t119 + (t41 + (t204 - t208) * qJD(3) + (t144 - t66) * qJD(1)) * t121) * t119) + (t119 * t21 + t121 * t20) * t179 + t119 * ((t119 * t41 + (-t20 + t129) * qJD(1)) * t119 + (t21 * qJD(1) + (-t176 * t70 + t178 * t68 + t190) * t121 + (t42 + (t203 - t207) * qJD(3) + t145 * qJD(1)) * t119) * t121) + (t11 * t3 - t28 * t53 + t29 * t52) * t222 + t134; m(5) * (-(t119 * t39 - t121 * t40) * t65 + t124 * t81) + t132; -m(5) * t154; m(5) * (t10 * t11 + t31 * t3 - (t119 * t52 + t121 * t53) * t65 + t123 * t81) + t134; (t31 * t10 - t154 * t81) * t222 + t134;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
