% Calculate time derivative of joint inertia matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:13:01
% DurationCPUTime: 4.65s
% Computational Cost: add. (3253->268), mult. (3606->378), div. (0->0), fcn. (2662->6), ass. (0->157)
t143 = cos(qJ(3));
t141 = sin(qJ(3));
t216 = Icges(5,4) * t141;
t218 = Icges(4,4) * t141;
t256 = t216 + t218 + (Icges(4,2) + Icges(5,2)) * t143;
t215 = Icges(5,4) * t143;
t217 = Icges(4,4) * t143;
t255 = t215 + t217 + (Icges(4,1) + Icges(5,1)) * t141;
t166 = -Icges(5,2) * t141 + t215;
t167 = -Icges(4,2) * t141 + t217;
t168 = Icges(5,1) * t143 - t216;
t169 = Icges(4,1) * t143 - t218;
t254 = (t166 + t167) * t143 + (t168 + t169) * t141;
t139 = qJ(1) + qJ(2);
t136 = cos(t139);
t138 = qJD(1) + qJD(2);
t211 = t136 * t138;
t251 = t256 * t141 - t255 * t143;
t250 = rSges(5,1) + pkin(3);
t134 = pkin(3) * t143 + pkin(2);
t135 = sin(t139);
t140 = -qJ(4) - pkin(6);
t209 = t136 * t141;
t126 = t135 * rSges(5,3);
t208 = t136 * t143;
t237 = rSges(5,1) * t208 + t126;
t50 = -rSges(5,2) * t209 + t136 * t134 - t135 * t140 + t237;
t108 = Icges(5,5) * t141 + Icges(5,6) * t143;
t109 = Icges(4,5) * t141 + Icges(4,6) * t143;
t249 = t108 + t109;
t201 = qJD(3) * t143;
t206 = t138 * t141;
t246 = -t135 * t201 - t136 * t206;
t213 = t135 * t141;
t245 = rSges(5,2) * t213 + (rSges(5,3) - t140) * t136;
t164 = Icges(5,5) * t143 - Icges(5,6) * t141;
t165 = Icges(4,5) * t143 - Icges(4,6) * t141;
t244 = t251 * t138 + (t164 + t165) * qJD(3);
t132 = t136 * pkin(6);
t212 = t135 * t143;
t243 = t132 + t135 * (-pkin(2) + t134) + rSges(5,1) * t212 - t245;
t203 = -t136 * pkin(2) - t135 * pkin(6);
t242 = t203 + t50;
t159 = t167 * t136;
t66 = Icges(4,6) * t135 + t159;
t161 = t169 * t136;
t70 = Icges(4,5) * t135 + t161;
t172 = t141 * t66 - t143 * t70;
t241 = t135 * t172;
t158 = t166 * t136;
t64 = Icges(5,6) * t135 + t158;
t160 = t168 * t136;
t68 = Icges(5,5) * t135 + t160;
t176 = t141 * t64 - t143 * t68;
t240 = t135 * t176;
t65 = -Icges(4,6) * t136 + t135 * t167;
t69 = -Icges(4,5) * t136 + t135 * t169;
t174 = t141 * t65 - t143 * t69;
t239 = t136 * t174;
t63 = -Icges(5,6) * t136 + t135 * t166;
t67 = -Icges(5,5) * t136 + t135 * t168;
t178 = t141 * t63 - t143 * t67;
t238 = t136 * t178;
t196 = t135 * t206;
t236 = rSges(5,2) * t196 + rSges(5,3) * t211 + qJD(4) * t135;
t202 = qJD(3) * t141;
t193 = t135 * t202;
t207 = t138 * t140;
t235 = -rSges(5,2) * t246 + qJD(4) * t136 + t135 * t207 + t250 * t193;
t148 = t254 * qJD(3) + t255 * t201 - t256 * t202;
t232 = 2 * m(3);
t231 = 2 * m(4);
t230 = 2 * m(5);
t224 = rSges(4,1) * t143;
t100 = (-rSges(4,2) * t141 + t224) * qJD(3);
t227 = m(4) * t100;
t119 = rSges(4,1) * t141 + rSges(4,2) * t143;
t226 = m(4) * t119;
t142 = sin(qJ(1));
t225 = pkin(1) * t142;
t223 = rSges(5,1) * t143;
t222 = rSges(5,2) * t143;
t221 = pkin(1) * qJD(1);
t127 = t135 * rSges(4,3);
t219 = rSges(4,2) * t196 + rSges(4,3) * t211;
t214 = t135 * t138;
t204 = rSges(4,2) * t213 + t136 * rSges(4,3);
t200 = t142 * t221;
t144 = cos(qJ(1));
t199 = t144 * t221;
t197 = -rSges(4,1) * t193 + rSges(4,2) * t246;
t191 = t136 * t202;
t190 = t136 * t201;
t187 = -pkin(2) - t224;
t118 = rSges(5,1) * t141 + t222;
t186 = -pkin(3) * t141 - t118;
t185 = -t134 - t223;
t88 = t136 * rSges(3,1) - rSges(3,2) * t135;
t78 = -rSges(3,1) * t211 + rSges(3,2) * t214;
t87 = -rSges(3,1) * t135 - rSges(3,2) * t136;
t179 = t141 * t67 + t143 * t63;
t177 = t141 * t68 + t143 * t64;
t175 = t141 * t69 + t143 * t65;
t173 = t141 * t70 + t143 * t66;
t171 = t185 * t135;
t74 = rSges(4,1) * t208 - rSges(4,2) * t209 + t127;
t77 = t87 * t138;
t157 = t165 * t136;
t156 = t164 * t136;
t155 = -t138 * t212 - t191;
t54 = t74 - t203;
t53 = t135 * t187 + t132 + t204;
t150 = Icges(4,3) * t138 - qJD(3) * t109;
t149 = Icges(5,3) * t138 - qJD(3) * t108;
t49 = t171 + t245;
t147 = (t244 * t135 + (-t172 - t176) * qJD(3) - t254 * t214) * t135 / 0.2e1 - (-t244 * t136 + (-t174 - t178) * qJD(3)) * t136 / 0.2e1 + (-t135 * t251 - t136 * t249 + t175 + t179) * t214 / 0.2e1 + (t135 * t249 - t136 * t251 + t173 + t177) * t211 / 0.2e1 - ((t158 + t159) * t143 + (t160 + t161) * t141) * t211 / 0.2e1;
t24 = (t187 * t136 + (-rSges(4,3) - pkin(6)) * t135) * t138 - t197;
t12 = (t136 * t185 - t126) * t138 + t235;
t117 = pkin(6) * t211;
t23 = rSges(4,1) * t155 - rSges(4,2) * t190 - pkin(2) * t214 + t117 + t219;
t11 = t138 * t171 + (-t207 + (-t250 * t141 - t222) * qJD(3)) * t136 + t236;
t137 = t144 * pkin(1);
t99 = (-rSges(5,2) * t141 + t223) * qJD(3);
t80 = t137 + t88;
t79 = t87 - t225;
t76 = t186 * t136;
t75 = t186 * t135;
t72 = rSges(4,1) * t212 - t204;
t62 = Icges(4,3) * t135 + t157;
t61 = -Icges(4,3) * t136 + t135 * t165;
t60 = Icges(5,3) * t135 + t156;
t59 = -Icges(5,3) * t136 + t135 * t164;
t58 = t78 - t199;
t57 = t77 - t200;
t52 = t137 + t54;
t51 = t53 - t225;
t48 = t137 + t50;
t47 = t49 - t225;
t38 = t135 * t150 + t138 * t157;
t37 = t136 * t150 - t165 * t214;
t36 = t135 * t149 + t138 * t156;
t35 = t136 * t149 - t164 * t214;
t34 = pkin(3) * t246 - t118 * t211 - t135 * t99;
t33 = t118 * t214 - t136 * t99 + (-t190 + t196) * pkin(3);
t22 = t24 - t199;
t21 = t23 - t200;
t20 = t135 * t62 - t136 * t172;
t19 = t135 * t61 - t239;
t18 = t135 * t60 - t136 * t176;
t17 = t135 * t59 - t238;
t16 = -t136 * t62 - t241;
t15 = -t135 * t174 - t136 * t61;
t14 = -t136 * t60 - t240;
t13 = -t135 * t178 - t136 * t59;
t10 = t12 - t199;
t9 = t11 - t200;
t1 = [(t57 * t80 + t58 * t79) * t232 + (t21 * t52 + t22 * t51) * t231 + (t10 * t47 + t48 * t9) * t230 + t148; m(3) * (t57 * t88 + t58 * t87 + t77 * t80 + t78 * t79) + m(4) * (t21 * t54 + t53 * t22 + t23 * t52 + t24 * t51) + m(5) * (t10 * t49 + t11 * t48 + t12 * t47 + t50 * t9) + t148; (t77 * t88 + t78 * t87) * t232 + (t23 * t54 + t53 * t24) * t231 + (t11 * t50 + t12 * t49) * t230 + t148; t147 + ((-t138 * t52 - t22) * t136 + (t138 * t51 - t21) * t135) * t226 + (-t135 * t52 - t136 * t51) * t227 + m(5) * (t10 * t76 + t33 * t47 + t34 * t48 + t75 * t9); t147 + m(5) * (t11 * t75 + t12 * t76 + t33 * t49 + t34 * t50) + ((-t138 * t54 - t24) * t136 + (t138 * t53 - t23) * t135) * t226 + (-t135 * t54 - t136 * t53) * t227; ((t135 * t72 + t136 * t74) * (((-t74 + t127) * t138 + t197) * t135 + (-qJD(3) * t119 * t136 + t138 * t72 + t219) * t136) + (t135 ^ 2 + t136 ^ 2) * t119 * t100) * t231 + (t135 * t20 - t136 * t19) * t211 + t135 * ((t135 * t37 + (t19 + t241) * t138) * t135 + (t20 * t138 + (t201 * t65 + t202 * t69) * t136 + (-t173 * qJD(3) - t174 * t138 - t38) * t135) * t136) + (t135 * t16 - t136 * t15) * t214 - t136 * ((t136 * t38 + (t16 + t239) * t138) * t136 + (t15 * t138 + (-t201 * t66 - t202 * t70) * t135 + (t175 * qJD(3) - t138 * t172 - t37) * t136) * t135) + (t76 * t33 + t75 * t34 + (t135 * t243 + t242 * t136) * (-t242 * t214 + t243 * t211 + (rSges(5,1) * t155 - rSges(5,2) * t190 - pkin(3) * t191 - t136 * t207 - t117 + t236) * t136 + (-pkin(6) * t214 + t138 * t237 - t235) * t135)) * t230 + (t135 * t18 - t136 * t17) * t211 + t135 * ((t135 * t35 + (t17 + t240) * t138) * t135 + (t18 * t138 + (t201 * t63 + t202 * t67) * t136 + (-t177 * qJD(3) - t178 * t138 - t36) * t135) * t136) + (-t13 * t136 + t135 * t14) * t214 - t136 * ((t136 * t36 + (t14 + t238) * t138) * t136 + (t13 * t138 + (-t201 * t64 - t202 * t68) * t135 + (t179 * qJD(3) - t138 * t176 - t35) * t136) * t135); m(5) * ((t138 * t47 - t9) * t136 + (t138 * t48 + t10) * t135); m(5) * ((t138 * t49 - t11) * t136 + (t138 * t50 + t12) * t135); m(5) * ((t138 * t76 - t34) * t136 + (t138 * t75 + t33) * t135); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
