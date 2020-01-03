% Calculate time derivative of joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:16
% DurationCPUTime: 4.32s
% Computational Cost: add. (4968->273), mult. (3832->376), div. (0->0), fcn. (2806->8), ass. (0->152)
t153 = cos(qJ(4));
t151 = sin(qJ(4));
t234 = Icges(6,4) * t151;
t236 = Icges(5,4) * t151;
t273 = t234 + t236 + (Icges(5,2) + Icges(6,2)) * t153;
t233 = Icges(6,4) * t153;
t235 = Icges(5,4) * t153;
t272 = t233 + t235 + (Icges(5,1) + Icges(6,1)) * t151;
t180 = -Icges(6,2) * t151 + t233;
t181 = -Icges(5,2) * t151 + t235;
t182 = Icges(6,1) * t153 - t234;
t183 = Icges(5,1) * t153 - t236;
t271 = (t180 + t181) * t153 + (t182 + t183) * t151;
t149 = qJ(1) + pkin(8);
t212 = pkin(2) * cos(t149) + cos(qJ(1)) * pkin(1);
t267 = t273 * t151 - t272 * t153;
t148 = qJD(1) + qJD(3);
t178 = Icges(6,5) * t153 - Icges(6,6) * t151;
t179 = Icges(5,5) * t153 - Icges(5,6) * t151;
t270 = (-t178 - t179) * qJD(4) + (-t267 + t271) * t148;
t118 = Icges(6,5) * t151 + Icges(6,6) * t153;
t119 = Icges(5,5) * t151 + Icges(5,6) * t153;
t269 = -t118 - t119;
t145 = qJ(3) + t149;
t139 = sin(t145);
t140 = cos(t145);
t141 = pkin(4) * t153 + pkin(3);
t220 = t140 * t153;
t221 = t140 * t151;
t150 = -qJ(5) - pkin(7);
t225 = t139 * t150;
t54 = rSges(6,1) * t220 - rSges(6,2) * t221 + t139 * rSges(6,3) + t140 * t141 - t225;
t124 = t140 * t150;
t223 = t139 * t153;
t224 = t139 * t151;
t53 = rSges(6,1) * t223 - rSges(6,2) * t224 - rSges(6,3) * t140 + t139 * t141 + t124;
t68 = -Icges(5,6) * t139 - t140 * t181;
t72 = -Icges(5,5) * t139 - t140 * t183;
t185 = t151 * t68 - t153 * t72;
t265 = t139 * t185;
t66 = -Icges(6,6) * t139 - t140 * t180;
t70 = -Icges(6,5) * t139 - t140 * t182;
t189 = t151 * t66 - t153 * t70;
t264 = t139 * t189;
t263 = -t139 * rSges(4,1) - t140 * rSges(4,2);
t262 = -pkin(2) * sin(t149) - sin(qJ(1)) * pkin(1);
t260 = -Icges(6,3) * t148 + qJD(4) * t118;
t259 = -Icges(5,3) * t148 + qJD(4) * t119;
t210 = qJD(4) * t153;
t211 = qJD(4) * t151;
t157 = t271 * qJD(4) + t272 * t210 - t273 * t211;
t252 = 2 * m(4);
t251 = 2 * m(5);
t250 = 2 * m(6);
t241 = rSges(5,2) * t151;
t243 = rSges(5,1) * t153;
t105 = (-t241 + t243) * qJD(4);
t247 = m(5) * t105;
t126 = rSges(5,1) * t151 + rSges(5,2) * t153;
t246 = m(5) * t126;
t135 = t139 * pkin(3);
t245 = pkin(7) * t140 - t135 + t53;
t214 = t140 * pkin(3) + t139 * pkin(7);
t244 = -t214 + t54;
t242 = rSges(6,1) * t153;
t240 = rSges(6,2) * t151;
t239 = rSges(6,2) * t153;
t218 = t148 * t153;
t206 = t140 * t218;
t226 = t139 * t148;
t238 = rSges(5,1) * t206 + rSges(5,3) * t226;
t219 = t148 * t151;
t207 = t139 * t219;
t222 = t140 * t148;
t237 = rSges(5,2) * t207 + rSges(5,3) * t222;
t216 = pkin(3) * t222 + pkin(7) * t226;
t215 = t212 * qJD(1);
t209 = rSges(6,1) * t206 + rSges(6,3) * t226 + t141 * t222;
t208 = rSges(6,2) * t207 + rSges(6,3) * t222 + qJD(5) * t139;
t205 = t140 * t210;
t125 = rSges(6,1) * t151 + t239;
t202 = pkin(4) * t151 + t125;
t201 = -t141 - t242;
t84 = t140 * rSges(4,1) - rSges(4,2) * t139;
t82 = rSges(4,1) * t222 - rSges(4,2) * t226;
t199 = rSges(5,1) * t223 - rSges(5,2) * t224;
t65 = -Icges(6,6) * t140 + t139 * t180;
t69 = -Icges(6,5) * t140 + t139 * t182;
t192 = t151 * t69 + t153 * t65;
t191 = t151 * t65 - t153 * t69;
t190 = t151 * t70 + t153 * t66;
t67 = -Icges(5,6) * t140 + t139 * t181;
t71 = -Icges(5,5) * t140 + t139 * t183;
t188 = t151 * t71 + t153 * t67;
t187 = t151 * t67 - t153 * t71;
t186 = t151 * t72 + t153 * t68;
t76 = -rSges(5,1) * t220 + rSges(5,2) * t221 - t139 * rSges(5,3);
t173 = t262 * qJD(1);
t81 = t263 * t148;
t172 = t191 * t140;
t171 = t187 * t140;
t170 = qJD(4) * t126;
t165 = t179 * t148;
t164 = t178 * t148;
t163 = -t139 * t210 - t140 * t219;
t56 = -t76 + t214;
t160 = (-t239 + (-rSges(6,1) - pkin(4)) * t151) * qJD(4);
t55 = t135 + (-rSges(5,3) - pkin(7)) * t140 + t199;
t156 = -t148 * t150 + t160;
t155 = -((-t185 - t189) * qJD(4) + t270 * t139) * t139 / 0.2e1 - ((-t187 - t191) * qJD(4) + t270 * t140) * t140 / 0.2e1 + (-t267 * t139 + t269 * t140 + t188 + t192) * t226 / 0.2e1 - (t269 * t139 + t267 * t140 + t186 + t190) * t222 / 0.2e1;
t26 = -rSges(5,1) * t139 * t211 + rSges(5,2) * t163 + t216 + t238;
t112 = pkin(7) * t222;
t25 = -rSges(5,2) * t205 - pkin(3) * t226 + t112 + (-t139 * t218 - t140 * t211) * rSges(5,1) + t237;
t22 = (-rSges(6,2) * t219 - qJD(5)) * t140 + t156 * t139 + t209;
t21 = t140 * t156 + t201 * t226 + t208;
t104 = (-t240 + t242) * qJD(4);
t80 = t84 + t212;
t79 = -t262 - t263;
t78 = t202 * t140;
t77 = t202 * t139;
t74 = -rSges(5,3) * t140 + t199;
t64 = -Icges(5,3) * t139 - t140 * t179;
t63 = -Icges(5,3) * t140 + t139 * t179;
t62 = -Icges(6,3) * t139 - t140 * t178;
t61 = -Icges(6,3) * t140 + t139 * t178;
t60 = t82 + t215;
t59 = t81 + t173;
t52 = t56 + t212;
t51 = t55 - t262;
t50 = t54 + t212;
t49 = t53 - t262;
t40 = -t259 * t139 + t140 * t165;
t39 = t139 * t165 + t259 * t140;
t38 = -t260 * t139 + t140 * t164;
t37 = t139 * t164 + t260 * t140;
t36 = pkin(4) * t163 - t104 * t139 - t125 * t222;
t35 = -t125 * t226 + t104 * t140 + (t205 - t207) * pkin(4);
t24 = t26 + t215;
t23 = t173 + t25;
t20 = -t139 * t64 + t185 * t140;
t19 = -t139 * t63 + t171;
t18 = -t139 * t62 + t189 * t140;
t17 = -t139 * t61 + t172;
t16 = -t140 * t64 - t265;
t15 = -t187 * t139 - t140 * t63;
t14 = -t140 * t62 - t264;
t13 = -t191 * t139 - t140 * t61;
t12 = t22 + t215;
t11 = t173 + t21;
t2 = (-t140 * t170 + t148 * t74 + t237) * t140 + (-t139 * t170 + (t76 + (-t241 - t243) * t140) * t148 + t238) * t139;
t1 = (-t112 + (-t124 + t245) * t148 + t140 * t160 + t208) * t140 + (-qJD(5) * t140 + t139 * t160 + (-t225 + (pkin(3) + t201 - t240) * t140 - t244) * t148 + t209 - t216) * t139;
t3 = [(t59 * t80 + t60 * t79) * t252 + (t23 * t52 + t24 * t51) * t251 + (t11 * t50 + t12 * t49) * t250 + t157; 0; 0; m(4) * (-t263 * t60 + t59 * t84 + t79 * t82 + t80 * t81) + m(5) * (t23 * t56 + t24 * t55 + t25 * t52 + t26 * t51) + m(6) * (t11 * t54 + t12 * t53 + t21 * t50 + t22 * t49) + t157; 0; (t21 * t54 + t22 * t53) * t250 + (t25 * t56 + t26 * t55) * t251 + (-t263 * t82 + t81 * t84) * t252 + t157; m(6) * (-t11 * t77 + t12 * t78 + t35 * t49 + t36 * t50) + ((-t148 * t52 + t24) * t140 + (-t148 * t51 - t23) * t139) * t246 + (-t139 * t52 + t140 * t51) * t247 + t155; m(5) * t2 + m(6) * t1; ((-t148 * t56 + t26) * t140 + (-t148 * t55 - t25) * t139) * t246 + (-t139 * t56 + t140 * t55) * t247 + m(6) * (-t21 * t77 + t22 * t78 + t35 * t53 + t36 * t54) + t155; ((t139 * t74 - t140 * t76) * t2 + (t139 ^ 2 + t140 ^ 2) * t126 * t105) * t251 - t140 * ((t140 * t40 + (-t16 + t171) * t148) * t140 + (t15 * t148 + (t210 * t68 + t211 * t72) * t139 + (t188 * qJD(4) + t185 * t148 + t39) * t140) * t139) - t139 * ((t139 * t39 + (t19 + t265) * t148) * t139 + (-t20 * t148 + (-t210 * t67 - t211 * t71) * t140 + (-t186 * qJD(4) + t148 * t187 + t40) * t139) * t140) + ((t139 * t245 + t140 * t244) * t1 - t77 * t36 + t78 * t35) * t250 - t140 * ((t140 * t38 + (-t14 + t172) * t148) * t140 + (t13 * t148 + (t210 * t66 + t211 * t70) * t139 + (t192 * qJD(4) + t148 * t189 + t37) * t140) * t139) - t139 * ((t139 * t37 + (t17 + t264) * t148) * t139 + (-t18 * t148 + (-t210 * t65 - t211 * t69) * t140 + (-t190 * qJD(4) + t148 * t191 + t38) * t139) * t140) + ((-t13 - t15) * t140 + (-t14 - t16) * t139) * t226 + ((t17 + t19) * t140 + (t18 + t20) * t139) * t222; m(6) * ((-t148 * t49 - t11) * t140 + (t148 * t50 - t12) * t139); 0; m(6) * ((-t148 * t53 - t21) * t140 + (t148 * t54 - t22) * t139); m(6) * ((-t148 * t78 - t36) * t140 + (-t148 * t77 - t35) * t139); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
