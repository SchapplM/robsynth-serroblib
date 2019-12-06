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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:31
% DurationCPUTime: 3.52s
% Computational Cost: add. (4968->272), mult. (3832->370), div. (0->0), fcn. (2806->8), ass. (0->159)
t149 = sin(qJ(4));
t151 = cos(qJ(4));
t232 = Icges(6,4) * t151;
t179 = -Icges(6,2) * t149 + t232;
t234 = Icges(5,4) * t151;
t180 = -Icges(5,2) * t149 + t234;
t233 = Icges(6,4) * t149;
t181 = Icges(6,1) * t151 - t233;
t235 = Icges(5,4) * t149;
t182 = Icges(5,1) * t151 - t235;
t281 = (t179 + t180) * t151 + (t181 + t182) * t149;
t280 = t233 + t235 + (Icges(5,2) + Icges(6,2)) * t151;
t279 = t232 + t234 + (Icges(5,1) + Icges(6,1)) * t149;
t272 = t280 * t149 - t279 * t151;
t146 = qJD(1) + qJD(3);
t177 = Icges(6,5) * t151 - Icges(6,6) * t149;
t178 = Icges(5,5) * t151 - Icges(5,6) * t149;
t276 = (t177 + t178) * qJD(4) + (t272 - t281) * t146;
t147 = qJ(1) + pkin(8);
t275 = pkin(2) * sin(t147) + sin(qJ(1)) * pkin(1);
t124 = Icges(6,5) * t149 + Icges(6,6) * t151;
t125 = Icges(5,5) * t149 + Icges(5,6) * t151;
t274 = t124 + t125;
t145 = qJ(3) + t147;
t140 = cos(t145);
t213 = qJD(4) * t149;
t205 = t140 * t213;
t139 = sin(t145);
t223 = t139 * t151;
t271 = -t146 * t223 - t205;
t212 = qJD(4) * t151;
t217 = t146 * t149;
t270 = t139 * t212 + t140 * t217;
t154 = t281 * qJD(4) + t279 * t212 - t280 * t213;
t67 = Icges(5,6) * t140 - t139 * t180;
t71 = Icges(5,5) * t140 - t139 * t182;
t186 = t149 * t67 - t151 * t71;
t268 = t140 * t186;
t65 = Icges(6,6) * t140 - t139 * t179;
t69 = Icges(6,5) * t140 - t139 * t181;
t190 = t149 * t65 - t151 * t69;
t267 = t140 * t190;
t148 = -qJ(5) - pkin(7);
t220 = t140 * t149;
t266 = rSges(6,2) * t220 + t139 * t148;
t224 = t139 * t149;
t265 = rSges(6,2) * t224 + (rSges(6,3) - t148) * t140;
t264 = -Icges(6,3) * t146 + qJD(4) * t124;
t263 = -Icges(5,3) * t146 + qJD(4) * t125;
t256 = 2 * m(4);
t255 = 2 * m(5);
t254 = 2 * m(6);
t251 = -rSges(5,3) - pkin(7);
t250 = -pkin(7) + rSges(6,3);
t239 = rSges(5,2) * t149;
t241 = rSges(5,1) * t151;
t105 = (-t239 + t241) * qJD(4);
t249 = m(5) * t105;
t132 = rSges(5,1) * t149 + rSges(5,2) * t151;
t248 = m(5) * t132;
t141 = pkin(4) * t151 + pkin(3);
t244 = pkin(3) - t141;
t138 = t140 * pkin(7);
t243 = rSges(6,1) * t223 - t139 * t244 + t138 - t265;
t219 = t140 * t151;
t236 = t139 * rSges(6,3);
t242 = rSges(6,1) * t219 - t139 * pkin(7) - t140 * t244 + t236 - t266;
t240 = rSges(6,1) * t151;
t238 = rSges(6,2) * t149;
t237 = t139 * rSges(5,3);
t137 = t140 * rSges(5,3);
t225 = t139 * t146;
t222 = t140 * t146;
t218 = t146 * t148;
t81 = rSges(4,1) * t225 + rSges(4,2) * t222;
t121 = rSges(5,2) * t224;
t215 = t121 + t137;
t214 = t275 * qJD(1);
t133 = qJD(5) * t140;
t207 = t139 * t213;
t211 = rSges(5,1) * t207 + t270 * rSges(5,2);
t204 = t140 * t212;
t210 = t271 * rSges(5,1) - rSges(5,2) * t204;
t201 = -pkin(3) - t241;
t131 = rSges(6,1) * t149 + rSges(6,2) * t151;
t200 = pkin(4) * t149 + t131;
t84 = -t140 * rSges(4,1) + t139 * rSges(4,2);
t199 = -t141 - t240;
t82 = -rSges(4,1) * t222 + rSges(4,2) * t225;
t197 = -pkin(2) * cos(t147) - cos(qJ(1)) * pkin(1);
t83 = -rSges(4,1) * t139 - rSges(4,2) * t140;
t191 = t149 * t69 + t151 * t65;
t66 = Icges(6,6) * t139 + t140 * t179;
t70 = Icges(6,5) * t139 + t140 * t181;
t189 = t149 * t70 + t151 * t66;
t188 = t149 * t66 - t151 * t70;
t187 = t149 * t71 + t151 * t67;
t68 = Icges(5,6) * t139 + t140 * t180;
t72 = Icges(5,5) * t139 + t140 * t182;
t185 = t149 * t72 + t151 * t68;
t184 = t149 * t68 - t151 * t72;
t172 = t271 * rSges(6,1) - rSges(6,2) * t204 - pkin(4) * t205 - t140 * t218 - t141 * t225;
t171 = t139 * t218 + t133 + (rSges(6,1) + pkin(4)) * t207 + t270 * rSges(6,2);
t170 = t197 * qJD(1);
t169 = t188 * t139;
t168 = t184 * t139;
t163 = t178 * t146;
t162 = t177 * t146;
t158 = t140 * t199 - t236;
t55 = t139 * t201 + t138 + t215;
t155 = t139 * t251 + t140 * t201;
t123 = rSges(5,2) * t220;
t56 = t123 + t155;
t54 = t158 + t266;
t53 = t139 * t199 + t265;
t153 = ((-t184 - t188) * qJD(4) + t276 * t139) * t139 / 0.2e1 + ((-t186 - t190) * qJD(4) + t276 * t140) * t140 / 0.2e1 - (t272 * t139 + t274 * t140 + t187 + t191) * t225 / 0.2e1 + (t274 * t139 - t272 * t140 + t185 + t189) * t222 / 0.2e1;
t119 = pkin(3) * t225;
t25 = t119 + (t140 * t251 - t121) * t146 - t210;
t26 = t146 * t155 + t211;
t21 = -rSges(6,3) * t222 + (-rSges(6,2) * t217 - qJD(5)) * t139 - t172;
t22 = t146 * t158 + t171;
t104 = (-t238 + t240) * qJD(4);
t80 = t197 + t84;
t79 = -t275 + t83;
t78 = t200 * t140;
t77 = t200 * t139;
t76 = rSges(5,1) * t219 - t123 + t237;
t74 = -rSges(5,1) * t223 + t215;
t64 = Icges(5,3) * t139 + t140 * t178;
t63 = Icges(5,3) * t140 - t139 * t178;
t62 = Icges(6,3) * t139 + t140 * t177;
t61 = Icges(6,3) * t140 - t139 * t177;
t60 = t170 + t82;
t59 = t214 + t81;
t52 = t197 + t56;
t51 = -t275 + t55;
t50 = t197 + t54;
t49 = -t275 + t53;
t40 = t263 * t139 - t140 * t163;
t39 = -t139 * t163 - t263 * t140;
t38 = t264 * t139 - t140 * t162;
t37 = -t139 * t162 - t264 * t140;
t36 = t270 * pkin(4) + t104 * t139 + t131 * t222;
t35 = t131 * t225 - t104 * t140 + (t139 * t217 - t204) * pkin(4);
t24 = t170 + t26;
t23 = t25 + t214;
t20 = t139 * t64 - t184 * t140;
t19 = t139 * t63 - t268;
t18 = t139 * t62 - t188 * t140;
t17 = t139 * t61 - t267;
t16 = t140 * t64 + t168;
t15 = t186 * t139 + t140 * t63;
t14 = t140 * t62 + t169;
t13 = t190 * t139 + t140 * t61;
t12 = t170 + t22;
t11 = t21 + t214;
t2 = t140 * t210 - t139 * t211 + ((-t74 + t137) * t140 + (t237 - t76 + (t239 + t241) * t140) * t139) * t146;
t1 = (t119 + t172) * t140 + (-t171 + t133) * t139 + ((t250 * t139 - t242) * t139 + (t250 * t140 + (-pkin(3) - t199 + t238) * t139 + t243) * t140) * t146;
t3 = [(t59 * t80 + t60 * t79) * t256 + (t23 * t52 + t24 * t51) * t255 + (t11 * t50 + t12 * t49) * t254 + t154; 0; 0; m(4) * (t59 * t84 + t60 * t83 + t79 * t82 + t80 * t81) + m(5) * (t23 * t56 + t24 * t55 + t25 * t52 + t26 * t51) + m(6) * (t11 * t54 + t12 * t53 + t21 * t50 + t22 * t49) + t154; 0; (t21 * t54 + t22 * t53) * t254 + (t25 * t56 + t26 * t55) * t255 + (t81 * t84 + t82 * t83) * t256 + t154; m(6) * (t11 * t77 - t12 * t78 + t35 * t49 + t36 * t50) + ((t146 * t52 - t24) * t140 + (t146 * t51 + t23) * t139) * t248 + (t139 * t52 - t140 * t51) * t249 + t153; m(5) * t2 + m(6) * t1; m(6) * (t21 * t77 - t22 * t78 + t35 * t53 + t36 * t54) + ((t146 * t56 - t26) * t140 + (t146 * t55 + t25) * t139) * t248 + (t139 * t56 - t140 * t55) * t249 + t153; ((-t139 * t74 + t140 * t76) * t2 + (t139 ^ 2 + t140 ^ 2) * t132 * t105) * t255 + t140 * ((t140 * t40 + (t16 + t268) * t146) * t140 + (-t15 * t146 + (t212 * t68 + t213 * t72) * t139 + (t187 * qJD(4) + t146 * t184 + t39) * t140) * t139) + t139 * ((t139 * t39 + (-t19 + t168) * t146) * t139 + (t20 * t146 + (-t212 * t67 - t213 * t71) * t140 + (-t185 * qJD(4) + t146 * t186 + t40) * t139) * t140) + ((t139 * t243 + t140 * t242) * t1 + t77 * t36 - t78 * t35) * t254 + t140 * ((t140 * t38 + (t14 + t267) * t146) * t140 + (-t13 * t146 + (t212 * t66 + t213 * t70) * t139 + (t191 * qJD(4) + t146 * t188 + t37) * t140) * t139) + t139 * ((t139 * t37 + (-t17 + t169) * t146) * t139 + (t18 * t146 + (-t212 * t65 - t213 * t69) * t140 + (-t189 * qJD(4) + t146 * t190 + t38) * t139) * t140) + ((-t13 - t15) * t140 + (-t14 - t16) * t139) * t225 + ((t17 + t19) * t140 + (t18 + t20) * t139) * t222; m(6) * ((t146 * t49 + t11) * t140 + (-t146 * t50 + t12) * t139); 0; m(6) * ((t146 * t53 + t21) * t140 + (-t146 * t54 + t22) * t139); m(6) * ((-t146 * t78 + t36) * t140 + (-t146 * t77 + t35) * t139); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
