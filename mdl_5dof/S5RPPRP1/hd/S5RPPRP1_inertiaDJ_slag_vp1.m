% Calculate time derivative of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:09
% EndTime: 2020-01-03 11:25:22
% DurationCPUTime: 4.68s
% Computational Cost: add. (5212->389), mult. (7260->548), div. (0->0), fcn. (6751->8), ass. (0->176)
t144 = qJ(1) + pkin(7);
t141 = cos(t144);
t150 = cos(qJ(4));
t140 = sin(t144);
t146 = cos(pkin(8));
t148 = sin(qJ(4));
t192 = t146 * t148;
t175 = t140 * t192;
t108 = -t141 * t150 - t175;
t191 = t146 * t150;
t195 = t141 * t148;
t109 = t140 * t191 - t195;
t145 = sin(pkin(8));
t198 = t140 * t145;
t65 = Icges(6,5) * t109 + Icges(6,6) * t108 + Icges(6,3) * t198;
t67 = Icges(5,5) * t109 + Icges(5,6) * t108 + Icges(5,3) * t198;
t228 = t65 + t67;
t110 = -t140 * t150 + t141 * t192;
t197 = t140 * t148;
t155 = t141 * t191 + t197;
t196 = t141 * t145;
t66 = -Icges(6,5) * t155 + Icges(6,6) * t110 - Icges(6,3) * t196;
t68 = -Icges(5,5) * t155 + Icges(5,6) * t110 - Icges(5,3) * t196;
t227 = t66 + t68;
t69 = Icges(6,4) * t109 + Icges(6,2) * t108 + Icges(6,6) * t198;
t71 = Icges(5,4) * t109 + Icges(5,2) * t108 + Icges(5,6) * t198;
t226 = t71 + t69;
t70 = -Icges(6,4) * t155 + Icges(6,2) * t110 - Icges(6,6) * t196;
t72 = -Icges(5,4) * t155 + Icges(5,2) * t110 - Icges(5,6) * t196;
t225 = t72 + t70;
t74 = -Icges(6,1) * t155 + Icges(6,4) * t110 - Icges(6,5) * t196;
t76 = -Icges(5,1) * t155 + Icges(5,4) * t110 - Icges(5,5) * t196;
t224 = t74 + t76;
t73 = Icges(6,1) * t109 + Icges(6,4) * t108 + Icges(6,5) * t198;
t75 = Icges(5,1) * t109 + Icges(5,4) * t108 + Icges(5,5) * t198;
t223 = t75 + t73;
t200 = Icges(6,4) * t150;
t102 = -Icges(6,6) * t146 + (-Icges(6,2) * t148 + t200) * t145;
t202 = Icges(5,4) * t150;
t103 = -Icges(5,6) * t146 + (-Icges(5,2) * t148 + t202) * t145;
t243 = -t102 - t103;
t201 = Icges(6,4) * t148;
t104 = -Icges(6,5) * t146 + (Icges(6,1) * t150 - t201) * t145;
t203 = Icges(5,4) * t148;
t105 = -Icges(5,5) * t146 + (Icges(5,1) * t150 - t203) * t145;
t242 = -t104 - t105;
t181 = qJD(4) * t145;
t114 = (-Icges(6,2) * t150 - t201) * t181;
t115 = (-Icges(5,2) * t150 - t203) * t181;
t241 = (-t114 - t115) * t148;
t112 = (-Icges(6,5) * t148 - Icges(6,6) * t150) * t181;
t113 = (-Icges(5,5) * t148 - Icges(5,6) * t150) * t181;
t116 = (-Icges(6,1) * t148 - t200) * t181;
t117 = (-Icges(5,1) * t148 - t202) * t181;
t239 = (t116 + t117) * t145 * t150 + (-t112 - t113) * t146;
t240 = (t239 + ((t242 * t148 + t243 * t150) * qJD(4) + t241) * t145) * t146;
t238 = t226 * t108 + t223 * t109 + t228 * t198;
t237 = t225 * t108 + t224 * t109 + t227 * t198;
t236 = -t225 * t110 + t224 * t155 + t227 * t196;
t235 = -t226 * t110 + t223 * t155 + t228 * t196;
t182 = qJD(1) * t145;
t171 = t141 * t182;
t89 = -qJD(1) * t110 - qJD(4) * t109;
t90 = qJD(1) * t155 + qJD(4) * t108;
t42 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t171;
t44 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t171;
t233 = t44 + t42;
t172 = t140 * t182;
t87 = qJD(1) * t108 + qJD(4) * t155;
t154 = t110 * qJD(4);
t88 = qJD(1) * t109 + t154;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t172;
t47 = Icges(5,4) * t88 + Icges(5,2) * t87 + Icges(5,6) * t172;
t232 = t45 + t47;
t46 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t171;
t48 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t171;
t231 = t48 + t46;
t49 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t172;
t51 = Icges(5,1) * t88 + Icges(5,4) * t87 + Icges(5,5) * t172;
t230 = t49 + t51;
t50 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t171;
t52 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t171;
t229 = t52 + t50;
t138 = pkin(4) * t150 + pkin(3);
t199 = t138 * t146;
t222 = t109 * rSges(6,1) + t108 * rSges(6,2) + rSges(6,3) * t198 + t140 * t199;
t178 = pkin(4) * t197;
t179 = qJD(5) * t145;
t183 = qJD(1) * t141;
t206 = pkin(4) * qJD(4);
t221 = t90 * rSges(6,1) + t89 * rSges(6,2) + rSges(6,3) * t171 + qJD(1) * t178 + t140 * t179 - t175 * t206 + t183 * t199;
t220 = rSges(6,1) * t155 - rSges(6,2) * t110 + t178;
t218 = -rSges(6,1) * t88 - rSges(6,2) * t87 + t141 * t179;
t41 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t172;
t43 = Icges(5,5) * t88 + Icges(5,6) * t87 + Icges(5,3) * t172;
t217 = (-t43 - t41) * t141;
t147 = -qJ(5) - pkin(6);
t216 = (rSges(6,3) - t147) * t145 + t199;
t215 = 2 * m(5);
t214 = 2 * m(6);
t213 = pkin(3) * t146;
t142 = sin(qJ(1)) * pkin(1);
t143 = cos(qJ(1)) * pkin(1);
t212 = -pkin(3) + t138;
t211 = pkin(6) + t147;
t167 = t211 * t145;
t177 = pkin(4) * t195;
t210 = rSges(6,3) * t172 + pkin(4) * t154 + (-t177 + (t146 * t212 - t167) * t140) * qJD(1) - t218;
t156 = -t167 - t213;
t180 = qJD(4) * t150;
t176 = pkin(4) * t180;
t209 = (qJD(1) * t156 - t176) * t141 + t221;
t208 = t140 * t156 - t177 + t222;
t187 = pkin(6) * t196 + t141 * t213;
t194 = t145 * t147;
t207 = (t194 - t199) * t141 + t187 - rSges(6,3) * t196 - t220;
t204 = (t211 - rSges(6,3)) * t146 + (rSges(6,1) * t150 - rSges(6,2) * t148 + t212) * t145;
t170 = t148 * t181;
t188 = -(-rSges(6,1) * t148 - rSges(6,2) * t150) * t181 + pkin(4) * t170 + qJD(5) * t146;
t186 = qJ(3) * t183 + qJD(3) * t140;
t185 = t140 * pkin(2) + t142;
t184 = qJD(1) * t140;
t56 = t90 * rSges(5,1) + t89 * rSges(5,2) + rSges(5,3) * t171;
t80 = t109 * rSges(5,1) + t108 * rSges(5,2) + rSges(5,3) * t198;
t174 = pkin(2) * t183 + qJ(3) * t184 + qJD(1) * t143;
t173 = t141 * pkin(2) + t140 * qJ(3) + t143;
t166 = t204 * t145;
t163 = -rSges(5,1) * t88 - rSges(5,2) * t87;
t161 = pkin(6) * t145 + t213;
t158 = rSges(4,1) * t146 - rSges(4,2) * t145;
t82 = -rSges(5,1) * t155 + rSges(5,2) * t110 - rSges(5,3) * t196;
t152 = rSges(4,3) * t140 + t141 * t158;
t119 = (-rSges(5,1) * t148 - rSges(5,2) * t150) * t181;
t107 = -rSges(5,3) * t146 + (rSges(5,1) * t150 - rSges(5,2) * t148) * t145;
t101 = -Icges(5,3) * t146 + (Icges(5,5) * t150 - Icges(5,6) * t148) * t145;
t100 = -Icges(6,3) * t146 + (Icges(6,5) * t150 - Icges(6,6) * t148) * t145;
t92 = t152 + t173;
t91 = (-rSges(4,3) - qJ(3)) * t141 + t158 * t140 + t185;
t64 = qJD(1) * t152 - qJD(3) * t141 + t174;
t63 = (rSges(4,3) * t141 - t142 + (-pkin(2) - t158) * t140) * qJD(1) + t186;
t62 = -t107 * t196 + t146 * t82;
t61 = -t107 * t198 - t146 * t80;
t60 = -t82 + t173 + t187;
t59 = -qJ(3) * t141 + t140 * t161 + t185 + t80;
t54 = rSges(5,3) * t172 - t163;
t40 = t141 * t216 + t173 + t220;
t39 = -t140 * t194 + (-pkin(4) * t148 - qJ(3)) * t141 + t185 + t222;
t38 = -t101 * t196 + t103 * t110 - t105 * t155;
t37 = -t100 * t196 + t102 * t110 - t104 * t155;
t36 = t101 * t198 + t103 * t108 + t105 * t109;
t35 = t100 * t198 + t102 * t108 + t104 * t109;
t32 = (qJD(1) * t161 - qJD(3)) * t141 + t174 + t56;
t31 = (-t142 + (-t213 - pkin(2) + (-rSges(5,3) - pkin(6)) * t145) * t140) * qJD(1) + t163 + t186;
t30 = -t146 * t56 + (-t107 * t183 - t119 * t140) * t145;
t29 = t146 * t54 + (t107 * t184 - t119 * t141) * t145;
t28 = -t141 * t166 + t146 * t207;
t27 = -t140 * t166 - t146 * t208;
t26 = -t146 * t68 + (-t148 * t72 + t150 * t76) * t145;
t25 = -t146 * t67 + (-t148 * t71 + t150 * t75) * t145;
t24 = -t146 * t66 + (-t148 * t70 + t150 * t74) * t145;
t23 = -t146 * t65 + (-t148 * t69 + t150 * t73) * t145;
t22 = (-t147 * t182 - qJD(3) - t176) * t141 + t174 + t221;
t21 = -t110 * t206 + (t177 - t142 + (-pkin(2) - t216) * t140) * qJD(1) + t186 + t218;
t12 = t103 * t89 + t105 * t90 + t108 * t115 + t109 * t117 + (t101 * t183 + t113 * t140) * t145;
t11 = t102 * t89 + t104 * t90 + t108 * t114 + t109 * t116 + (t100 * t183 + t112 * t140) * t145;
t10 = t103 * t87 + t105 * t88 + t110 * t115 - t155 * t117 + (t101 * t184 - t113 * t141) * t145;
t9 = t102 * t87 + t104 * t88 + t110 * t114 - t155 * t116 + (t100 * t184 - t112 * t141) * t145;
t8 = -t209 * t146 + (t140 * t188 - t183 * t204) * t145;
t7 = t210 * t146 + (t141 * t188 + t184 * t204) * t145;
t6 = (t140 * t54 + t141 * t56 + (-t140 * t80 + t141 * t82) * qJD(1)) * t145;
t5 = -t146 * t43 + (-t148 * t47 + t150 * t51 + (-t148 * t76 - t150 * t72) * qJD(4)) * t145;
t4 = -t146 * t44 + (-t148 * t48 + t150 * t52 + (-t148 * t75 - t150 * t71) * qJD(4)) * t145;
t3 = -t146 * t41 + (-t148 * t45 + t150 * t49 + (-t148 * t74 - t150 * t70) * qJD(4)) * t145;
t2 = -t146 * t42 + (-t148 * t46 + t150 * t50 + (-t148 * t73 - t150 * t69) * qJD(4)) * t145;
t1 = (t209 * t141 + t210 * t140 + (-t140 * t208 + t141 * t207) * qJD(1)) * t145;
t13 = [0.2e1 * m(4) * (t63 * t92 + t64 * t91) + (t31 * t60 + t32 * t59) * t215 + (t21 * t40 + t22 * t39) * t214 + t242 * t170 + t239 + (t243 * t180 + t241) * t145; 0; 0; m(4) * (-t140 * t64 - t141 * t63 + (t140 * t92 - t141 * t91) * qJD(1)) + m(5) * (-t140 * t32 - t141 * t31 + (t140 * t60 - t141 * t59) * qJD(1)) + m(6) * (-t140 * t22 - t141 * t21 + (t140 * t40 - t141 * t39) * qJD(1)); 0; 0; -t240 + m(5) * (t29 * t60 + t30 * t59 + t31 * t62 + t32 * t61) + m(6) * (t21 * t28 + t22 * t27 + t39 * t8 + t40 * t7) + ((-t5 / 0.2e1 - t3 / 0.2e1 - t10 / 0.2e1 - t9 / 0.2e1) * t141 + (t4 / 0.2e1 + t2 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t140 + ((t25 / 0.2e1 + t23 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1) * t141 + (t26 / 0.2e1 + t24 / 0.2e1 + t38 / 0.2e1 + t37 / 0.2e1) * t140) * qJD(1)) * t145; m(5) * t6 + m(6) * t1; m(5) * (-t140 * t30 - t141 * t29 + (t140 * t62 - t141 * t61) * qJD(1)) + m(6) * (-t140 * t8 - t141 * t7 + (t140 * t28 - t141 * t27) * qJD(1)); (t62 * t29 + t61 * t30) * t215 + (t27 * t8 + t28 * t7) * t214 + ((t140 * t82 + t141 * t80) * t6 * t215 + (t140 * t207 + t141 * t208) * t1 * t214 + (-t140 * t235 + t141 * t236) * t172 + (t140 * t238 - t237 * t141) * t171 + (t237 * t184 + t238 * t183 + (-t145 * t183 * t227 - t108 * t232 - t109 * t230 - t224 * t90 - t225 * t89) * t141 + (t223 * t90 + t226 * t89 + (t140 * t233 + t183 * t228 + t217) * t145 + t229 * t109 + t231 * t108) * t140) * t198 + (t236 * t184 + t235 * t183 + (t224 * t88 + t225 * t87 + (t184 * t227 + t217) * t145 - t230 * t155 + t232 * t110) * t141 + ((t141 * t233 - t184 * t228) * t145 - t223 * t88 - t226 * t87 + t229 * t155 - t231 * t110) * t140) * t196) * t145 + (t240 + (-t12 - t11) * t198 + (t10 + t9) * t196 + (-t38 - t37) * t172 + (-t36 - t35) * t171 + ((t3 + t5) * t141 + (-t2 - t4) * t140 + ((-t23 - t25) * t141 + (-t24 - t26) * t140) * qJD(1)) * t145) * t146; m(6) * (t140 * t21 - t141 * t22 + (t140 * t39 + t141 * t40) * qJD(1)) * t145; 0; 0; m(6) * (-t1 * t146 + (t140 * t7 - t141 * t8 + (t140 * t27 + t141 * t28) * qJD(1)) * t145); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
