% Calculate time derivative of joint inertia matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:49
% DurationCPUTime: 2.84s
% Computational Cost: add. (2326->225), mult. (2848->312), div. (0->0), fcn. (2112->6), ass. (0->135)
t251 = Icges(5,5) + Icges(6,5);
t250 = Icges(5,6) + Icges(6,6);
t249 = Icges(5,3) + Icges(6,3);
t103 = sin(qJ(4));
t105 = cos(qJ(4));
t248 = t103 * t251 + t105 * t250;
t165 = Icges(6,4) * t105;
t167 = Icges(5,4) * t105;
t247 = t165 + t167 + (-Icges(5,2) - Icges(6,2)) * t103;
t166 = Icges(6,4) * t103;
t168 = Icges(5,4) * t103;
t246 = -t166 - t168 + (Icges(5,1) + Icges(6,1)) * t105;
t101 = qJ(1) + pkin(7);
t98 = sin(t101);
t99 = cos(t101);
t245 = t248 * t98 + t249 * t99;
t131 = Icges(5,1) * t103 + t167;
t195 = -Icges(5,5) * t98 + t131 * t99;
t127 = Icges(5,2) * t105 + t168;
t197 = -Icges(5,6) * t98 + t127 * t99;
t133 = -t103 * t195 - t105 * t197;
t129 = Icges(6,1) * t103 + t165;
t196 = -Icges(6,5) * t98 + t129 * t99;
t125 = Icges(6,2) * t105 + t166;
t198 = -Icges(6,6) * t98 + t125 * t99;
t135 = -t103 * t196 - t105 * t198;
t244 = t133 + t135;
t190 = rSges(6,1) + pkin(4);
t243 = t98 * t99;
t41 = Icges(6,6) * t99 + t125 * t98;
t43 = Icges(5,6) * t99 + t127 * t98;
t241 = t43 + t41;
t45 = Icges(6,5) * t99 + t129 * t98;
t47 = Icges(5,5) * t99 + t131 * t98;
t240 = -t47 - t45;
t159 = qJD(4) * t99;
t239 = t159 * t98;
t134 = t103 * t47 + t105 * t43;
t136 = t103 * t45 + t105 * t41;
t202 = t134 + t136;
t238 = t202 * t98;
t237 = t247 * qJD(4);
t236 = t246 * qJD(4);
t235 = -t103 * t250 + t105 * t251;
t233 = -t196 - t195;
t232 = t197 + t198;
t231 = -t248 * t99 + t249 * t98;
t116 = t134 * t99;
t118 = t136 * t99;
t226 = t245 * t98 - t116 - t118;
t96 = t98 ^ 2;
t97 = t99 ^ 2;
t223 = t96 + t97;
t222 = t244 * t98;
t155 = qJD(4) * t105;
t161 = qJD(1) * t99;
t221 = t103 * t161 + t155 * t98;
t156 = qJD(4) * t103;
t147 = t98 * t156;
t157 = qJD(1) * t105;
t148 = t99 * t157;
t220 = t147 - t148;
t219 = t245 * qJD(1);
t173 = t103 * rSges(6,2);
t218 = qJD(5) * t98 - (t105 * t190 - t173) * t159;
t217 = t233 * t98 + t240 * t99;
t216 = t232 * t98 + t241 * t99;
t215 = -t245 * t99 - t238;
t214 = -t231 * t99 - t222;
t102 = -qJ(5) - pkin(6);
t162 = qJD(1) * t98;
t107 = t220 * rSges(6,2) - qJD(5) * t99 - t102 * t162 - t190 * t221;
t171 = t105 * rSges(6,2);
t120 = t103 * t190 + t171;
t169 = t102 - rSges(6,3);
t141 = -rSges(6,1) * t103 - t171;
t184 = -pkin(6) - t102;
t182 = (-pkin(4) * t103 + t141) * t99 + (rSges(6,3) + t184) * t98;
t170 = t105 * t98;
t172 = t103 * t98;
t203 = -rSges(6,2) * t170 - t99 * rSges(6,3) - t172 * t190;
t183 = -t184 * t99 + t203;
t1 = (((rSges(6,3) - pkin(6)) * t98 - t182) * qJD(1) + t107) * t98 + (((-pkin(6) - t169) * t99 + t120 * t98 + t183) * qJD(1) + t218) * t99;
t192 = 2 * m(6);
t213 = t1 * t192;
t193 = 2 * m(5);
t142 = rSges(5,1) * t103 + rSges(5,2) * t105;
t119 = t99 * t142;
t152 = rSges(5,1) * t221 + rSges(5,2) * t148;
t180 = rSges(5,2) * t103;
t185 = t98 * rSges(5,3);
t94 = t99 * rSges(5,3);
t50 = rSges(5,1) * t172 + rSges(5,2) * t170 + t94;
t52 = t185 - t119;
t81 = rSges(5,1) * t105 - t180;
t2 = -t98 * t152 + (t180 * t96 - t81 * t97) * qJD(4) + ((-t50 + t94) * t99 + (t119 - t52 + t185) * t98) * qJD(1);
t212 = t193 * t2;
t210 = t231 * t98 - t244 * t99;
t209 = -t159 * t235 + t219;
t208 = qJD(4) * t235 * t98 - qJD(1) * t231;
t154 = -rSges(5,3) - pkin(2) - pkin(6);
t188 = sin(qJ(1)) * pkin(1);
t114 = t154 * t98 - t188;
t181 = qJ(3) * t161 + qJD(3) * t98;
t13 = -rSges(5,2) * t147 + qJD(1) * t114 + t152 + t181;
t100 = cos(qJ(1)) * pkin(1);
t90 = qJD(3) * t99;
t14 = t90 + t81 * t159 + (-t100 + t154 * t99 + (-qJ(3) - t142) * t98) * qJD(1);
t201 = -t13 * t99 + t14 * t98;
t191 = m(5) * t81;
t189 = rSges(4,2) - pkin(2);
t153 = -pkin(2) + t169;
t151 = pkin(2) * t99 + qJ(3) * t98 + t100;
t64 = t142 * qJD(4);
t145 = t223 * t64;
t80 = rSges(6,1) * t105 - t173;
t144 = pkin(4) * t105 + t80;
t108 = rSges(4,3) * t99 + t189 * t98 - t188;
t92 = t99 * qJ(3);
t63 = t141 * qJD(4);
t54 = t144 * t99;
t53 = t144 * t98;
t36 = -rSges(4,2) * t99 + rSges(4,3) * t98 + t151;
t35 = t92 + t108;
t34 = t90 + (-t100 + t189 * t99 + (-rSges(4,3) - qJ(3)) * t98) * qJD(1);
t33 = qJD(1) * t108 + t181;
t32 = pkin(6) * t99 + t151 + t50;
t31 = t92 + t119 + t114;
t18 = -pkin(4) * t220 + t161 * t80 + t63 * t98;
t17 = t80 * t162 - t63 * t99 + (t156 * t99 + t157 * t98) * pkin(4);
t16 = -t102 * t99 + t151 - t203;
t15 = t120 * t99 + t153 * t98 - t188 + t92;
t4 = t90 + (-t100 + t153 * t99 + (-qJ(3) - t120) * t98) * qJD(1) - t218;
t3 = (-t188 + (-rSges(6,3) - pkin(2)) * t98) * qJD(1) - t107 + t181;
t5 = [0.2e1 * m(4) * (t33 * t36 + t34 * t35) + (t13 * t32 + t14 * t31) * t193 + (t15 * t4 + t16 * t3) * t192 + (t127 + t125) * t156 + (-t131 - t129) * t155 - t237 * t105 - t236 * t103; 0; 0; m(4) * (-t33 * t99 + t34 * t98 + (t35 * t99 + t36 * t98) * qJD(1)) + m(5) * ((t31 * t99 + t32 * t98) * qJD(1) + t201) + m(6) * (-t3 * t99 + t4 * t98 + (t15 * t99 + t16 * t98) * qJD(1)); 0; 0; m(5) * (t201 * t81 - (t31 * t98 - t32 * t99) * t64) + m(6) * (t15 * t18 + t16 * t17 - t3 * t54 + t4 * t53) + (t236 * t243 / 0.2e1 - t246 * t239 / 0.2e1 + (-t233 * t99 - t240 * t98) * qJD(1) / 0.2e1) * t105 + (-t237 * t243 / 0.2e1 + t247 * t239 / 0.2e1 - (t232 * t99 + t241 * t98) * qJD(1) / 0.2e1) * t103 + ((t32 * t191 + (-t47 / 0.2e1 - t45 / 0.2e1) * t105 + (t43 / 0.2e1 + t41 / 0.2e1) * t103) * t98 + (t31 * t191 + (-t195 / 0.2e1 - t196 / 0.2e1) * t105 + (t197 / 0.2e1 + t198 / 0.2e1) * t103) * t99) * qJD(1) + (-t116 / 0.2e1 - t118 / 0.2e1 + (-t133 / 0.2e1 - t135 / 0.2e1) * t98 - t223 * t248 / 0.2e1) * qJD(4); m(5) * t2 + m(6) * t1; m(6) * (-t17 * t99 + t18 * t98 + (t53 * t99 - t54 * t98) * qJD(1)) - m(5) * t145; -t81 * t145 * t193 + (-t54 * t17 + t53 * t18) * t192 + (t182 * t213 + t52 * t212 + t208 * t97 + (t202 * t99 - t214) * t161) * t99 + (-t50 * t212 + t183 * t213 + (t209 * t98 + (t222 - t226) * qJD(1)) * t98 + (t208 * t98 + (t209 + t219) * t99 + t216 * t156 + t217 * t155 + (-t103 * t216 - t105 * t217) * qJD(4) + (t238 + (t244 - t245) * t99 + t210 + t215) * qJD(1)) * t99) * t98 + (t214 * t98 + t215 * t99) * t162 + (t210 * t98 + t226 * t99) * t161; m(6) * (t3 * t98 + t4 * t99 + (-t15 * t98 + t16 * t99) * qJD(1)); 0; 0; m(6) * (t17 * t98 + t18 * t99 + (-t53 * t98 - t54 * t99) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
