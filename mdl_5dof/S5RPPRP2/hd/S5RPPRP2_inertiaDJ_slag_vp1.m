% Calculate time derivative of joint inertia matrix for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:11
% DurationCPUTime: 2.64s
% Computational Cost: add. (3535->212), mult. (3162->288), div. (0->0), fcn. (2381->8), ass. (0->130)
t233 = Icges(6,4) + Icges(5,5);
t232 = -Icges(5,6) + Icges(6,6);
t231 = Icges(6,2) + Icges(5,3);
t101 = pkin(8) + qJ(4);
t96 = sin(t101);
t98 = cos(t101);
t229 = t232 * t96 + t233 * t98;
t190 = rSges(6,1) + pkin(4);
t230 = t190 * t98;
t102 = qJ(1) + pkin(7);
t97 = sin(t102);
t99 = cos(t102);
t226 = t229 * t99 + t231 * t97;
t171 = Icges(6,5) * t98;
t128 = Icges(6,3) * t96 + t171;
t41 = Icges(6,6) * t97 + t128 * t99;
t175 = Icges(5,4) * t98;
t132 = -Icges(5,2) * t96 + t175;
t47 = Icges(5,6) * t97 + t132 * t99;
t228 = t41 - t47;
t227 = t229 * t97 - t231 * t99;
t40 = -Icges(6,6) * t99 + t128 * t97;
t46 = -Icges(5,6) * t99 + t132 * t97;
t225 = t46 - t40;
t172 = Icges(6,5) * t96;
t134 = Icges(6,1) * t98 + t172;
t48 = -Icges(6,4) * t99 + t134 * t97;
t176 = Icges(5,4) * t96;
t136 = Icges(5,1) * t98 - t176;
t50 = -Icges(5,5) * t99 + t136 * t97;
t224 = t48 + t50;
t49 = Icges(6,4) * t97 + t134 * t99;
t51 = Icges(5,5) * t97 + t136 * t99;
t223 = t51 + t49;
t137 = t47 * t96 - t51 * t98;
t139 = t41 * t96 + t49 * t98;
t198 = t137 - t139;
t222 = t198 * t99;
t157 = qJD(4) * t98;
t163 = rSges(6,3) + qJ(5);
t221 = t163 * t157;
t220 = (t232 * t98 - t233 * t96) * qJD(4);
t138 = t46 * t96 - t50 * t98;
t140 = t40 * t96 + t48 * t98;
t219 = t138 - t140;
t179 = -t163 * t98 + t190 * t96;
t166 = qJ(5) * t96;
t216 = rSges(6,3) * t96 + t166 + t230;
t215 = t226 * qJD(1);
t120 = t137 * t97;
t121 = t138 * t99;
t122 = t139 * t97;
t123 = t140 * t99;
t214 = t226 * t99 - t227 * t97 + t120 + t121 - t122 - t123;
t213 = t223 * t97 + t224 * t99;
t212 = t225 * t99 - t228 * t97;
t94 = t97 ^ 2;
t95 = t99 ^ 2;
t211 = t219 * t97 + t227 * t99;
t195 = 2 * m(5);
t75 = rSges(5,1) * t96 + rSges(5,2) * t98;
t119 = t75 * qJD(4);
t159 = qJD(1) * t99;
t160 = qJD(1) * t97;
t188 = rSges(5,2) * t96;
t109 = rSges(5,3) * t159 - t119 * t99 + t160 * t188;
t184 = t96 * t99;
t91 = t97 * rSges(5,3);
t199 = -rSges(5,2) * t184 + t91;
t200 = t97 * t119;
t189 = rSges(5,1) * t98;
t149 = -t188 + t189;
t53 = -rSges(5,3) * t99 + t149 * t97;
t183 = t98 * t99;
t55 = rSges(5,1) * t183 + t199;
t2 = (qJD(1) * t53 + t109) * t99 + (-t200 + (-t55 + t199) * qJD(1)) * t97;
t209 = t195 * t2;
t207 = t226 * t97 - t222;
t206 = -t227 * qJD(1) + t220 * t99;
t205 = -t220 * t97 - t215;
t105 = -pkin(6) - qJ(3);
t104 = cos(pkin(8));
t93 = pkin(3) * t104 + pkin(2);
t126 = -t149 - t93;
t187 = sin(qJ(1)) * pkin(1);
t34 = -t187 + (rSges(5,3) - t105) * t99 + t126 * t97;
t100 = cos(qJ(1)) * pkin(1);
t151 = -t105 * t97 + t99 * t93 + t100;
t35 = t151 + t55;
t196 = t34 * t99 + t35 * t97;
t125 = rSges(4,1) * t104 - rSges(4,2) * sin(pkin(8)) + pkin(2);
t164 = rSges(4,3) + qJ(3);
t37 = t125 * t99 + t164 * t97 + t100;
t194 = 2 * m(6);
t191 = m(5) * t75;
t92 = t97 * rSges(6,2);
t182 = -rSges(6,2) * t99 + t216 * t97;
t181 = rSges(6,3) * t184 + t99 * t166 + t183 * t190 + t92;
t180 = -qJD(4) * t216 + qJD(5) * t98;
t90 = qJD(3) * t99;
t178 = t105 * t160 + t90;
t177 = t94 + t95;
t165 = t105 * t99;
t158 = qJD(4) * t96;
t156 = qJD(5) * t96;
t154 = m(6) * t158;
t39 = t179 * t99;
t152 = t190 * qJD(4);
t111 = -t163 * t96 - t230 - t93;
t108 = t111 * t97 - t187;
t18 = (rSges(6,2) - t105) * t99 + t108;
t19 = t151 + t181;
t146 = t18 * t99 + t19 * t97;
t38 = t179 * t97;
t141 = -t38 * t97 - t39 * t99;
t110 = rSges(6,2) * t159 - t152 * t184 + (t156 + t221) * t99;
t36 = -t125 * t97 + t164 * t99 - t187;
t89 = qJD(3) * t97;
t66 = t149 * qJD(4);
t33 = -qJD(1) * t37 + t90;
t32 = qJD(1) * t36 + t89;
t17 = -qJD(1) * t39 + t180 * t97;
t16 = t160 * t179 + t180 * t99;
t15 = t200 + (t126 * t99 - t100 - t91) * qJD(1) + t178;
t14 = t89 + (-t187 - t165 + (-t93 - t189) * t97) * qJD(1) + t109;
t5 = t181 * t99 + t182 * t97;
t4 = (t179 * qJD(4) - t156) * t97 + (t111 * t99 - t100 - t92) * qJD(1) + t178;
t3 = t89 + (t108 - t165) * qJD(1) + t110;
t1 = (qJD(1) * t182 + t110) * t99 + ((-t181 + t92) * qJD(1) + (t221 + (qJD(5) - t152) * t96) * t97) * t97;
t6 = [0.2e1 * m(4) * (t32 * t37 + t33 * t36) + (t14 * t35 + t15 * t34) * t195 + (t18 * t4 + t19 * t3) * t194 + (t134 + t136 + t172 - t176 + (-Icges(5,2) - Icges(6,3)) * t98) * t158 + (-t128 + t132 - t171 + t175 + (Icges(5,1) + Icges(6,1)) * t96) * t157; 0; 0; m(4) * (-t32 * t99 + t33 * t97 + (t36 * t99 + t37 * t97) * qJD(1)) + m(5) * (t196 * qJD(1) - t14 * t99 + t15 * t97) + m(6) * (qJD(1) * t146 - t3 * t99 + t4 * t97); 0; 0; m(5) * ((-t14 * t97 - t15 * t99) * t75 - t196 * t66) + m(6) * (t16 * t18 + t17 * t19 - t3 * t38 - t39 * t4) + (-t120 / 0.2e1 + t122 / 0.2e1 + t121 / 0.2e1 - t123 / 0.2e1 + t229 * (t94 / 0.2e1 + t95 / 0.2e1)) * qJD(4) + ((-t35 * t191 + (t47 / 0.2e1 - t41 / 0.2e1) * t98 + (t51 / 0.2e1 + t49 / 0.2e1) * t96) * t99 + (t34 * t191 + (t46 / 0.2e1 - t40 / 0.2e1) * t98 + (t50 / 0.2e1 + t48 / 0.2e1) * t96) * t97 + (-t224 * t96 - t225 * t98) * t97 / 0.2e1 - (t223 * t96 - t228 * t98) * t99 / 0.2e1) * qJD(1); m(5) * t2 + m(6) * t1; m(6) * (qJD(1) * t141 + t16 * t97 - t17 * t99); t177 * t75 * t66 * t195 + (t1 * t5 - t16 * t39 - t17 * t38) * t194 + (t211 * t160 + t205 * t95 + t55 * t209 + (-t219 * t99 + t214) * t159) * t99 + (t53 * t209 + t206 * t94 + t207 * t159 + ((t205 - t215) * t97 + t206 * t99 + t213 * t158 + t212 * t157 + (-t212 * t98 - t213 * t96) * qJD(4) + ((-t219 + t226) * t97 + t222 + t207 + t211) * qJD(1)) * t99 + (t198 * t97 - t214) * t160) * t97; m(6) * (t146 * t157 + (t3 * t97 + t4 * t99 + (-t18 * t97 + t19 * t99) * qJD(1)) * t96); t154; 0; m(6) * ((qJD(4) * t141 - t1) * t98 + (qJD(4) * t5 + t16 * t99 + t17 * t97 + (-t38 * t99 + t39 * t97) * qJD(1)) * t96); 0.2e1 * (-0.1e1 + t177) * t98 * t154;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
