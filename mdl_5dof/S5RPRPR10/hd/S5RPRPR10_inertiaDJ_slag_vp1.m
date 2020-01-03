% Calculate time derivative of joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:58
% DurationCPUTime: 1.80s
% Computational Cost: add. (4391->229), mult. (5614->319), div. (0->0), fcn. (5344->8), ass. (0->129)
t114 = sin(qJ(5));
t180 = qJD(5) * t114;
t196 = cos(qJ(1));
t112 = t196 * pkin(2);
t195 = sin(qJ(1));
t181 = t196 * pkin(1) + t195 * qJ(2);
t219 = -t112 - t181;
t116 = cos(qJ(5));
t217 = qJD(1) - qJD(3);
t178 = qJ(3) + pkin(8);
t159 = sin(t178);
t160 = cos(t178);
t68 = -t195 * t159 - t196 * t160;
t54 = t217 * t68;
t69 = t196 * t159 - t195 * t160;
t140 = -t116 * t54 + t69 * t180;
t55 = t217 * t69;
t138 = t116 * t55 + t68 * t180;
t218 = -rSges(3,1) * t196 - rSges(3,3) * t195 - t181;
t115 = sin(qJ(3));
t168 = t196 * qJD(1);
t169 = qJD(3) * t196;
t117 = cos(qJ(3));
t171 = t195 * t117;
t216 = qJD(3) * t171 + (t168 - t169) * t115;
t179 = qJD(5) * t116;
t139 = -t114 * t55 + t179 * t68;
t141 = t114 * t54 + t179 * t69;
t201 = 2 * m(6);
t189 = rSges(6,2) * t114;
t190 = rSges(6,1) * t116;
t77 = (t189 - t190) * qJD(5);
t88 = -rSges(6,1) * t114 - rSges(6,2) * t116;
t215 = -(Icges(6,5) * t140 + Icges(6,6) * t141 + Icges(6,3) * t55) * t68 + (Icges(6,5) * t138 + Icges(6,6) * t139 + Icges(6,3) * t54) * t69 + t88 * t77 * t201;
t183 = Icges(6,4) * t116;
t145 = Icges(6,2) * t114 - t183;
t40 = -Icges(6,6) * t68 + t145 * t69;
t184 = Icges(6,4) * t114;
t146 = -Icges(6,1) * t116 + t184;
t42 = -Icges(6,5) * t68 + t146 * t69;
t149 = t114 * t40 - t116 * t42;
t212 = t149 * t68;
t41 = Icges(6,6) * t69 + t145 * t68;
t43 = Icges(6,5) * t69 + t146 * t68;
t148 = t114 * t41 - t116 * t43;
t213 = t148 * t69;
t144 = -Icges(6,5) * t116 + Icges(6,6) * t114;
t38 = -Icges(6,3) * t68 + t144 * t69;
t39 = Icges(6,3) * t69 + t144 * t68;
t214 = t38 * t69 - t39 * t68 + t212 + t213;
t109 = qJD(2) * t196;
t211 = -qJD(1) * t181 + t109;
t167 = t195 * qJD(1);
t182 = qJ(2) * t168 + qJD(2) * t195;
t210 = -pkin(1) * t167 + t182;
t132 = t138 * rSges(6,1) + t54 * rSges(6,3);
t194 = pkin(3) * t117;
t107 = pkin(2) + t194;
t164 = t216 * pkin(3) - t107 * t167;
t209 = -t139 * rSges(6,2) - t55 * pkin(4) - t54 * pkin(7) - t132 - t164;
t133 = t140 * rSges(6,1) + t55 * rSges(6,3);
t172 = t195 * t115;
t102 = pkin(3) * t172;
t165 = -t217 * t102 - t107 * t168 + t169 * t194;
t208 = t141 * rSges(6,2) - t54 * pkin(4) + t55 * pkin(7) + t133 - t165;
t170 = -pkin(4) + t189;
t173 = t196 * t115;
t185 = -pkin(3) * t173 + t195 * t107;
t192 = t68 * rSges(6,3) + t69 * t190;
t207 = t68 * pkin(7) - t170 * t69 - t185 + t192;
t186 = t196 * t107 + t102;
t191 = t69 * rSges(6,3) - t68 * t190;
t206 = -pkin(7) * t69 - t170 * t68 - t186 - t191;
t166 = (Icges(6,2) * t116 + t146 + t184) * t180;
t73 = t145 * qJD(5);
t83 = -Icges(6,1) * t114 - t183;
t205 = -(qJD(5) * t83 + t73) * t116 - t166;
t203 = 2 * m(4);
t202 = 2 * m(5);
t177 = t195 * pkin(1);
t176 = t195 * pkin(2);
t158 = pkin(2) * t168;
t156 = pkin(2) * t167;
t111 = t196 * qJ(2);
t152 = t111 - t177;
t126 = t117 * t196 + t172;
t56 = t217 * t126;
t57 = -t117 * t167 + t216;
t35 = -t57 * rSges(4,1) - t56 * rSges(4,2);
t34 = t56 * rSges(4,1) - t57 * rSges(4,2);
t76 = -t173 + t171;
t58 = t76 * rSges(4,1) - rSges(4,2) * t126;
t59 = -rSges(4,1) * t126 - rSges(4,2) * t76;
t131 = rSges(5,1) * t68 + rSges(5,2) * t69 - t186;
t130 = -t69 * rSges(5,1) + t68 * rSges(5,2) + t185;
t129 = -t177 - t176;
t125 = -t55 * rSges(5,1) + t54 * rSges(5,2) - t164;
t124 = -t54 * rSges(5,1) - t55 * rSges(5,2) - t165;
t122 = -rSges(3,1) * t195 + rSges(3,3) * t196 - t177;
t72 = t144 * qJD(5);
t118 = (-t114 * t43 - t116 * t41) * t54 / 0.2e1 + (-t114 * t42 - t116 * t40) * t55 / 0.2e1 - (qJD(5) * t149 - t114 * (Icges(6,1) * t140 + Icges(6,4) * t141 + Icges(6,5) * t55) - t116 * (Icges(6,4) * t140 + Icges(6,2) * t141 + Icges(6,6) * t55) - t68 * t72) * t68 / 0.2e1 + (qJD(5) * t148 - t114 * (Icges(6,1) * t138 + Icges(6,4) * t139 + Icges(6,5) * t54) - t116 * (Icges(6,4) * t138 + Icges(6,2) * t139 + Icges(6,6) * t54) + t69 * t72) * t69 / 0.2e1 + (t69 * t54 - t68 * t55) * (-Icges(6,5) * t114 - Icges(6,6) * t116);
t80 = t88 ^ 2;
t66 = t111 + t122;
t61 = t218 * qJD(1) + t109;
t60 = qJD(1) * t122 + t182;
t51 = -t59 - t219;
t50 = t111 + t129 - t58;
t47 = t112 + t131;
t46 = -t176 + t130;
t45 = t189 * t68 + t191;
t44 = t189 * t69 - t192;
t37 = -t131 + t181;
t36 = -t130 + t152;
t33 = t219 * qJD(1) + t109 - t34;
t32 = qJD(1) * t129 + t182 - t35;
t29 = t112 + t206;
t28 = -t176 - t207;
t27 = -t158 + t124;
t26 = -t156 + t125;
t25 = t181 - t206;
t24 = t152 + t207;
t21 = -t124 + t211;
t20 = -t125 + t210;
t9 = -t158 + t208;
t8 = -t156 + t209;
t7 = -t208 + t211;
t6 = -t209 + t210;
t1 = t54 * t44 + t69 * t133 - t55 * t45 + t68 * t132 + (t139 * t68 + t141 * t69) * rSges(6,2);
t2 = [(t24 * t7 + t25 * t6) * t201 - t83 * t179 - t116 * t73 + (t20 * t37 + t21 * t36) * t202 + (t32 * t51 + t33 * t50) * t203 + 0.2e1 * m(3) * (-t218 * t60 + t61 * t66) - t166; m(6) * (t195 * t7 - t196 * t6 + (t195 * t25 + t196 * t24) * qJD(1)) + m(5) * (t195 * t21 - t196 * t20 + (t195 * t37 + t196 * t36) * qJD(1)) + m(4) * (t195 * t33 - t196 * t32 + (t195 * t51 + t196 * t50) * qJD(1)) + m(3) * (t195 * t61 - t196 * t60 + (-t195 * t218 + t196 * t66) * qJD(1)); 0; m(6) * (t24 * t9 + t25 * t8 + t28 * t7 + t29 * t6) + m(5) * (t20 * t47 + t21 * t46 + t26 * t37 + t27 * t36) + m(4) * (t32 * t59 + t33 * t58 + t34 * t50 + t35 * t51) - t205; m(4) * (t34 * t195 - t35 * t196 + (t195 * t59 + t196 * t58) * qJD(1)) + m(5) * (t27 * t195 - t26 * t196 + (t195 * t47 + t196 * t46) * qJD(1)) + m(6) * (t9 * t195 - t8 * t196 + (t195 * t29 + t196 * t28) * qJD(1)); (t28 * t9 + t29 * t8) * t201 + (t26 * t47 + t27 * t46) * t202 + (t34 * t58 + t35 * t59) * t203 + t205; 0; 0; 0; 0; m(6) * ((-t24 * t68 - t25 * t69) * t77 + (t24 * t55 - t25 * t54 - t6 * t69 - t68 * t7) * t88) + t118; m(6) * ((-t195 * t68 + t196 * t69) * t77 + (t195 * t55 + t196 * t54 + (-t195 * t69 - t196 * t68) * qJD(1)) * t88); m(6) * ((-t28 * t68 - t29 * t69) * t77 + (t28 * t55 - t29 * t54 - t68 * t9 - t69 * t8) * t88) - t118; -m(6) * t1; ((t44 * t1 + t80 * t54) * t201 + (-t213 + t214) * t55 + (0.3e1 * t54 * t39 + t215) * t69) * t69 + ((-t149 - t39) * t55 * t69 + (0.3e1 * t55 * t38 + t215) * t68 - (-t212 + t214 + (t38 - t148) * t69) * t54 + (t45 * t1 - t80 * t55) * t201) * t68;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
