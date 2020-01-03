% Calculate time derivative of joint inertia matrix for
% S4RPRP7
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:08
% DurationCPUTime: 3.72s
% Computational Cost: add. (1124->238), mult. (2944->361), div. (0->0), fcn. (2218->4), ass. (0->138)
t102 = cos(qJ(1));
t101 = cos(qJ(3));
t173 = rSges(5,3) + qJ(4);
t142 = t173 * t101;
t192 = rSges(5,1) + pkin(3);
t99 = sin(qJ(3));
t148 = t192 * t99;
t210 = t148 - t142;
t221 = t102 * t210;
t179 = Icges(5,5) * t99;
t180 = Icges(4,4) * t99;
t220 = t179 - t180 + (Icges(4,1) + Icges(5,1)) * t101;
t155 = qJD(4) * t101;
t213 = t173 * t99;
t186 = t101 * t192 + t213;
t219 = (t186 * qJD(3) - t155) * t102;
t218 = t220 * qJD(3);
t100 = sin(qJ(1));
t217 = t100 / 0.2e1;
t156 = qJD(3) * t102;
t216 = -t156 / 0.2e1;
t215 = -qJD(1) / 0.2e1;
t214 = qJD(1) / 0.2e1;
t97 = t100 ^ 2;
t98 = t102 ^ 2;
t182 = t97 + t98;
t157 = qJD(3) * t101;
t159 = qJD(1) * t102;
t211 = t100 * t157 + t99 * t159;
t175 = t100 * t99;
t94 = t102 * rSges(5,2);
t207 = -t192 * t175 - t94;
t181 = rSges(4,2) * t101;
t139 = rSges(4,1) * t99 + t181;
t112 = t102 * t139;
t121 = Icges(4,2) * t101 + t180;
t43 = Icges(4,6) * t102 + t100 * t121;
t172 = Icges(4,4) * t101;
t125 = Icges(4,1) * t99 + t172;
t47 = Icges(4,5) * t102 + t100 * t125;
t130 = t101 * t43 + t47 * t99;
t109 = t130 * t102;
t115 = -Icges(5,3) * t101 + t179;
t37 = Icges(5,6) * t102 + t100 * t115;
t169 = Icges(5,5) * t101;
t123 = Icges(5,1) * t99 - t169;
t45 = Icges(5,4) * t102 + t100 * t123;
t132 = t101 * t37 - t45 * t99;
t111 = t132 * t102;
t152 = rSges(4,1) * t211 + t159 * t181;
t194 = -pkin(1) - pkin(5);
t153 = -rSges(4,3) + t194;
t162 = qJD(3) * t99;
t184 = qJ(2) * t159 + qJD(2) * t100;
t15 = (-rSges(4,2) * t162 + qJD(1) * t153) * t100 + t152 + t184;
t190 = rSges(4,2) * t99;
t74 = rSges(4,1) * t101 - t190;
t90 = qJD(2) * t102;
t16 = t90 + t74 * t156 + (t153 * t102 + (-qJ(2) - t139) * t100) * qJD(1);
t206 = t100 * t16 - t102 * t15;
t205 = -Icges(5,6) * t100 + t102 * t115;
t117 = Icges(4,5) * t99 + Icges(4,6) * t101;
t204 = -Icges(4,3) * t100 + t102 * t117;
t119 = Icges(5,4) * t99 - Icges(5,6) * t101;
t203 = -Icges(5,2) * t100 + t102 * t119;
t202 = -Icges(4,6) * t100 + t102 * t121;
t201 = -Icges(5,4) * t100 + t102 * t123;
t200 = -Icges(4,5) * t100 + t102 * t125;
t143 = qJD(1) * t186;
t187 = -qJD(3) * t210 + qJD(4) * t99;
t13 = t100 * t143 - t102 * t187;
t14 = t100 * t187 + t102 * t143;
t35 = t186 * t100;
t36 = t186 * t102;
t199 = qJD(1) * (-t100 * t36 + t102 * t35) + t14 * t100 - t13 * t102;
t154 = -rSges(5,2) + t194;
t141 = t154 * t100;
t92 = t102 * qJ(2);
t17 = t141 + t92 + t221;
t183 = t102 * pkin(1) + t100 * qJ(2);
t151 = t102 * pkin(5) + t183;
t18 = -t100 * t142 + t151 - t207;
t158 = qJD(3) * t100;
t103 = t100 * t155 - t158 * t213 - t192 * t211;
t2 = (-t102 * t142 + t141) * qJD(1) - t103 + t184;
t3 = t90 + t219 + (t154 * t102 + (-qJ(2) - t210) * t100) * qJD(1);
t198 = qJD(1) * (t100 * t18 + t102 * t17) + t100 * t3 - t102 * t2;
t196 = 2 * m(4);
t195 = 2 * m(5);
t193 = m(4) * t74;
t191 = rSges(3,2) - pkin(1);
t161 = t100 * t101;
t189 = t161 * t173 + t207;
t178 = t100 * rSges(5,2);
t188 = t178 - t221;
t177 = t100 * rSges(4,3);
t93 = t102 * rSges(4,3);
t39 = Icges(4,3) * t102 + t100 * t117;
t164 = qJD(1) * t39;
t41 = Icges(5,2) * t102 + t100 * t119;
t163 = qJD(1) * t41;
t50 = rSges(4,1) * t175 + rSges(4,2) * t161 + t93;
t150 = m(5) * t162;
t65 = t139 * qJD(3);
t145 = t182 * t65;
t131 = -t101 * t205 + t201 * t99;
t129 = -t101 * t202 - t200 * t99;
t120 = -Icges(4,2) * t99 + t172;
t118 = Icges(5,4) * t101 + Icges(5,6) * t99;
t116 = Icges(4,5) * t101 - Icges(4,6) * t99;
t114 = Icges(5,3) * t99 + t169;
t113 = rSges(3,3) * t102 + t100 * t191;
t110 = t131 * t100;
t108 = t129 * t100;
t105 = qJD(3) * t120;
t104 = qJD(3) * t114;
t54 = -rSges(3,2) * t102 + t100 * rSges(3,3) + t183;
t53 = t92 + t113;
t52 = t177 - t112;
t34 = t90 + (t191 * t102 + (-rSges(3,3) - qJ(2)) * t100) * qJD(1);
t33 = qJD(1) * t113 + t184;
t32 = t151 + t50;
t31 = t100 * t153 + t112 + t92;
t24 = t203 * qJD(1) + t118 * t158;
t23 = -t118 * t156 + t163;
t22 = t204 * qJD(1) + t116 * t158;
t21 = -t116 * t156 + t164;
t12 = -t100 * t204 - t102 * t129;
t11 = t100 * t39 - t109;
t10 = -t100 * t203 + t102 * t131;
t9 = t100 * t41 + t111;
t8 = -t102 * t204 + t108;
t7 = t130 * t100 + t102 * t39;
t6 = -t102 * t203 - t110;
t5 = -t132 * t100 + t102 * t41;
t4 = t100 * t189 + t102 * t188;
t1 = ((t178 - t188) * qJD(1) + t103) * t100 + (-t219 + (t100 * t148 + t189 + t94) * qJD(1)) * t102;
t19 = [0.2e1 * m(3) * (t33 * t54 + t53 * t34) + (t15 * t32 + t16 * t31) * t196 + (t17 * t3 + t18 * t2) * t195 - t218 * t99 + (t121 - t115) * t162 + (-t125 - t123) * t157 + (-t105 + t104) * t101; m(3) * (t100 * t34 - t102 * t33 + (t100 * t54 + t102 * t53) * qJD(1)) + m(4) * ((t100 * t32 + t102 * t31) * qJD(1) + t206) + m(5) * t198; 0; m(4) * (t206 * t74 - (t100 * t31 - t102 * t32) * t65) + m(5) * (t13 * t18 + t14 * t17 - t2 * t36 + t3 * t35) + ((t202 * t215 - t100 * t105 / 0.2e1 + t205 * t214 + t104 * t217) * t102 + (t43 * t215 + t120 * t156 / 0.2e1 + t37 * t214 + t114 * t216) * t100) * t99 + (t218 * t102 * t217 + t220 * t100 * t216 + ((t200 + t201) * t102 + (t45 + t47) * t100) * t214) * t101 + ((t32 * t193 + (t43 / 0.2e1 - t37 / 0.2e1) * t99 + (-t47 / 0.2e1 - t45 / 0.2e1) * t101) * t100 + (t31 * t193 + (t202 / 0.2e1 - t205 / 0.2e1) * t99 + (-t200 / 0.2e1 - t201 / 0.2e1) * t101) * t102) * qJD(1) + (-t109 / 0.2e1 + t111 / 0.2e1 + (-t129 / 0.2e1 + t131 / 0.2e1) * t100 - t182 * (t117 + t119) / 0.2e1) * qJD(3); -m(4) * t145 + m(5) * t199; ((-t100 * t50 + t102 * t52) * (-t100 * t152 + (t190 * t97 - t74 * t98) * qJD(3) + ((-t50 + t93) * t102 + (t112 - t52 + t177) * t100) * qJD(1)) - t74 * t145) * t196 + t102 * ((t102 * t22 + (t8 + t109) * qJD(1)) * t102 + (-t7 * qJD(1) + (-t157 * t200 + t162 * t202) * t100 + (t21 + (t101 * t47 - t43 * t99) * qJD(3) + (t129 - t39) * qJD(1)) * t102) * t100) + t100 * ((t100 * t21 + (-t11 + t108) * qJD(1)) * t100 + (t12 * qJD(1) + (-t157 * t47 + t162 * t43 + t164) * t102 + (t22 + (t101 * t200 - t202 * t99) * qJD(3) + t130 * qJD(1)) * t100) * t102) + (t1 * t4 - t13 * t36 + t14 * t35) * t195 + t102 * ((t102 * t24 + (t6 - t111) * qJD(1)) * t102 + (-t5 * qJD(1) + (-t157 * t201 - t162 * t205) * t100 + (t23 + (t101 * t45 + t37 * t99) * qJD(3) + (-t131 - t41) * qJD(1)) * t102) * t100) + t100 * ((t100 * t23 + (-t9 - t110) * qJD(1)) * t100 + (t10 * qJD(1) + (-t157 * t45 - t162 * t37 + t163) * t102 + (t24 + (t101 * t201 + t205 * t99) * qJD(3) - t132 * qJD(1)) * t100) * t102) + ((-t5 - t7) * t102 + (-t6 - t8) * t100) * qJD(1) * t100 + ((t11 + t9) * t102 + (t10 + t12) * t100) * t159; m(5) * ((t100 * t17 - t102 * t18) * t162 - t198 * t101); t182 * t150; m(5) * ((t1 + (t100 * t35 + t102 * t36) * qJD(3)) * t99 + (qJD(3) * t4 - t199) * t101); 0.2e1 * (0.1e1 - t182) * t101 * t150;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t19(1), t19(2), t19(4), t19(7); t19(2), t19(3), t19(5), t19(8); t19(4), t19(5), t19(6), t19(9); t19(7), t19(8), t19(9), t19(10);];
Mq = res;
