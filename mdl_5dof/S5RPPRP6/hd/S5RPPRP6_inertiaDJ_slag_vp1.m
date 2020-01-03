% Calculate time derivative of joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:06
% DurationCPUTime: 4.38s
% Computational Cost: add. (2444->276), mult. (3472->408), div. (0->0), fcn. (2594->6), ass. (0->160)
t116 = cos(qJ(1));
t109 = pkin(7) + qJ(4);
t100 = cos(t109);
t200 = rSges(6,3) + qJ(5);
t164 = t200 * t100;
t221 = rSges(6,1) + pkin(4);
t99 = sin(t109);
t171 = t221 * t99;
t239 = t171 - t164;
t250 = t116 * t239;
t209 = Icges(6,5) * t99;
t210 = Icges(5,4) * t99;
t249 = t209 - t210 + (Icges(5,1) + Icges(6,1)) * t100;
t178 = qJD(5) * t100;
t242 = t200 * t99;
t212 = t100 * t221 + t242;
t248 = (t212 * qJD(4) - t178) * t116;
t247 = t249 * qJD(4);
t115 = sin(qJ(1));
t246 = t115 / 0.2e1;
t179 = qJD(4) * t116;
t245 = -t179 / 0.2e1;
t244 = -qJD(1) / 0.2e1;
t243 = qJD(1) / 0.2e1;
t110 = t115 ^ 2;
t111 = t116 ^ 2;
t184 = t110 + t111;
t180 = qJD(4) * t115;
t182 = qJD(1) * t116;
t240 = t100 * t180 + t99 * t182;
t107 = t116 * rSges(6,2);
t204 = t115 * t99;
t236 = -t221 * t204 - t107;
t211 = rSges(5,2) * t100;
t158 = rSges(5,1) * t99 + t211;
t127 = t116 * t158;
t137 = Icges(5,2) * t100 + t210;
t47 = Icges(5,6) * t116 + t115 * t137;
t197 = Icges(5,4) * t100;
t141 = Icges(5,1) * t99 + t197;
t51 = Icges(5,5) * t116 + t115 * t141;
t149 = t100 * t47 + t51 * t99;
t124 = t149 * t116;
t131 = -Icges(6,3) * t100 + t209;
t41 = Icges(6,6) * t116 + t115 * t131;
t194 = Icges(6,5) * t100;
t139 = Icges(6,1) * t99 - t194;
t49 = Icges(6,4) * t116 + t115 * t139;
t151 = t100 * t41 - t49 * t99;
t126 = t151 * t116;
t105 = t116 * qJ(2);
t114 = -pkin(6) - qJ(3);
t112 = sin(pkin(7));
t216 = pkin(3) * t112;
t96 = t116 * t216;
t174 = t115 * t114 + t105 + t96;
t218 = -rSges(5,3) - pkin(1);
t33 = t115 * t218 + t127 + t174;
t185 = t116 * pkin(1) + t115 * qJ(2);
t129 = -t114 * t116 + t115 * t216 + t185;
t106 = t116 * rSges(5,3);
t186 = t100 * t115;
t55 = rSges(5,1) * t204 + rSges(5,2) * t186 + t106;
t34 = t129 + t55;
t235 = -t115 * t33 + t116 * t34;
t199 = qJ(2) * t182 + qJD(2) * t115;
t172 = qJD(3) * t116 + t199;
t183 = qJD(1) * t115;
t161 = qJD(1) * t96 + t114 * t183 + t172;
t176 = rSges(5,1) * t240 + t182 * t211;
t187 = qJD(4) * t99;
t13 = (-rSges(5,2) * t187 + qJD(1) * t218) * t115 + t161 + t176;
t103 = qJD(2) * t116;
t162 = -qJD(3) * t115 + t103;
t159 = t114 * t182 + t162;
t166 = -qJ(2) - t216;
t217 = rSges(5,2) * t99;
t78 = rSges(5,1) * t100 - t217;
t14 = t78 * t179 + (t218 * t116 + (-t158 + t166) * t115) * qJD(1) + t159;
t234 = t115 * t14 - t116 * t13;
t233 = -Icges(6,6) * t115 + t116 * t131;
t133 = Icges(5,5) * t99 + Icges(5,6) * t100;
t232 = -Icges(5,3) * t115 + t116 * t133;
t135 = Icges(6,4) * t99 - Icges(6,6) * t100;
t231 = -Icges(6,2) * t115 + t116 * t135;
t230 = -Icges(5,6) * t115 + t116 * t137;
t229 = -Icges(6,4) * t115 + t116 * t139;
t228 = -Icges(5,5) * t115 + t116 * t141;
t165 = qJD(1) * t212;
t215 = -qJD(4) * t239 + qJD(5) * t99;
t15 = t115 * t165 - t116 * t215;
t16 = t115 * t215 + t116 * t165;
t35 = t212 * t115;
t36 = t212 * t116;
t227 = qJD(1) * (-t115 * t36 + t116 * t35) + t16 * t115 - t15 * t116;
t219 = -rSges(6,2) - pkin(1);
t168 = t219 * t115;
t17 = t168 + t174 + t250;
t18 = -t115 * t164 + t129 - t236;
t118 = t115 * t178 - t180 * t242 - t221 * t240;
t2 = (-t116 * t164 + t168) * qJD(1) - t118 + t161;
t3 = t248 + (t219 * t116 + (t166 - t239) * t115) * qJD(1) + t159;
t226 = qJD(1) * (t115 * t18 + t116 * t17) + t115 * t3 - t116 * t2;
t224 = 2 * m(5);
t223 = 2 * m(6);
t222 = m(5) * t78;
t220 = rSges(3,2) - pkin(1);
t214 = t186 * t200 + t236;
t208 = t115 * rSges(6,2);
t213 = t208 - t250;
t207 = t115 * rSges(5,3);
t201 = rSges(4,3) + qJ(3);
t43 = Icges(5,3) * t116 + t115 * t133;
t189 = qJD(1) * t43;
t45 = Icges(6,2) * t116 + t115 * t135;
t188 = qJD(1) * t45;
t181 = qJD(4) * t100;
t177 = -pkin(1) - t201;
t175 = m(6) * t187;
t69 = t158 * qJD(4);
t163 = t184 * t69;
t150 = -t100 * t233 + t229 * t99;
t148 = -t100 * t230 - t228 * t99;
t147 = rSges(4,1) * t112 + rSges(4,2) * cos(pkin(7));
t144 = t115 * t17 - t116 * t18;
t142 = t115 * t35 + t116 * t36;
t136 = -Icges(5,2) * t99 + t197;
t134 = Icges(6,4) * t100 + Icges(6,6) * t99;
t132 = Icges(5,5) * t100 - Icges(5,6) * t99;
t130 = Icges(6,3) * t99 + t194;
t128 = rSges(3,3) * t116 + t115 * t220;
t125 = t150 * t115;
t123 = t148 * t115;
t120 = qJD(4) * t136;
t119 = qJD(4) * t130;
t117 = t115 * t177 + t116 * t147;
t61 = -rSges(3,2) * t116 + t115 * rSges(3,3) + t185;
t60 = t105 + t128;
t57 = t207 - t127;
t40 = t103 + (t220 * t116 + (-rSges(3,3) - qJ(2)) * t115) * qJD(1);
t39 = qJD(1) * t128 + t199;
t38 = t115 * t147 + t116 * t201 + t185;
t37 = t105 + t117;
t26 = t231 * qJD(1) + t134 * t180;
t25 = -t134 * t179 + t188;
t24 = t232 * qJD(1) + t132 * t180;
t23 = -t132 * t179 + t189;
t20 = (t177 * t116 + (-qJ(2) - t147) * t115) * qJD(1) + t162;
t19 = qJD(1) * t117 + t172;
t12 = -t115 * t232 - t116 * t148;
t11 = t115 * t43 - t124;
t10 = -t115 * t231 + t116 * t150;
t9 = t115 * t45 + t126;
t8 = -t116 * t232 + t123;
t7 = t149 * t115 + t116 * t43;
t6 = -t116 * t231 - t125;
t5 = -t151 * t115 + t116 * t45;
t4 = t115 * t214 + t116 * t213;
t1 = ((t208 - t213) * qJD(1) + t118) * t115 + (-t248 + (t115 * t171 + t107 + t214) * qJD(1)) * t116;
t21 = [(t17 * t3 + t18 * t2) * t223 + (t13 * t34 + t14 * t33) * t224 + 0.2e1 * m(4) * (t19 * t38 + t20 * t37) + 0.2e1 * m(3) * (t39 * t61 + t40 * t60) - t247 * t99 + (t137 - t131) * t187 + (-t141 - t139) * t181 + (-t120 + t119) * t100; m(6) * t226 + m(5) * ((t115 * t34 + t116 * t33) * qJD(1) + t234) + m(4) * (t115 * t20 - t116 * t19 + (t115 * t38 + t116 * t37) * qJD(1)) + m(3) * (t115 * t40 - t116 * t39 + (t115 * t61 + t116 * t60) * qJD(1)); 0; m(6) * (-qJD(1) * t144 + t115 * t2 + t116 * t3) + m(5) * (t235 * qJD(1) + t115 * t13 + t116 * t14) + m(4) * (t115 * t19 + t116 * t20 + (-t115 * t37 + t116 * t38) * qJD(1)); 0; 0; m(6) * (t15 * t18 + t16 * t17 - t2 * t36 + t3 * t35) + m(5) * (t234 * t78 + t235 * t69) + ((t230 * t244 - t115 * t120 / 0.2e1 + t233 * t243 + t119 * t246) * t116 + (t47 * t244 + t136 * t179 / 0.2e1 + t41 * t243 + t130 * t245) * t115) * t99 + (t247 * t116 * t246 + t249 * t115 * t245 + ((t228 + t229) * t116 + (t49 + t51) * t115) * t243) * t100 + ((t34 * t222 + (t47 / 0.2e1 - t41 / 0.2e1) * t99 + (-t51 / 0.2e1 - t49 / 0.2e1) * t100) * t115 + (t33 * t222 + (t230 / 0.2e1 - t233 / 0.2e1) * t99 + (-t228 / 0.2e1 - t229 / 0.2e1) * t100) * t116) * qJD(1) + (-t124 / 0.2e1 + t126 / 0.2e1 + (-t148 / 0.2e1 + t150 / 0.2e1) * t115 - t184 * (t133 + t135) / 0.2e1) * qJD(4); -m(5) * t163 + m(6) * t227; m(6) * (-qJD(1) * t142 + t15 * t115 + t116 * t16); ((-t115 * t55 + t116 * t57) * (-t115 * t176 + (t110 * t217 - t111 * t78) * qJD(4) + ((-t55 + t106) * t116 + (t127 - t57 + t207) * t115) * qJD(1)) - t78 * t163) * t224 + t116 * ((t116 * t24 + (t8 + t124) * qJD(1)) * t116 + (-t7 * qJD(1) + (-t181 * t228 + t187 * t230) * t115 + (t23 + (t100 * t51 - t47 * t99) * qJD(4) + (t148 - t43) * qJD(1)) * t116) * t115) + t115 * ((t115 * t23 + (-t11 + t123) * qJD(1)) * t115 + (t12 * qJD(1) + (-t181 * t51 + t187 * t47 + t189) * t116 + (t24 + (t100 * t228 - t230 * t99) * qJD(4) + t149 * qJD(1)) * t115) * t116) + (t1 * t4 - t15 * t36 + t16 * t35) * t223 + t116 * ((t116 * t26 + (t6 - t126) * qJD(1)) * t116 + (-t5 * qJD(1) + (-t181 * t229 - t187 * t233) * t115 + (t25 + (t100 * t49 + t41 * t99) * qJD(4) + (-t150 - t45) * qJD(1)) * t116) * t115) + t115 * ((t115 * t25 + (-t9 - t125) * qJD(1)) * t115 + (t10 * qJD(1) + (-t181 * t49 - t187 * t41 + t188) * t116 + (t26 + (t100 * t229 + t233 * t99) * qJD(4) - t151 * qJD(1)) * t115) * t116) + ((-t5 - t7) * t116 + (-t6 - t8) * t115) * t183 + ((t11 + t9) * t116 + (t10 + t12) * t115) * t182; m(6) * (-t226 * t100 + t144 * t187); t184 * t175; 0; m(6) * ((qJD(4) * t142 + t1) * t99 + (qJD(4) * t4 - t227) * t100); 0.2e1 * (0.1e1 - t184) * t100 * t175;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t21(1), t21(2), t21(4), t21(7), t21(11); t21(2), t21(3), t21(5), t21(8), t21(12); t21(4), t21(5), t21(6), t21(9), t21(13); t21(7), t21(8), t21(9), t21(10), t21(14); t21(11), t21(12), t21(13), t21(14), t21(15);];
Mq = res;
