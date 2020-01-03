% Calculate time derivative of joint inertia matrix for
% S5RPPRP5
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:28
% DurationCPUTime: 3.55s
% Computational Cost: add. (2483->270), mult. (6910->390), div. (0->0), fcn. (6687->6), ass. (0->128)
t226 = Icges(5,1) + Icges(6,1);
t225 = Icges(5,4) - Icges(6,5);
t222 = Icges(6,4) + Icges(5,5);
t224 = Icges(5,2) + Icges(6,3);
t221 = Icges(6,2) + Icges(5,3);
t229 = Icges(5,6) - Icges(6,6);
t125 = cos(qJ(1));
t122 = sin(pkin(7));
t123 = cos(pkin(7));
t180 = sin(qJ(4));
t181 = cos(qJ(4));
t101 = t122 * t181 - t123 * t180;
t124 = sin(qJ(1));
t94 = t101 * t124;
t100 = t122 * t180 + t123 * t181;
t95 = t100 * t124;
t202 = -t229 * t125 - t224 * t94 - t225 * t95;
t96 = t101 * t125;
t97 = t100 * t125;
t201 = t229 * t124 - t224 * t96 - t225 * t97;
t228 = -t221 * t125 - t222 * t95 - t229 * t94;
t227 = -t221 * t124 + t222 * t97 + t229 * t96;
t200 = -t222 * t125 - t225 * t94 - t226 * t95;
t199 = t222 * t124 - t225 * t96 - t226 * t97;
t208 = rSges(6,3) + qJ(5);
t218 = rSges(6,1) + pkin(4);
t223 = t208 * t94 - t218 * t95;
t219 = t227 * t125 - t199 * t95 - t201 * t94;
t217 = t228 * t124 - t200 * t97 - t202 * t96;
t156 = qJD(1) * t125;
t150 = qJD(4) * t180;
t151 = qJD(4) * t181;
t160 = t123 * t125;
t161 = t122 * t125;
t194 = t101 * qJD(1);
t62 = -t124 * t194 - t150 * t161 - t151 * t160;
t131 = qJD(4) * t101;
t157 = qJD(1) * t124;
t63 = -t100 * t157 + t125 * t131;
t216 = -t156 * t229 + t224 * t62 + t225 * t63;
t64 = qJD(4) * t95 - t125 * t194;
t65 = qJD(1) * t97 + t124 * t131;
t215 = -t157 * t229 - t224 * t64 + t225 * t65;
t214 = t222 * t156 - t225 * t62 - t226 * t63;
t213 = t222 * t157 + t225 * t64 - t226 * t65;
t98 = t100 * qJD(4);
t99 = t122 * t151 - t123 * t150;
t212 = -t222 * t98 - t229 * t99;
t211 = -t224 * t99 - t225 * t98;
t210 = -t225 * t99 - t226 * t98;
t209 = -t224 * t100 + t225 * t101;
t197 = -t225 * t100 + t226 * t101;
t207 = t228 * t125 + t200 * t95 + t202 * t94;
t191 = 2 * m(5);
t145 = t65 * rSges(5,1) - t64 * rSges(5,2);
t176 = t63 * rSges(5,1) + t62 * rSges(5,2);
t195 = t124 ^ 2 + t125 ^ 2;
t144 = -t95 * rSges(5,1) - t94 * rSges(5,2);
t52 = rSges(5,3) * t125 - t144;
t173 = t97 * rSges(5,1) + t96 * rSges(5,2);
t54 = -rSges(5,3) * t124 + t173;
t2 = -t124 * t145 - t125 * t176 + (rSges(5,3) * t195 + t124 * t54 - t125 * t52) * qJD(1);
t206 = t191 * t2;
t205 = -t227 * t124 - t199 * t97 - t201 * t96;
t204 = t221 * t156 - t222 * t63 - t229 * t62;
t203 = -t221 * t157 + t222 * t65 - t229 * t64;
t196 = t208 * t96 - t218 * t97;
t118 = t125 * qJ(2);
t136 = -qJ(3) * t122 - pkin(1) + (-pkin(2) - pkin(3)) * t123;
t134 = t136 * t124;
t182 = -rSges(5,3) - pkin(6);
t127 = t125 * t182 + t134;
t32 = t118 + t127 + t144;
t158 = t125 * pkin(1) + t124 * qJ(2);
t148 = pkin(2) * t160 + qJ(3) * t161 + t158;
t139 = pkin(3) * t160 + t148;
t33 = t124 * t182 + t139 + t173;
t193 = t124 * t33 + t125 * t32;
t155 = qJD(3) * t122;
t159 = qJ(2) * t156 + qJD(2) * t124;
t154 = t125 * t155 + t159;
t16 = qJD(1) * t127 + t154 + t176;
t133 = t136 * t125;
t116 = qJD(2) * t125;
t143 = -t124 * t155 + t116;
t137 = pkin(6) * t157 + t143;
t17 = ((rSges(5,3) - qJ(2)) * t124 + t133) * qJD(1) + t137 - t145;
t192 = t124 * t16 + t125 * t17;
t190 = 2 * m(6);
t189 = m(5) / 0.2e1;
t188 = m(6) / 0.2e1;
t183 = -rSges(6,2) - pkin(6);
t179 = qJD(5) * t100 + t208 * t99 - t218 * t98;
t171 = rSges(6,2) * t125;
t178 = -t171 + t223;
t172 = rSges(6,2) * t124;
t177 = t172 + t196;
t175 = t208 * t100 + t218 * t101;
t35 = t175 * t125;
t77 = -rSges(5,1) * t98 - rSges(5,2) * t99;
t149 = t195 * t77;
t142 = rSges(3,1) * t123 - rSges(3,2) * t122;
t138 = -pkin(1) - t142;
t135 = t96 * qJD(5) + t208 * t62 - t218 * t63;
t132 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t123 + (-rSges(4,3) - qJ(3)) * t122;
t130 = t94 * qJD(5) - t208 * t64 - t218 * t65;
t129 = rSges(3,3) * t125 + t124 * t138;
t128 = t125 * t183 + t134;
t126 = rSges(4,2) * t125 + t124 * t132;
t88 = t124 * rSges(3,3) + t125 * t142 + t158;
t87 = t118 + t129;
t86 = rSges(5,1) * t101 - rSges(5,2) * t100;
t69 = t116 + ((-rSges(3,3) - qJ(2)) * t124 + t138 * t125) * qJD(1);
t68 = qJD(1) * t129 + t159;
t61 = t124 * rSges(4,2) + (rSges(4,1) * t123 + rSges(4,3) * t122) * t125 + t148;
t60 = t118 + t126;
t37 = ((-rSges(4,2) - qJ(2)) * t124 + t132 * t125) * qJD(1) + t143;
t36 = qJD(1) * t126 + t154;
t34 = t175 * t124;
t19 = t124 * t183 + t139 - t196;
t18 = t118 + t128 + t223;
t15 = qJD(1) * t35 + t124 * t179;
t14 = t125 * t179 - t157 * t175;
t5 = t124 * t178 + t125 * t177;
t4 = ((rSges(6,2) - qJ(2)) * t124 + t133) * qJD(1) + t130 + t137;
t3 = qJD(1) * t128 - t135 + t154;
t1 = t135 * t125 + t130 * t124 + ((t171 + t178) * t125 + (t172 - t177) * t124) * qJD(1);
t6 = [(t18 * t4 + t19 * t3) * t190 + (t16 * t33 + t17 * t32) * t191 + 0.2e1 * m(3) * (t68 * t88 + t69 * t87) + 0.2e1 * m(4) * (t36 * t61 + t37 * t60) - t209 * t99 - t197 * t98 + t210 * t101 - t211 * t100; m(6) * (t124 * t4 - t125 * t3 + (t124 * t19 + t125 * t18) * qJD(1)) + m(5) * (qJD(1) * t193 + t124 * t17 - t125 * t16) + m(3) * (t124 * t69 - t125 * t68 + (t124 * t88 + t125 * t87) * qJD(1)) + m(4) * (t124 * t37 - t125 * t36 + (t124 * t61 + t125 * t60) * qJD(1)); 0; 0.2e1 * ((t124 * t3 + t125 * t4 + t156 * t19 - t157 * t18) * t188 + (t156 * t33 - t157 * t32 + t192) * t189 + m(4) * (t124 * t36 + t125 * t37 + t156 * t61 - t157 * t60) / 0.2e1) * t122; 0; 0; m(6) * (t14 * t18 + t15 * t19 + t3 * t34 + t35 * t4) + m(5) * (t192 * t86 + t193 * t77) + (m(5) * (-t124 * t32 + t125 * t33) * t86 - (t201 * t100 - t199 * t101 + t197 * t97 + t209 * t96) * t125 / 0.2e1) * qJD(1) - ((t197 * t95 + t209 * t94) * qJD(1) + t201 * t99 + t199 * t98 + t210 * t97 + t211 * t96 + t197 * t63 + t209 * t62 - t212 * t124 + (-qJD(1) * t200 - t214) * t101 + (qJD(1) * t202 - t216) * t100) * t124 / 0.2e1 + (-t215 * t100 - t213 * t101 + t212 * t125 + t197 * t65 + t200 * t98 + t202 * t99 - t209 * t64 + t210 * t95 + t211 * t94) * t125 / 0.2e1; m(6) * (t14 * t124 - t125 * t15 + (t124 * t34 + t125 * t35) * qJD(1)); 0.2e1 * (-m(5) * t2 / 0.2e1 - m(6) * t1 / 0.2e1) * t123 + 0.2e1 * ((t124 * t15 + t125 * t14 + (-t124 * t35 + t125 * t34) * qJD(1)) * t188 + t149 * t189) * t122; t86 * t149 * t191 + (t1 * t5 + t14 * t35 + t15 * t34) * t190 + (-t54 * t206 + (-t219 * qJD(1) + t125 * t203 - t200 * t65 + t202 * t64 - t213 * t95 + t215 * t94) * t125 + t207 * t157 - t217 * t156) * t125 + (-t52 * t206 + (t217 * qJD(1) + t124 * t204 - t199 * t63 - t201 * t62 - t214 * t97 + t216 * t96) * t124 + t219 * t157 + t205 * t156 + (t203 * t124 + t204 * t125 + t213 * t97 - t215 * t96 + t214 * t95 - t216 * t94 + t199 * t65 - t201 * t64 + t200 * t63 + t202 * t62 + (t205 + t207) * qJD(1)) * t125) * t124; m(6) * (-t18 * t62 + t19 * t64 - t3 * t94 - t4 * t96); m(6) * (-t62 * t124 - t64 * t125 + (-t124 * t94 - t125 * t96) * qJD(1)); m(6) * (-t99 * t123 + (t124 * t64 - t125 * t62 + (t96 * t124 - t94 * t125) * qJD(1)) * t122); m(6) * (t1 * t100 - t14 * t96 - t15 * t94 + t34 * t64 - t35 * t62 + t5 * t99); (t100 * t99 + t62 * t96 - t64 * t94) * t190;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
