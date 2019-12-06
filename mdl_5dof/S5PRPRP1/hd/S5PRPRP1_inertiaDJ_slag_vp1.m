% Calculate time derivative of joint inertia matrix for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:26
% DurationCPUTime: 2.67s
% Computational Cost: add. (3485->204), mult. (3080->286), div. (0->0), fcn. (2329->6), ass. (0->128)
t227 = Icges(6,4) + Icges(5,5);
t226 = -Icges(5,6) + Icges(6,6);
t225 = Icges(6,2) + Icges(5,3);
t100 = pkin(8) + qJ(4);
t96 = sin(t100);
t98 = cos(t100);
t223 = t226 * t96 + t227 * t98;
t184 = rSges(6,1) + pkin(4);
t224 = t184 * t98;
t101 = pkin(7) + qJ(2);
t97 = sin(t101);
t99 = cos(t101);
t221 = t223 * t99 + t225 * t97;
t222 = t223 * t97 - t225 * t99;
t166 = Icges(6,5) * t98;
t124 = Icges(6,3) * t96 + t166;
t40 = -Icges(6,6) * t99 + t124 * t97;
t170 = Icges(5,4) * t98;
t128 = -Icges(5,2) * t96 + t170;
t46 = -Icges(5,6) * t99 + t128 * t97;
t220 = t46 - t40;
t41 = Icges(6,6) * t97 + t124 * t99;
t47 = Icges(5,6) * t97 + t128 * t99;
t219 = t47 - t41;
t167 = Icges(6,5) * t96;
t130 = Icges(6,1) * t98 + t167;
t48 = -Icges(6,4) * t99 + t130 * t97;
t171 = Icges(5,4) * t96;
t132 = Icges(5,1) * t98 - t171;
t50 = -Icges(5,5) * t99 + t132 * t97;
t218 = t50 + t48;
t49 = Icges(6,4) * t97 + t130 * t99;
t51 = Icges(5,5) * t97 + t132 * t99;
t217 = t51 + t49;
t133 = t47 * t96 - t51 * t98;
t135 = t41 * t96 + t49 * t98;
t192 = t133 - t135;
t216 = t192 * t99;
t152 = qJD(4) * t98;
t158 = rSges(6,3) + qJ(5);
t215 = t158 * t152;
t214 = (t226 * t98 - t227 * t96) * qJD(4);
t134 = t46 * t96 - t50 * t98;
t136 = t40 * t96 + t48 * t98;
t213 = t134 - t136;
t174 = -t158 * t98 + t184 * t96;
t161 = qJ(5) * t96;
t210 = rSges(6,3) * t96 + t161 + t224;
t209 = t221 * qJD(2);
t116 = t133 * t97;
t117 = t134 * t99;
t118 = t135 * t97;
t119 = t136 * t99;
t208 = t221 * t99 - t222 * t97 + t116 + t117 - t118 - t119;
t207 = t217 * t97 + t218 * t99;
t206 = t219 * t97 + t220 * t99;
t94 = t97 ^ 2;
t95 = t99 ^ 2;
t205 = t213 * t97 + t222 * t99;
t189 = 2 * m(5);
t75 = rSges(5,1) * t96 + rSges(5,2) * t98;
t115 = t75 * qJD(4);
t154 = qJD(2) * t99;
t155 = qJD(2) * t97;
t182 = rSges(5,2) * t96;
t106 = rSges(5,3) * t154 - t115 * t99 + t155 * t182;
t179 = t96 * t99;
t91 = t97 * rSges(5,3);
t193 = -rSges(5,2) * t179 + t91;
t194 = t97 * t115;
t183 = rSges(5,1) * t98;
t145 = -t182 + t183;
t53 = -t99 * rSges(5,3) + t145 * t97;
t178 = t98 * t99;
t55 = rSges(5,1) * t178 + t193;
t2 = (qJD(2) * t53 + t106) * t99 + (-t194 + (-t55 + t193) * qJD(2)) * t97;
t203 = t189 * t2;
t201 = t221 * t97 - t216;
t200 = -qJD(2) * t222 + t214 * t99;
t199 = -t214 * t97 - t209;
t104 = -pkin(6) - qJ(3);
t103 = cos(pkin(8));
t93 = pkin(3) * t103 + pkin(2);
t122 = -t145 - t93;
t34 = (rSges(5,3) - t104) * t99 + t122 * t97;
t146 = -t104 * t97 + t93 * t99;
t35 = t146 + t55;
t190 = t34 * t99 + t35 * t97;
t121 = rSges(4,1) * t103 - rSges(4,2) * sin(pkin(8)) + pkin(2);
t159 = rSges(4,3) + qJ(3);
t37 = t121 * t99 + t159 * t97;
t188 = 2 * m(6);
t185 = m(5) * t75;
t92 = t97 * rSges(6,2);
t177 = -t99 * rSges(6,2) + t210 * t97;
t176 = rSges(6,3) * t179 + t99 * t161 + t178 * t184 + t92;
t175 = -qJD(4) * t210 + qJD(5) * t98;
t90 = qJD(3) * t99;
t173 = t104 * t155 + t90;
t172 = t94 + t95;
t160 = t104 * t99;
t153 = qJD(4) * t96;
t151 = qJD(5) * t96;
t149 = m(6) * t153;
t39 = t174 * t99;
t147 = t184 * qJD(4);
t108 = -t158 * t96 - t224 - t93;
t105 = t108 * t97;
t18 = (rSges(6,2) - t104) * t99 + t105;
t19 = t146 + t176;
t142 = t18 * t99 + t19 * t97;
t38 = t174 * t97;
t137 = -t38 * t97 - t39 * t99;
t107 = rSges(6,2) * t154 - t147 * t179 + (t151 + t215) * t99;
t36 = -t121 * t97 + t159 * t99;
t89 = qJD(3) * t97;
t66 = t145 * qJD(4);
t33 = -qJD(2) * t37 + t90;
t32 = qJD(2) * t36 + t89;
t17 = t194 + (t122 * t99 - t91) * qJD(2) + t173;
t16 = t89 + (-t160 + (-t93 - t183) * t97) * qJD(2) + t106;
t15 = -qJD(2) * t39 + t175 * t97;
t14 = t155 * t174 + t175 * t99;
t5 = t176 * t99 + t177 * t97;
t4 = (qJD(4) * t174 - t151) * t97 + (t108 * t99 - t92) * qJD(2) + t173;
t3 = t89 + (t105 - t160) * qJD(2) + t107;
t1 = (qJD(2) * t177 + t107) * t99 + ((-t176 + t92) * qJD(2) + (t215 + (qJD(5) - t147) * t96) * t97) * t97;
t6 = [0; 0; 0.2e1 * m(4) * (t32 * t37 + t33 * t36) + (t16 * t35 + t17 * t34) * t189 + (t18 * t4 + t19 * t3) * t188 + (t130 + t132 + t167 - t171 + (-Icges(5,2) - Icges(6,3)) * t98) * t153 + (-t124 + t128 - t166 + t170 + (Icges(5,1) + Icges(6,1)) * t96) * t152; 0; m(4) * (-t32 * t99 + t33 * t97 + (t36 * t99 + t37 * t97) * qJD(2)) + m(5) * (qJD(2) * t190 - t16 * t99 + t17 * t97) + m(6) * (qJD(2) * t142 - t3 * t99 + t4 * t97); 0; m(5) * t2 + m(6) * t1; m(5) * ((-t16 * t97 - t17 * t99) * t75 - t190 * t66) + m(6) * (t14 * t18 + t15 * t19 - t3 * t38 - t39 * t4) + (-t116 / 0.2e1 + t118 / 0.2e1 + t117 / 0.2e1 - t119 / 0.2e1 + t223 * (t94 / 0.2e1 + t95 / 0.2e1)) * qJD(4) + ((-t35 * t185 + (t47 / 0.2e1 - t41 / 0.2e1) * t98 + (t51 / 0.2e1 + t49 / 0.2e1) * t96) * t99 + (t34 * t185 + (t46 / 0.2e1 - t40 / 0.2e1) * t98 + (t50 / 0.2e1 + t48 / 0.2e1) * t96) * t97 + (-t218 * t96 - t220 * t98) * t97 / 0.2e1 - (t217 * t96 + t219 * t98) * t99 / 0.2e1) * qJD(2); m(6) * (qJD(2) * t137 + t14 * t97 - t15 * t99); t172 * t75 * t66 * t189 + (t1 * t5 - t14 * t39 - t15 * t38) * t188 + (t205 * t155 + t199 * t95 + t55 * t203 + (-t213 * t99 + t208) * t154) * t99 + (t53 * t203 + t200 * t94 + t201 * t154 + ((t199 - t209) * t97 + t200 * t99 + t207 * t153 + t206 * t152 + (-t206 * t98 - t207 * t96) * qJD(4) + ((-t213 + t221) * t97 + t216 + t201 + t205) * qJD(2)) * t99 + (t192 * t97 - t208) * t155) * t97; t149; m(6) * (t142 * t152 + (t3 * t97 + t4 * t99 + (-t18 * t97 + t19 * t99) * qJD(2)) * t96); 0; m(6) * ((qJD(4) * t137 - t1) * t98 + (qJD(4) * t5 + t14 * t99 + t15 * t97 + (-t38 * t99 + t39 * t97) * qJD(2)) * t96); 0.2e1 * (-0.1e1 + t172) * t98 * t149;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
