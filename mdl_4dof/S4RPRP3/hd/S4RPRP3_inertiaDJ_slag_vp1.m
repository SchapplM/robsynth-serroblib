% Calculate time derivative of joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:36
% DurationCPUTime: 2.56s
% Computational Cost: add. (2004->193), mult. (2482->261), div. (0->0), fcn. (1876->6), ass. (0->119)
t205 = Icges(4,5) + Icges(5,5);
t204 = -Icges(4,6) - Icges(5,6);
t203 = Icges(4,3) + Icges(5,3);
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t202 = t204 * t86 + t205 * t88;
t84 = qJ(1) + pkin(6);
t81 = sin(t84);
t82 = cos(t84);
t200 = t202 * t82 + t203 * t81;
t201 = t202 * t81 - t203 * t82;
t146 = Icges(5,4) * t88;
t111 = -Icges(5,2) * t86 + t146;
t39 = -Icges(5,6) * t82 + t111 * t81;
t148 = Icges(4,4) * t88;
t113 = -Icges(4,2) * t86 + t148;
t41 = -Icges(4,6) * t82 + t113 * t81;
t199 = t41 + t39;
t40 = Icges(5,6) * t81 + t111 * t82;
t42 = Icges(4,6) * t81 + t113 * t82;
t198 = t42 + t40;
t147 = Icges(5,4) * t86;
t115 = Icges(5,1) * t88 - t147;
t43 = -Icges(5,5) * t82 + t115 * t81;
t149 = Icges(4,4) * t86;
t117 = Icges(4,1) * t88 - t149;
t45 = -Icges(4,5) * t82 + t117 * t81;
t197 = t45 + t43;
t44 = Icges(5,5) * t81 + t115 * t82;
t46 = Icges(4,5) * t81 + t117 * t82;
t196 = t46 + t44;
t195 = (t204 * t88 - t205 * t86) * qJD(3);
t118 = t42 * t86 - t46 * t88;
t120 = t40 * t86 - t44 * t88;
t174 = t118 + t120;
t194 = t174 * t82;
t119 = t41 * t86 - t45 * t88;
t121 = t39 * t86 - t43 * t88;
t193 = t119 + t121;
t160 = rSges(5,1) * t88;
t126 = -rSges(5,2) * t86 + t160;
t80 = pkin(3) * t88 + pkin(2);
t106 = -t126 - t80;
t191 = t200 * qJD(1);
t190 = t196 * t81 + t197 * t82;
t189 = t198 * t81 + t199 * t82;
t169 = t81 ^ 2;
t168 = t82 ^ 2;
t188 = t193 * t81 + t201 * t82;
t101 = t118 * t81;
t103 = t120 * t81;
t187 = -t200 * t82 - t101 - t103;
t102 = t119 * t82;
t104 = t121 * t82;
t186 = -t201 * t81 + t102 + t104;
t185 = t200 * t81 - t194;
t135 = qJD(1) * t86;
t85 = -qJ(4) - pkin(5);
t162 = -pkin(5) - t85;
t155 = t82 * t88;
t156 = t82 * t86;
t176 = rSges(5,1) * t155 - rSges(5,2) * t156 + t81 * rSges(5,3) + t82 * t80;
t151 = -t82 * pkin(2) + t162 * t81 + t176;
t157 = t82 * t85;
t79 = t82 * pkin(5);
t152 = -rSges(5,3) * t82 + t157 + t79 + (-pkin(2) - t106) * t81;
t153 = rSges(5,3) - t85;
t130 = t81 * t135;
t136 = qJD(1) * t82;
t154 = t88 * rSges(5,2);
t172 = (rSges(5,1) + pkin(3)) * t86 + t154;
t92 = qJD(3) * t172;
t90 = rSges(5,2) * t130 + rSges(5,3) * t136 + qJD(4) * t81 - t82 * t92;
t1 = ((t162 * t82 + t152) * qJD(1) + t90) * t82 + (-t81 * t92 + ((-pkin(5) + t153) * t81 - t151) * qJD(1) + (-rSges(5,2) * t135 - qJD(4)) * t82) * t81;
t170 = 2 * m(5);
t184 = t1 * t170;
t171 = 2 * m(4);
t71 = rSges(4,1) * t86 + rSges(4,2) * t88;
t100 = qJD(3) * t71;
t131 = rSges(4,2) * t156;
t159 = rSges(4,2) * t86;
t161 = rSges(4,1) * t88;
t127 = -t159 + t161;
t158 = rSges(4,3) * t82;
t48 = t127 * t81 - t158;
t78 = t81 * rSges(4,3);
t150 = rSges(4,1) * t155 + t78;
t50 = -t131 + t150;
t91 = rSges(4,2) * t130 + rSges(4,3) * t136 - t100 * t82;
t2 = (qJD(1) * t48 + t91) * t82 + (-t81 * t100 + (-t131 - t50 + t78) * qJD(1)) * t81;
t183 = t171 * t2;
t182 = -t201 * qJD(1) + t195 * t82;
t181 = -t195 * t81 - t191;
t164 = sin(qJ(1)) * pkin(1);
t175 = t79 - t164;
t137 = qJD(1) * t81;
t165 = m(4) * t71;
t83 = cos(qJ(1)) * pkin(1);
t134 = qJD(3) * t81;
t133 = qJD(3) * t86;
t132 = qJD(3) * t88;
t70 = rSges(5,1) * t86 + t154;
t129 = -pkin(3) * t86 - t70;
t107 = -pkin(2) - t127;
t60 = t127 * qJD(3);
t59 = t126 * qJD(3);
t52 = t129 * t82;
t51 = t129 * t81;
t32 = pkin(5) * t81 + t83 + (pkin(2) - t159) * t82 + t150;
t31 = t107 * t81 + t158 + t175;
t30 = -t81 * t85 + t176 + t83;
t29 = t106 * t81 + t153 * t82 - t164;
t16 = -t70 * t136 - t59 * t81 + (-t132 * t81 - t135 * t82) * pkin(3);
t15 = t70 * t137 - t59 * t82 + (-t132 * t82 + t130) * pkin(3);
t14 = t71 * t134 + (-t83 + (-rSges(4,3) - pkin(5)) * t81 + t107 * t82) * qJD(1);
t13 = ((-pkin(2) - t161) * t81 + t175) * qJD(1) + t91;
t12 = qJD(4) * t82 + t172 * t134 + (t106 * t82 - t153 * t81 - t83) * qJD(1);
t11 = (-t164 - t157 + (-t80 - t160) * t81) * qJD(1) + t90;
t3 = [(t13 * t32 + t14 * t31) * t171 + (t11 * t30 + t12 * t29) * t170 + (t115 + t117 - t147 - t149 + (-Icges(4,2) - Icges(5,2)) * t88) * t133 + (t111 + t113 + t146 + t148 + (Icges(4,1) + Icges(5,1)) * t86) * t132; 0; 0; m(4) * ((-t13 * t81 - t14 * t82) * t71 + (-t31 * t82 - t32 * t81) * t60) + m(5) * (t11 * t51 + t12 * t52 + t15 * t29 + t16 * t30) + (-t101 / 0.2e1 - t103 / 0.2e1 + t102 / 0.2e1 + t104 / 0.2e1 + t202 * (t169 / 0.2e1 + t168 / 0.2e1)) * qJD(3) + ((-t32 * t165 + (t42 / 0.2e1 + t40 / 0.2e1) * t88 + (t46 / 0.2e1 + t44 / 0.2e1) * t86) * t82 + (t31 * t165 + (t41 / 0.2e1 + t39 / 0.2e1) * t88 + (t45 / 0.2e1 + t43 / 0.2e1) * t86) * t81 + (-t197 * t86 - t199 * t88) * t81 / 0.2e1 - (t196 * t86 + t198 * t88) * t82 / 0.2e1) * qJD(1); m(4) * t2 + m(5) * t1; (t168 + t169) * t71 * t60 * t171 + (t52 * t15 + t51 * t16) * t170 + (t151 * t184 + t50 * t183 + t181 * t168 + (-t193 * t82 - t187) * t136) * t82 + (t48 * t183 + t152 * t184 + t182 * t169 + (t174 * t81 - t186) * t137 + ((t181 - t191) * t81 + t182 * t82 + t190 * t133 + t189 * t132 + (-t189 * t88 - t190 * t86) * qJD(3) + ((-t193 + t200) * t81 + t194 + t185 + t188) * qJD(1)) * t82) * t81 + (t187 * t81 + t188 * t82) * t137 + (t185 * t81 + t186 * t82) * t136; m(5) * (-t11 * t82 + t12 * t81 + (t29 * t82 + t30 * t81) * qJD(1)); 0; m(5) * (t15 * t81 - t16 * t82 + (t51 * t81 + t52 * t82) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
