% Calculate time derivative of joint inertia matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:48
% DurationCPUTime: 2.46s
% Computational Cost: add. (1972->186), mult. (2426->258), div. (0->0), fcn. (1842->4), ass. (0->116)
t198 = Icges(4,5) + Icges(5,5);
t197 = -Icges(4,6) - Icges(5,6);
t196 = Icges(4,3) + Icges(5,3);
t85 = sin(qJ(3));
t86 = cos(qJ(3));
t195 = t197 * t85 + t198 * t86;
t83 = pkin(6) + qJ(2);
t81 = sin(t83);
t82 = cos(t83);
t193 = t195 * t82 + t196 * t81;
t194 = t195 * t81 - t196 * t82;
t141 = Icges(5,4) * t86;
t106 = -Icges(5,2) * t85 + t141;
t39 = -Icges(5,6) * t82 + t106 * t81;
t143 = Icges(4,4) * t86;
t108 = -Icges(4,2) * t85 + t143;
t41 = -Icges(4,6) * t82 + t108 * t81;
t192 = t41 + t39;
t40 = Icges(5,6) * t81 + t106 * t82;
t42 = Icges(4,6) * t81 + t108 * t82;
t191 = t42 + t40;
t142 = Icges(5,4) * t85;
t110 = Icges(5,1) * t86 - t142;
t43 = -Icges(5,5) * t82 + t110 * t81;
t144 = Icges(4,4) * t85;
t112 = Icges(4,1) * t86 - t144;
t45 = -Icges(4,5) * t82 + t112 * t81;
t190 = t43 + t45;
t44 = Icges(5,5) * t81 + t110 * t82;
t46 = Icges(4,5) * t81 + t112 * t82;
t189 = t46 + t44;
t188 = (t197 * t86 - t198 * t85) * qJD(3);
t84 = -qJ(4) - pkin(5);
t148 = rSges(5,3) - t84;
t187 = t148 * t82;
t113 = t42 * t85 - t46 * t86;
t115 = t40 * t85 - t44 * t86;
t167 = t113 + t115;
t186 = t167 * t82;
t114 = t41 * t85 - t45 * t86;
t116 = t39 * t85 - t43 * t86;
t185 = t114 + t116;
t154 = rSges(5,1) * t86;
t121 = -rSges(5,2) * t85 + t154;
t80 = pkin(3) * t86 + pkin(2);
t101 = -t121 - t80;
t183 = t193 * qJD(2);
t182 = t189 * t81 + t190 * t82;
t181 = t191 * t81 + t192 * t82;
t162 = t81 ^ 2;
t161 = t82 ^ 2;
t180 = t185 * t81 + t194 * t82;
t96 = t113 * t81;
t98 = t115 * t81;
t179 = -t193 * t82 - t96 - t98;
t97 = t114 * t82;
t99 = t116 * t82;
t178 = -t194 * t81 + t97 + t99;
t177 = t193 * t81 - t186;
t130 = qJD(2) * t85;
t124 = t81 * t130;
t131 = qJD(2) * t82;
t125 = rSges(5,2) * t124 + rSges(5,3) * t131 + qJD(4) * t81;
t156 = -pkin(5) - t84;
t150 = t82 * t86;
t151 = t82 * t85;
t168 = rSges(5,1) * t150 - rSges(5,2) * t151 + t81 * rSges(5,3) + t82 * t80;
t146 = -t82 * pkin(2) + t156 * t81 + t168;
t79 = t82 * pkin(5);
t147 = t79 - t187 + (-pkin(2) - t101) * t81;
t149 = t86 * rSges(5,2);
t165 = (rSges(5,1) + pkin(3)) * t85 + t149;
t88 = t165 * qJD(3);
t1 = (-t82 * t88 + (t156 * t82 + t147) * qJD(2) + t125) * t82 + (-t81 * t88 + ((-pkin(5) + t148) * t81 - t146) * qJD(2) + (-rSges(5,2) * t130 - qJD(4)) * t82) * t81;
t163 = 2 * m(5);
t176 = t1 * t163;
t164 = 2 * m(4);
t126 = rSges(4,2) * t151;
t153 = rSges(4,2) * t85;
t155 = rSges(4,1) * t86;
t122 = -t153 + t155;
t152 = t82 * rSges(4,3);
t48 = t122 * t81 - t152;
t78 = t81 * rSges(4,3);
t145 = rSges(4,1) * t150 + t78;
t50 = -t126 + t145;
t71 = t85 * rSges(4,1) + rSges(4,2) * t86;
t95 = t71 * qJD(3);
t87 = rSges(4,2) * t124 + rSges(4,3) * t131 - t82 * t95;
t2 = (qJD(2) * t48 + t87) * t82 + (-t81 * t95 + (-t126 - t50 + t78) * qJD(2)) * t81;
t175 = t164 * t2;
t174 = -t194 * qJD(2) + t188 * t82;
t173 = -t188 * t81 - t183;
t132 = qJD(2) * t81;
t158 = m(4) * t71;
t129 = qJD(3) * t81;
t128 = qJD(3) * t85;
t127 = qJD(3) * t86;
t70 = t85 * rSges(5,1) + t149;
t123 = -pkin(3) * t85 - t70;
t102 = -pkin(2) - t122;
t60 = t122 * qJD(3);
t59 = t121 * qJD(3);
t52 = t123 * t82;
t51 = t123 * t81;
t32 = pkin(5) * t81 + (pkin(2) - t153) * t82 + t145;
t31 = t102 * t81 + t152 + t79;
t30 = -t81 * t84 + t168;
t29 = t101 * t81 + t187;
t16 = -t70 * t131 - t81 * t59 + (-t127 * t81 - t130 * t82) * pkin(3);
t15 = t70 * t132 - t82 * t59 + (-t127 * t82 + t124) * pkin(3);
t14 = t71 * t129 + ((-rSges(4,3) - pkin(5)) * t81 + t102 * t82) * qJD(2);
t13 = (t79 + (-pkin(2) - t155) * t81) * qJD(2) + t87;
t12 = qJD(4) * t82 + t165 * t129 + (t101 * t82 - t148 * t81) * qJD(2);
t11 = (-t80 - t154) * t132 + (-qJD(2) * t84 - t88) * t82 + t125;
t3 = [0; 0; (t13 * t32 + t14 * t31) * t164 + (t11 * t30 + t12 * t29) * t163 + (t110 + t112 - t142 - t144 + (-Icges(4,2) - Icges(5,2)) * t86) * t128 + (t106 + t108 + t141 + t143 + (Icges(4,1) + Icges(5,1)) * t85) * t127; m(4) * t2 + m(5) * t1; m(4) * ((-t13 * t81 - t14 * t82) * t71 + (-t31 * t82 - t32 * t81) * t60) + m(5) * (t11 * t51 + t12 * t52 + t15 * t29 + t16 * t30) + (-t96 / 0.2e1 - t98 / 0.2e1 + t97 / 0.2e1 + t99 / 0.2e1 + t195 * (t162 / 0.2e1 + t161 / 0.2e1)) * qJD(3) + ((-t32 * t158 + (t42 / 0.2e1 + t40 / 0.2e1) * t86 + (t46 / 0.2e1 + t44 / 0.2e1) * t85) * t82 + (t31 * t158 + (t41 / 0.2e1 + t39 / 0.2e1) * t86 + (t45 / 0.2e1 + t43 / 0.2e1) * t85) * t81 + (-t190 * t85 - t192 * t86) * t81 / 0.2e1 - (t189 * t85 + t191 * t86) * t82 / 0.2e1) * qJD(2); (t161 + t162) * t71 * t60 * t164 + (t52 * t15 + t51 * t16) * t163 + (t146 * t176 + t50 * t175 + t173 * t161 + (-t185 * t82 - t179) * t131) * t82 + (t48 * t175 + t147 * t176 + t174 * t162 + (t167 * t81 - t178) * t132 + ((t173 - t183) * t81 + t174 * t82 + t182 * t128 + t181 * t127 + (-t181 * t86 - t182 * t85) * qJD(3) + ((-t185 + t193) * t81 + t186 + t177 + t180) * qJD(2)) * t82) * t81 + (t179 * t81 + t180 * t82) * t132 + (t177 * t81 + t178 * t82) * t131; 0; m(5) * (-t11 * t82 + t12 * t81 + (t29 * t82 + t30 * t81) * qJD(2)); m(5) * (t15 * t81 - t16 * t82 + (t51 * t81 + t52 * t82) * qJD(2)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
