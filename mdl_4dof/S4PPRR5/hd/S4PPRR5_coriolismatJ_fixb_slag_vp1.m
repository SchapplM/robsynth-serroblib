% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:42
% EndTime: 2019-12-31 16:19:46
% DurationCPUTime: 2.53s
% Computational Cost: add. (3783->272), mult. (10011->449), div. (0->0), fcn. (10925->6), ass. (0->163)
t207 = Icges(4,1) - Icges(4,2);
t206 = 2 * Icges(4,4);
t126 = sin(qJ(3));
t205 = t207 * t126;
t128 = cos(qJ(3));
t123 = sin(pkin(6));
t121 = t123 ^ 2;
t124 = cos(pkin(6));
t122 = t124 ^ 2;
t151 = t121 + t122;
t204 = t151 * (-t126 * rSges(4,1) - rSges(4,2) * t128);
t202 = t126 * t206 - t207 * t128;
t127 = cos(qJ(4));
t125 = sin(qJ(4));
t163 = Icges(5,4) * t125;
t140 = Icges(5,1) * t127 - t163;
t96 = Icges(5,5) * t126 + t140 * t128;
t171 = (-Icges(5,2) * t127 - t163) * t128 + t96;
t162 = Icges(5,4) * t127;
t138 = -Icges(5,2) * t125 + t162;
t94 = Icges(5,6) * t126 + t138 * t128;
t170 = (-Icges(5,1) * t125 - t162) * t128 - t94;
t153 = t126 * (-Icges(5,5) * t125 - Icges(5,6) * t127) * t128;
t136 = Icges(5,5) * t127 - Icges(5,6) * t125;
t143 = -t125 * t94 + t127 * t96;
t129 = (Icges(5,3) * t128 - t136 * t126 - t143) * t126;
t154 = t125 * t126;
t103 = t123 * t127 + t124 * t154;
t155 = t124 * t128;
t152 = t126 * t127;
t104 = -t123 * t125 + t124 * t152;
t164 = Icges(5,4) * t104;
t60 = Icges(5,2) * t103 + Icges(5,6) * t155 - t164;
t100 = Icges(5,4) * t103;
t62 = -Icges(5,1) * t104 + Icges(5,5) * t155 + t100;
t144 = -t125 * t60 + t127 * t62;
t92 = Icges(5,3) * t126 + t136 * t128;
t132 = -t92 * t124 - t144;
t58 = -Icges(5,5) * t104 + Icges(5,6) * t103 + Icges(5,3) * t155;
t179 = t126 * t58;
t200 = t132 * t128 - t179;
t101 = -t123 * t154 + t124 * t127;
t156 = t123 * t128;
t102 = t123 * t152 + t124 * t125;
t165 = Icges(5,4) * t102;
t59 = Icges(5,2) * t101 - Icges(5,6) * t156 + t165;
t99 = Icges(5,4) * t101;
t61 = Icges(5,1) * t102 - Icges(5,5) * t156 + t99;
t145 = -t125 * t59 + t127 * t61;
t133 = t92 * t123 - t145;
t57 = Icges(5,5) * t102 + Icges(5,6) * t101 - Icges(5,3) * t156;
t180 = t126 * t57;
t199 = t133 * t128 - t180;
t147 = rSges(5,1) * t127 - rSges(5,2) * t125;
t98 = t126 * rSges(5,3) + t147 * t128;
t198 = (-Icges(4,6) * t123 + t202 * t124) * t126 + (Icges(4,5) * t123 + (-0.2e1 * Icges(4,4) * t128 - t205) * t124) * t128;
t197 = -(-Icges(4,6) * t124 - t202 * t123) * t126 - (Icges(4,5) * t124 + t123 * t205 + t156 * t206) * t128;
t194 = 2 * qJD(3);
t193 = m(4) / 0.2e1;
t192 = m(5) / 0.2e1;
t64 = -t104 * rSges(5,1) + t103 * rSges(5,2) + rSges(5,3) * t155;
t177 = t126 * t64;
t63 = t102 * rSges(5,1) + t101 * rSges(5,2) - rSges(5,3) * t156;
t178 = t126 * t63;
t85 = t98 * t123;
t86 = t98 * t124;
t39 = (-t128 * t85 + t178) * t124 + (t128 * t86 + t177) * t123;
t97 = rSges(5,3) * t128 - t147 * t126;
t142 = -t126 * t98 + t128 * t97;
t43 = t142 * t123 + t126 * t85 + t128 * t63;
t44 = t142 * t124 + t126 * t86 - t128 * t64;
t46 = (-t123 * t64 - t124 * t63) * t128;
t53 = t98 * t156 + t178;
t54 = t98 * t155 - t177;
t191 = m(5) * (t39 * t46 + t43 * t53 + t44 * t54);
t190 = -t123 / 0.2e1;
t189 = t123 / 0.2e1;
t188 = -t124 / 0.2e1;
t187 = t124 / 0.2e1;
t186 = -Icges(5,1) * t101 + t165 + t59;
t185 = -Icges(5,1) * t103 - t164 + t60;
t184 = -Icges(5,2) * t102 + t61 + t99;
t183 = Icges(5,2) * t104 + t100 + t62;
t182 = m(5) * qJD(4);
t176 = t126 * t92;
t117 = -t126 * pkin(3) + pkin(5) * t128;
t169 = t117 + t97;
t118 = pkin(3) * t128 + t126 * pkin(5);
t168 = t118 + t98;
t150 = qJD(4) * t128;
t67 = Icges(5,5) * t101 - Icges(5,6) * t102;
t24 = t184 * t101 - t186 * t102 - t67 * t156;
t68 = Icges(5,5) * t103 + Icges(5,6) * t104;
t25 = t183 * t101 - t185 * t102 - t68 * t156;
t11 = t123 * t25 + t124 * t24;
t81 = t94 * t123;
t83 = t96 * t123;
t19 = t101 * t81 + t102 * t83 - t199 * t123;
t82 = t94 * t124;
t84 = t96 * t124;
t20 = -t101 * t82 - t102 * t84 - t200 * t123;
t35 = t101 * t59 + t102 * t61 - t57 * t156;
t36 = t101 * t60 + t102 * t62 - t58 * t156;
t47 = t101 * t94 + t102 * t96 - t92 * t156;
t93 = Icges(5,6) * t128 - t138 * t126;
t95 = Icges(5,5) * t128 - t140 * t126;
t3 = (t101 * t93 + t102 * t95) * t126 + t47 * t128 + (-t36 * t126 + t20 * t128) * t124 + ((t35 + t176) * t126 + (-t19 - t129) * t128) * t123;
t149 = -t11 / 0.2e1 + t3 / 0.2e1;
t26 = t184 * t103 + t186 * t104 + t67 * t155;
t27 = t183 * t103 + t185 * t104 + t68 * t155;
t12 = t123 * t27 + t124 * t26;
t21 = t103 * t81 - t104 * t83 + t199 * t124;
t22 = -t103 * t82 + t104 * t84 + t200 * t124;
t37 = t103 * t59 - t104 * t61 + t57 * t155;
t38 = t103 * t60 - t104 * t62 + t58 * t155;
t48 = t103 * t94 - t104 * t96 + t92 * t155;
t4 = (t103 * t93 - t104 * t95) * t126 + t48 * t128 + (t37 * t126 - t21 * t128) * t123 + ((-t38 - t176) * t126 + (t22 + t129) * t128) * t124;
t148 = t12 / 0.2e1 - t4 / 0.2e1;
t116 = rSges(4,1) * t128 - t126 * rSges(4,2);
t40 = t145 * t128 + t180;
t41 = t144 * t128 + t179;
t146 = -t40 * t123 + t41 * t124;
t137 = Icges(4,5) * t128 - Icges(4,6) * t126;
t130 = m(5) * (t44 * t123 - t43 * t124);
t114 = (-rSges(5,1) * t125 - rSges(5,2) * t127) * t128;
t134 = t151 * t114 * t192;
t18 = -t130 / 0.2e1 + t134;
t77 = rSges(5,1) * t101 - rSges(5,2) * t102;
t78 = rSges(5,1) * t103 + rSges(5,2) * t104;
t51 = -t123 * t77 + t124 * t78;
t28 = 0.2e1 * (t39 / 0.4e1 - t51 / 0.4e1) * m(5);
t135 = t28 * qJD(1) - t18 * qJD(2);
t106 = t137 * t124;
t105 = t137 * t123;
t76 = t168 * t124;
t75 = t169 * t124;
t74 = t168 * t123;
t73 = t169 * t123;
t66 = t151 * t116;
t56 = t114 * t155 - t126 * t78;
t55 = t114 * t156 + t126 * t77;
t52 = t143 * t128 + t176;
t49 = (-t123 * t78 - t124 * t77) * t128;
t45 = (-t118 * t124 - t86) * t124 + (-t118 * t123 - t85) * t123;
t42 = (t117 * t124 + t64) * t124 + (t117 * t123 - t63) * t123;
t33 = t126 * t68 + (-t183 * t125 - t185 * t127) * t128;
t32 = t126 * t67 + (-t184 * t125 - t186 * t127) * t128;
t31 = (t125 * t82 - t127 * t84 + t58) * t128 + t132 * t126;
t30 = (-t125 * t81 + t127 * t83 + t57) * t128 + t133 * t126;
t29 = (t39 + t51) * t192;
t17 = t130 / 0.2e1 + t134;
t16 = t42 * t51 + (t123 * t74 + t124 * t76) * t114;
t15 = t52 * t126 + t146 * t128;
t14 = t48 * t126 + (-t37 * t123 + t38 * t124) * t128;
t13 = t47 * t126 + (-t35 * t123 + t36 * t124) * t128;
t10 = t123 * t22 + t124 * t21;
t9 = t123 * t20 + t124 * t19;
t7 = (t171 * t103 - t170 * t104) * t126 + (-t26 * t123 + (t27 + t153) * t124) * t128;
t6 = (t171 * t101 + t170 * t102) * t126 + (t25 * t124 + (-t24 - t153) * t123) * t128;
t5 = (-t30 * t123 + t31 * t124 + t52) * t128 + (t129 + (-t125 * t93 + t127 * t95 + t92) * t128 - t146) * t126;
t2 = m(5) * t16 + t11 * t187 + t12 * t189;
t1 = t191 + (t3 * t190 + t4 * t187 + t15 / 0.2e1) * t128 + (t13 * t189 + t14 * t188 + t5 / 0.2e1) * t126;
t8 = [0, 0, t29 * qJD(4) + (t45 * t192 - t66 * t193) * t194, t29 * qJD(3) + t49 * t182; 0, 0, t17 * qJD(4) + ((t123 * t73 + t124 * t75) * t192 + t193 * t204) * t194, t17 * qJD(3) + (t123 * t56 - t124 * t55) * t182; -t28 * qJD(4), t18 * qJD(4), (m(5) * (t42 * t45 + t73 * t74 + t75 * t76) + m(4) * (t116 - t66) * t204 + (-t121 * t106 + (t197 * t124 + (t105 - t198) * t123) * t124 + t10) * t189 + (t122 * t105 + (t198 * t123 + (-t106 - t197) * t124) * t123 + t9) * t187) * qJD(3) + t2 * qJD(4), t2 * qJD(3) + (-t15 / 0.2e1 + t148 * t124 + t149 * t123) * t150 - t135 + (m(5) * (t42 * t49 + t46 * t51 - t55 * t76 + t56 * t74 + (t123 * t54 - t124 * t53) * t114) + t6 * t187 + t7 * t189 - t191 + (-t5 / 0.2e1 + (t32 / 0.2e1 + t14 / 0.2e1) * t124 + (t33 / 0.2e1 - t13 / 0.2e1) * t123) * t126) * qJD(4); t28 * qJD(3), -t18 * qJD(3), (((t10 / 0.2e1 + t40 / 0.2e1) * t128 + t149) * t124 + ((-t9 / 0.2e1 + t41 / 0.2e1) * t128 - t148) * t123 + ((t123 * t38 + t124 * t37) * t188 + t30 * t187 + (t123 * t36 + t124 * t35 + t31) * t189) * t126 + (t39 * t42 - t43 * t76 + t44 * t74 + t45 * t46 - t53 * t75 + t54 * t73 - t16) * m(5)) * qJD(3) + t1 * qJD(4) + t135, t1 * qJD(3) + (m(5) * (t46 * t49 + t53 * t55 + t54 * t56) + t126 ^ 2 * t153 / 0.2e1) * qJD(4) + (t6 * t190 + t7 * t187 + t126 * (-t32 * t123 + t33 * t124 + (-t171 * t125 + t170 * t127) * t126) / 0.2e1) * t150;];
Cq = t8;
