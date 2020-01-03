% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:54
% DurationCPUTime: 3.20s
% Computational Cost: add. (4128->276), mult. (3679->371), div. (0->0), fcn. (2774->8), ass. (0->163)
t130 = qJD(1) ^ 2;
t129 = cos(qJ(1));
t122 = t129 * pkin(1);
t124 = qJ(1) + pkin(6);
t119 = sin(t124);
t121 = cos(t124);
t84 = t121 * pkin(2) + t119 * qJ(3);
t173 = t84 + t122;
t211 = rSges(4,2) * sin(pkin(7));
t126 = cos(pkin(7));
t214 = rSges(4,1) * t126;
t232 = -t119 * rSges(4,3) + (t211 - t214) * t121;
t241 = -t232 + t173;
t103 = qJD(3) * t119;
t182 = qJD(4) * t121;
t123 = pkin(7) + qJ(4);
t118 = sin(t123);
t120 = cos(t123);
t210 = rSges(5,2) * t120;
t80 = rSges(5,1) * t118 + t210;
t176 = t80 * t182;
t156 = t103 - t176;
t106 = t121 * qJ(3);
t117 = pkin(3) * t126 + pkin(2);
t221 = pkin(2) - t117;
t127 = -pkin(5) - qJ(3);
t189 = t121 * t127;
t128 = sin(qJ(1));
t225 = pkin(1) * t128;
t235 = -t189 - t225;
t191 = t119 * t120;
t193 = t118 * t119;
t91 = rSges(5,2) * t193;
t202 = t121 * rSges(5,3) + t91;
t54 = rSges(5,1) * t191 - t202;
t160 = t221 * t119 - t106 + t235 - t54;
t224 = pkin(2) * t119;
t81 = -t106 + t224;
t20 = (t160 - t81) * qJD(1) + t156;
t104 = qJD(3) * t121;
t183 = qJD(4) * t119;
t111 = t119 * rSges(5,3);
t190 = t120 * t121;
t192 = t118 * t121;
t55 = rSges(5,1) * t190 - rSges(5,2) * t192 + t111;
t237 = t121 * t117 - t119 * t127 + t55;
t21 = -t80 * t183 - t104 + (t173 - t84 + t237) * qJD(1);
t65 = t80 * t119;
t66 = t80 * t121;
t240 = t20 * t65 - t21 * t66;
t239 = 0.2e1 * qJD(4);
t110 = Icges(5,4) * t120;
t147 = -Icges(5,2) * t118 + t110;
t78 = Icges(5,1) * t118 + t110;
t162 = t121 * rSges(3,1) - rSges(3,2) * t119;
t236 = t122 + t162;
t94 = t119 * t211;
t234 = t121 * rSges(4,3) - t225 + t94;
t199 = Icges(5,4) * t118;
t76 = Icges(5,2) * t120 + t199;
t150 = t118 * t76 - t120 * t78;
t75 = Icges(5,5) * t120 - Icges(5,6) * t118;
t231 = qJD(1) * t150 + t75 * qJD(4);
t74 = Icges(5,5) * t118 + Icges(5,6) * t120;
t140 = qJD(4) * t74;
t51 = Icges(5,6) * t119 + t121 * t147;
t208 = t118 * t51;
t79 = Icges(5,1) * t120 - t199;
t53 = Icges(5,5) * t119 + t121 * t79;
t151 = -t120 * t53 + t208;
t195 = Icges(5,3) * t121;
t230 = -t121 * t140 + (-t119 * t75 + t151 + t195) * qJD(1);
t196 = Icges(5,6) * t121;
t50 = Icges(5,4) * t191 - Icges(5,2) * t193 - t196;
t209 = t118 * t50;
t198 = Icges(5,5) * t121;
t89 = Icges(5,4) * t193;
t52 = Icges(5,1) * t191 - t198 - t89;
t152 = -t120 * t52 + t209;
t49 = Icges(5,3) * t119 + t121 * t75;
t194 = qJD(1) * t49;
t229 = qJD(1) * t152 - t119 * t140 + t194;
t48 = Icges(5,5) * t191 - Icges(5,6) * t193 - t195;
t16 = -t119 * t152 - t121 * t48;
t228 = (-t76 * t121 + t53) * t119 - (-Icges(5,2) * t191 + t52 - t89) * t121;
t227 = t119 / 0.2e1;
t226 = -t121 / 0.2e1;
t222 = qJD(1) / 0.2e1;
t220 = -t119 * t48 - t52 * t190;
t219 = t119 * t49 + t53 * t190;
t216 = -t76 + t79;
t215 = t78 + t147;
t213 = rSges(5,1) * t120;
t207 = t119 * t74;
t205 = t121 * t74;
t184 = qJD(1) * t121;
t204 = rSges(5,3) * t184 + qJD(1) * t91;
t203 = rSges(4,3) * t184 + qJD(1) * t94;
t28 = -t119 * t150 - t205;
t188 = t28 * qJD(1);
t187 = t75 * qJD(1);
t100 = qJ(3) * t184;
t186 = t100 + t103;
t185 = qJD(1) * t119;
t181 = qJD(1) * qJD(3);
t180 = t130 * t225;
t179 = t130 * t122;
t178 = t119 * t214;
t175 = -pkin(2) - t214;
t174 = -t178 + t234;
t170 = -t183 / 0.2e1;
t167 = t182 / 0.2e1;
t166 = -t48 + t208;
t42 = t53 * t191;
t165 = t121 * t49 - t42;
t163 = -t117 - t213;
t159 = t121 * t181 - t179;
t82 = rSges(3,1) * t119 + rSges(3,2) * t121;
t153 = -rSges(5,2) * t118 + t213;
t23 = t118 * t52 + t120 * t50;
t24 = t118 * t53 + t120 * t51;
t149 = qJD(1) * (-pkin(2) * t185 + t186) + t119 * t181 - t180;
t17 = -t193 * t51 - t165;
t144 = (t119 * t17 - t121 * t16) * qJD(4);
t18 = -t192 * t50 - t220;
t19 = -t192 * t51 + t219;
t143 = (t119 * t19 - t121 * t18) * qJD(4);
t142 = qJD(4) * t78;
t141 = qJD(4) * t76;
t138 = -t119 * t51 + t121 * t50;
t137 = (-t215 * t118 + t216 * t120) * qJD(1);
t35 = t51 * qJD(1) - t119 * t141;
t37 = t53 * qJD(1) - t119 * t142;
t134 = qJD(1) * t48 - qJD(4) * t23 - t118 * t35 + t120 * t37;
t34 = -t121 * t141 + (-t119 * t147 + t196) * qJD(1);
t36 = -t121 * t142 + (-t119 * t79 + t198) * qJD(1);
t133 = -t24 * qJD(4) - t118 * t34 + t120 * t36 + t194;
t70 = t147 * qJD(4);
t71 = t79 * qJD(4);
t132 = qJD(1) * t74 - t118 * t70 + t120 * t71 + (-t118 * t78 - t120 * t76) * qJD(4);
t131 = -t118 * t228 + t138 * t120;
t73 = qJD(1) * t81;
t72 = t153 * qJD(4);
t64 = qJD(1) * t84 - t104;
t41 = qJD(1) * t241 - t104;
t40 = t103 + (t174 - t81) * qJD(1);
t39 = -qJD(4) * t65 + (t121 * t153 + t111) * qJD(1);
t38 = -t182 * t210 + (-t118 * t182 - t120 * t185) * rSges(5,1) + t204;
t29 = -t121 * t150 + t207;
t27 = t29 * qJD(1);
t26 = (qJD(1) * t232 - t64) * qJD(1) + t159;
t25 = qJD(1) * (-qJD(1) * t178 + t203) + t149;
t22 = qJD(2) + (t119 * t54 + t121 * t55) * qJD(4);
t13 = -t72 * t182 + (-t39 - t64 + t221 * t184 + (qJD(4) * t80 + (qJ(3) + t127) * qJD(1)) * t119) * qJD(1) + t159;
t12 = -t72 * t183 + (-t176 - t100 + t38 + (-t117 * t119 - t189 + t224) * qJD(1)) * qJD(1) + t149;
t11 = t132 * t119 - t121 * t231;
t10 = t119 * t231 + t132 * t121;
t9 = -qJD(4) * t151 + t118 * t36 + t120 * t34;
t8 = -qJD(4) * t152 + t118 * t37 + t120 * t35;
t7 = (t119 * t39 + t121 * t38 + (-t119 * t55 + t121 * t54) * qJD(1)) * qJD(4);
t6 = t27 + t143;
t5 = t144 + t188;
t1 = [(-qJD(4) * t150 + t118 * t71 + t120 * t70) * qJD(1) + (t27 + ((t17 - t42 + (t49 + t209) * t121 + t220) * t121 + t219 * t119) * qJD(4)) * t167 + m(3) * ((-t130 * t82 - t180) * t236 + (-t179 + (-0.2e1 * t162 - t122 + t236) * t130) * (-t82 - t225)) + (t5 - t188 + ((t121 * t166 + t19 - t219) * t121 + (t119 * t166 + t165 + t18) * t119) * qJD(4)) * t170 + (t9 + t10) * t183 / 0.2e1 + (t13 * (t119 * t163 + t202 + t235) + t20 * t104 + t12 * (t122 + t237) + t21 * (t103 + t204) + t240 * qJD(4) + ((-t128 * t21 - t20 * t129) * pkin(1) + (t20 * (-t117 - t153) - t21 * t127) * t121 + (t20 * (-rSges(5,3) + t127) + t21 * t163) * t119) * qJD(1) - (qJD(1) * t160 + t156 - t20 - t73) * t21) * m(5) + (t26 * (t119 * t175 + t106 + t234) + t40 * t104 + t25 * t241 + t41 * (t186 + t203) + ((-t128 * t41 - t40 * t129) * pkin(1) + t40 * (t175 + t211) * t121 + (t40 * (-rSges(4,3) - qJ(3)) + t41 * t175) * t119) * qJD(1) - (qJD(1) * t174 + t103 - t40 - t73) * t41) * m(4) - (t8 + t11 + t6) * t182 / 0.2e1 + ((t23 + t28) * t119 + (t24 + t29) * t121) * qJD(4) * t222; m(5) * t7; 0.2e1 * (t12 * t226 + t13 * t227) * m(5) + 0.2e1 * (t25 * t226 + t26 * t227) * m(4); (t9 * t119 - t8 * t121 + (t23 * t119 + t121 * t24) * qJD(1)) * t222 + ((-t183 * t205 + t187) * t119 + (t137 + (t119 * t207 + t131) * qJD(4)) * t121) * t170 + ((-t182 * t207 - t187) * t121 + (t137 + (t121 * t205 + t131) * qJD(4)) * t119) * t167 - qJD(1) * ((t216 * t118 + t215 * t120) * qJD(1) + (t138 * t118 + t120 * t228) * qJD(4)) / 0.2e1 + (qJD(1) * t10 + (-(t119 * t229 + t134 * t121) * t121 + t119 * (t119 * t230 + t133 * t121) + (t18 * t119 + t19 * t121) * qJD(1)) * t239) * t227 + (qJD(1) * t11 + (t119 * (t133 * t119 - t121 * t230) - t121 * (t134 * t119 - t121 * t229) + (t16 * t119 + t17 * t121) * qJD(1)) * t239) * t226 + (t5 + t144) * t185 / 0.2e1 + (t6 + t143) * t184 / 0.2e1 + ((-t20 * t72 + t7 * t55 + t22 * (qJD(1) * t54 + t38) + (-qJD(1) * t21 - t13) * t80) * t121 + (-t21 * t72 + t7 * t54 + t22 * (-qJD(1) * t55 + t39) + (qJD(1) * t20 - t12) * t80) * t119 - t240 * qJD(1) - (t22 * (-t119 * t65 - t121 * t66) + (-t21 * t119 - t20 * t121) * t153) * qJD(4)) * m(5);];
tauc = t1(:);
