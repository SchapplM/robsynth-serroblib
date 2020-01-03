% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:47
% EndTime: 2019-12-31 16:20:51
% DurationCPUTime: 2.14s
% Computational Cost: add. (4050->254), mult. (3505->350), div. (0->0), fcn. (2694->6), ass. (0->151)
t119 = pkin(6) + qJ(2);
t115 = sin(t119);
t100 = qJD(3) * t115;
t117 = cos(t119);
t165 = t117 * qJD(4);
t118 = pkin(7) + qJ(4);
t114 = sin(t118);
t116 = cos(t118);
t191 = rSges(5,2) * t116;
t79 = rSges(5,1) * t114 + t191;
t161 = t79 * t165;
t146 = t100 - t161;
t103 = t117 * qJ(3);
t122 = -pkin(5) - qJ(3);
t171 = t117 * t122;
t121 = cos(pkin(7));
t113 = pkin(3) * t121 + pkin(2);
t205 = pkin(2) - t113;
t173 = t115 * t116;
t175 = t114 * t115;
t90 = rSges(5,2) * t175;
t184 = t117 * rSges(5,3) + t90;
t54 = rSges(5,1) * t173 - t184;
t202 = t205 * t115 - t103 - t171 - t54;
t208 = pkin(2) * t115;
t80 = -t103 + t208;
t20 = (-t80 + t202) * qJD(2) + t146;
t101 = qJD(3) * t117;
t166 = qJD(4) * t115;
t108 = t115 * rSges(5,3);
t172 = t116 * t117;
t174 = t114 * t117;
t55 = rSges(5,1) * t172 - rSges(5,2) * t174 + t108;
t216 = t117 * t113 - t115 * t122 + t55;
t21 = t216 * qJD(2) - t79 * t166 - t101;
t65 = t79 * t115;
t66 = t79 * t117;
t220 = t20 * t65 - t21 * t66;
t219 = 0.2e1 * qJD(4);
t192 = rSges(4,2) * sin(pkin(7));
t194 = rSges(4,1) * t121;
t215 = -t115 * rSges(4,3) + (t192 - t194) * t117;
t83 = t117 * pkin(2) + t115 * qJ(3);
t217 = -t215 + t83;
t107 = Icges(5,4) * t116;
t139 = -Icges(5,2) * t114 + t107;
t77 = Icges(5,1) * t114 + t107;
t181 = Icges(5,4) * t114;
t75 = Icges(5,2) * t116 + t181;
t141 = t114 * t75 - t116 * t77;
t74 = Icges(5,5) * t116 - Icges(5,6) * t114;
t214 = t141 * qJD(2) + t74 * qJD(4);
t73 = Icges(5,5) * t114 + Icges(5,6) * t116;
t132 = qJD(4) * t73;
t51 = Icges(5,6) * t115 + t139 * t117;
t189 = t114 * t51;
t78 = Icges(5,1) * t116 - t181;
t53 = Icges(5,5) * t115 + t78 * t117;
t142 = -t116 * t53 + t189;
t177 = Icges(5,3) * t117;
t213 = -t117 * t132 + (-t115 * t74 + t142 + t177) * qJD(2);
t178 = Icges(5,6) * t117;
t50 = Icges(5,4) * t173 - Icges(5,2) * t175 - t178;
t190 = t114 * t50;
t180 = Icges(5,5) * t117;
t88 = Icges(5,4) * t175;
t52 = Icges(5,1) * t173 - t180 - t88;
t143 = -t116 * t52 + t190;
t49 = Icges(5,3) * t115 + t117 * t74;
t176 = qJD(2) * t49;
t212 = t143 * qJD(2) - t115 * t132 + t176;
t48 = Icges(5,5) * t173 - Icges(5,6) * t175 - t177;
t16 = -t143 * t115 - t117 * t48;
t211 = (-t75 * t117 + t53) * t115 - (-Icges(5,2) * t173 + t52 - t88) * t117;
t210 = t115 / 0.2e1;
t209 = -t117 / 0.2e1;
t206 = qJD(2) / 0.2e1;
t204 = -t115 * t48 - t52 * t172;
t203 = t115 * t49 + t53 * t172;
t164 = qJD(2) * qJD(3);
t168 = qJD(2) * t115;
t167 = qJD(2) * t117;
t97 = qJ(3) * t167;
t185 = t100 + t97;
t199 = qJD(2) * (-pkin(2) * t168 + t185) + t115 * t164;
t198 = -t75 + t78;
t197 = t77 + t139;
t196 = rSges(5,3) * t167 + qJD(2) * t90;
t93 = t115 * t192;
t195 = rSges(4,3) * t167 + qJD(2) * t93;
t193 = rSges(5,1) * t116;
t188 = t115 * t73;
t186 = t117 * t73;
t183 = t117 * rSges(4,3) + t93;
t26 = -t141 * t115 - t186;
t170 = t26 * qJD(2);
t169 = t74 * qJD(2);
t163 = t115 * t194;
t160 = -pkin(2) - t194;
t157 = -t166 / 0.2e1;
t154 = t165 / 0.2e1;
t153 = -t48 + t189;
t42 = t53 * t173;
t152 = t117 * t49 - t42;
t150 = -t113 - t193;
t144 = -rSges(5,2) * t114 + t193;
t23 = t114 * t52 + t116 * t50;
t24 = t114 * t53 + t116 * t51;
t17 = -t51 * t175 - t152;
t136 = (t115 * t17 - t117 * t16) * qJD(4);
t18 = -t50 * t174 - t204;
t19 = -t51 * t174 + t203;
t135 = (t115 * t19 - t117 * t18) * qJD(4);
t134 = t77 * qJD(4);
t133 = qJD(4) * t75;
t131 = -t51 * t115 + t50 * t117;
t130 = (-t197 * t114 + t198 * t116) * qJD(2);
t33 = t51 * qJD(2) - t115 * t133;
t35 = t53 * qJD(2) - t115 * t134;
t127 = qJD(2) * t48 - t23 * qJD(4) - t114 * t33 + t116 * t35;
t32 = -t117 * t133 + (-t139 * t115 + t178) * qJD(2);
t34 = -t117 * t134 + (-t78 * t115 + t180) * qJD(2);
t126 = -t24 * qJD(4) - t114 * t32 + t116 * t34 + t176;
t68 = t139 * qJD(4);
t69 = t78 * qJD(4);
t125 = qJD(2) * t73 - t114 * t68 + t116 * t69 + (-t114 * t77 - t116 * t75) * qJD(4);
t124 = -t211 * t114 + t131 * t116;
t96 = t117 * t164;
t72 = qJD(2) * t80;
t70 = t144 * qJD(4);
t64 = qJD(2) * t83 - t101;
t56 = t163 - t183;
t41 = t217 * qJD(2) - t101;
t40 = t100 + (-t56 - t80) * qJD(2);
t37 = -qJD(4) * t65 + (t144 * t117 + t108) * qJD(2);
t36 = -t165 * t191 + (-t114 * t165 - t116 * t168) * rSges(5,1) + t196;
t29 = t96 + (t215 * qJD(2) - t64) * qJD(2);
t28 = qJD(2) * (-qJD(2) * t163 + t195) + t199;
t27 = -t141 * t117 + t188;
t25 = t27 * qJD(2);
t22 = qJD(1) + (t115 * t54 + t117 * t55) * qJD(4);
t13 = -t70 * t165 + t96 + (-t37 - t64 + t205 * t167 + (qJD(4) * t79 + (qJ(3) + t122) * qJD(2)) * t115) * qJD(2);
t12 = -t70 * t166 + (-t161 + t36 - t97 + (-t113 * t115 - t171 + t208) * qJD(2)) * qJD(2) + t199;
t11 = t125 * t115 - t214 * t117;
t10 = t214 * t115 + t125 * t117;
t9 = -t142 * qJD(4) + t114 * t34 + t116 * t32;
t8 = -t143 * qJD(4) + t114 * t35 + t116 * t33;
t7 = (t115 * t37 + t117 * t36 + (-t115 * t55 + t117 * t54) * qJD(2)) * qJD(4);
t6 = t25 + t135;
t5 = t136 + t170;
t1 = [m(5) * t7; (-t141 * qJD(4) + t114 * t69 + t116 * t68) * qJD(2) + (t25 + ((t17 - t42 + (t49 + t190) * t117 + t204) * t117 + t203 * t115) * qJD(4)) * t154 + (t5 - t170 + ((t153 * t117 + t19 - t203) * t117 + (t153 * t115 + t152 + t18) * t115) * qJD(4)) * t157 + (t9 + t10) * t166 / 0.2e1 + (t13 * (t150 * t115 - t171 + t184) + t20 * t101 + t12 * t216 + t21 * (t100 + t196) + t220 * qJD(4) + ((t20 * (-t113 - t144) - t21 * t122) * t117 + (t20 * (-rSges(5,3) + t122) + t21 * t150) * t115) * qJD(2) - (t202 * qJD(2) + t146 - t20 - t72) * t21) * m(5) + (t29 * (t160 * t115 + t103 + t183) + t40 * t101 + t28 * t217 + t41 * (t185 + t195) + (t40 * (t160 + t192) * t117 + (t40 * (-rSges(4,3) - qJ(3)) + t41 * t160) * t115) * qJD(2) - (-qJD(2) * t56 + t100 - t40 - t72) * t41) * m(4) - (t8 + t11 + t6) * t165 / 0.2e1 + ((t23 + t26) * t115 + (t24 + t27) * t117) * qJD(4) * t206; 0.2e1 * (t12 * t209 + t13 * t210) * m(5) + 0.2e1 * (t28 * t209 + t29 * t210) * m(4); (t9 * t115 - t8 * t117 + (t23 * t115 + t117 * t24) * qJD(2)) * t206 + ((-t166 * t186 + t169) * t115 + (t130 + (t115 * t188 + t124) * qJD(4)) * t117) * t157 + ((-t165 * t188 - t169) * t117 + (t130 + (t117 * t186 + t124) * qJD(4)) * t115) * t154 - qJD(2) * ((t198 * t114 + t197 * t116) * qJD(2) + (t131 * t114 + t211 * t116) * qJD(4)) / 0.2e1 + (t10 * qJD(2) + (-(t212 * t115 + t127 * t117) * t117 + t115 * (t213 * t115 + t126 * t117) + (t18 * t115 + t19 * t117) * qJD(2)) * t219) * t210 + (t11 * qJD(2) + (t115 * (t126 * t115 - t213 * t117) - t117 * (t127 * t115 - t212 * t117) + (t16 * t115 + t17 * t117) * qJD(2)) * t219) * t209 + (t5 + t136) * t168 / 0.2e1 + (t6 + t135) * t167 / 0.2e1 + ((-t20 * t70 + t7 * t55 + t22 * (qJD(2) * t54 + t36) + (-qJD(2) * t21 - t13) * t79) * t117 + (-t21 * t70 + t7 * t54 + t22 * (-qJD(2) * t55 + t37) + (qJD(2) * t20 - t12) * t79) * t115 - t220 * qJD(2) - (t22 * (-t115 * t65 - t117 * t66) + (-t21 * t115 - t20 * t117) * t144) * qJD(4)) * m(5);];
tauc = t1(:);
