% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:37
% DurationCPUTime: 1.94s
% Computational Cost: add. (5389->229), mult. (4070->320), div. (0->0), fcn. (3160->6), ass. (0->149)
t126 = pkin(7) + qJ(2);
t122 = sin(t126);
t208 = pkin(2) * qJD(2);
t179 = t122 * t208;
t127 = qJD(2) + qJD(3);
t124 = qJ(3) + t126;
t120 = sin(t124);
t121 = cos(t124);
t80 = rSges(4,1) * t120 + rSges(4,2) * t121;
t206 = t127 * t80;
t65 = -t179 - t206;
t230 = 2 * qJD(4);
t129 = cos(qJ(4));
t190 = t121 * t129;
t229 = rSges(5,1) * t190 + t120 * rSges(5,3);
t228 = -t121 * pkin(3) - t120 * pkin(6);
t125 = Icges(5,4) * t129;
t128 = sin(qJ(4));
t155 = -Icges(5,2) * t128 + t125;
t106 = Icges(5,1) * t128 + t125;
t194 = t120 * t128;
t100 = rSges(5,2) * t194;
t209 = rSges(5,2) * t129;
t177 = qJD(4) * t209;
t192 = t121 * t127;
t153 = rSges(5,3) * t192 + t127 * t100 - t121 * t177;
t184 = qJD(4) * t128;
t175 = t121 * t184;
t189 = t121 * rSges(5,3) + t100;
t193 = t120 * t129;
t63 = rSges(5,1) * t193 - t189;
t56 = t127 * t63;
t117 = t121 * pkin(6);
t82 = pkin(3) * t120 - t117;
t96 = pkin(6) * t192;
t227 = -rSges(5,1) * t175 + t127 * t82 + t153 + t56 + t96;
t191 = t121 * t128;
t180 = rSges(5,2) * t191;
t64 = -t180 + t229;
t145 = t64 - t228;
t108 = rSges(5,1) * t128 + t209;
t186 = qJD(4) * t120;
t78 = t108 * t186;
t226 = t127 * t145 - t78;
t103 = Icges(5,5) * t129 - Icges(5,6) * t128;
t200 = Icges(5,4) * t128;
t104 = Icges(5,2) * t129 + t200;
t154 = t128 * t104 - t129 * t106;
t225 = t103 * qJD(4) + t127 * t154;
t102 = Icges(5,5) * t128 + Icges(5,6) * t129;
t139 = Icges(5,3) * t127 - qJD(4) * t102;
t148 = t155 * t121;
t60 = Icges(5,6) * t120 + t148;
t204 = t128 * t60;
t107 = Icges(5,1) * t129 - t200;
t149 = t107 * t121;
t62 = Icges(5,5) * t120 + t149;
t157 = -t129 * t62 + t204;
t195 = t120 * t127;
t224 = -t103 * t195 + t121 * t139 + t127 * t157;
t147 = t103 * t121;
t59 = Icges(5,4) * t193 - Icges(5,2) * t194 - Icges(5,6) * t121;
t205 = t128 * t59;
t99 = Icges(5,4) * t194;
t61 = Icges(5,1) * t193 - Icges(5,5) * t121 - t99;
t158 = -t129 * t61 + t205;
t223 = t120 * t139 + (t147 + t158) * t127;
t57 = Icges(5,5) * t193 - Icges(5,6) * t194 - Icges(5,3) * t121;
t24 = -t120 * t158 - t121 * t57;
t212 = -Icges(5,2) * t193 + t61 - t99;
t214 = t106 * t120 + t59;
t222 = -t128 * t212 - t129 * t214;
t219 = t127 / 0.2e1;
t218 = pkin(2) * t122;
t217 = pkin(2) * qJD(2) ^ 2;
t216 = -t120 * t57 - t61 * t190;
t58 = Icges(5,3) * t120 + t147;
t215 = t120 * t58 + t62 * t190;
t213 = -t106 * t121 - t60;
t211 = -t104 * t121 + t62;
t210 = rSges(5,1) * t129;
t185 = qJD(4) * t121;
t176 = t108 * t185;
t146 = -t176 - t179;
t30 = (-t63 - t82) * t127 + t146;
t203 = t30 * t121;
t197 = t102 * t121;
t47 = -t120 * t154 - t197;
t202 = t47 * t127;
t198 = t102 * t120;
t196 = t103 * t127;
t188 = -t104 + t107;
t187 = t106 + t155;
t183 = t127 * t180 + (rSges(5,1) * t184 + t177) * t120;
t182 = t122 * t217;
t123 = cos(t126);
t181 = t123 * t217;
t178 = t123 * t208;
t172 = -pkin(3) - t210;
t171 = -t186 / 0.2e1;
t168 = t185 / 0.2e1;
t49 = t62 * t193;
t166 = t121 * t58 - t49;
t44 = (-t127 * t193 - t175) * rSges(5,1) + t153;
t165 = t44 + t56;
t45 = t127 * t229 - t183;
t164 = -t127 * t64 + t45;
t163 = -t57 + t204;
t81 = t121 * rSges(4,1) - rSges(4,2) * t120;
t68 = rSges(4,1) * t192 - rSges(4,2) * t195;
t159 = -rSges(5,2) * t128 + t210;
t35 = t128 * t61 + t129 * t59;
t36 = t128 * t62 + t129 * t60;
t25 = -t194 * t60 - t166;
t151 = (t120 * t25 - t121 * t24) * qJD(4);
t26 = -t191 * t59 - t216;
t27 = -t191 * t60 + t215;
t150 = (t120 * t27 - t121 * t26) * qJD(4);
t144 = -t128 * t211 + t129 * t213;
t143 = t120 * t172 + t117 + t189;
t142 = (-t128 * t187 + t129 * t188) * t127;
t141 = Icges(5,5) * t127 - qJD(4) * t106;
t140 = Icges(5,6) * t127 - qJD(4) * t104;
t48 = -t121 * t154 + t198;
t46 = t48 * t127;
t10 = t46 + t150;
t14 = -qJD(4) * t158 + t128 * (t120 * t141 + t127 * t149) + t129 * (t120 * t140 + t127 * t148);
t15 = -qJD(4) * t157 + t128 * (-t107 * t195 + t121 * t141) + t129 * (t121 * t140 - t155 * t195);
t89 = t155 * qJD(4);
t90 = t107 * qJD(4);
t132 = t102 * t127 - t128 * t89 + t129 * t90 + (-t104 * t129 - t106 * t128) * qJD(4);
t18 = t120 * t225 + t132 * t121;
t19 = t132 * t120 - t121 * t225;
t9 = t151 + t202;
t138 = (t46 + ((t25 - t49 + (t58 + t205) * t121 + t216) * t121 + t215 * t120) * qJD(4)) * t168 + (-qJD(4) * t154 + t128 * t90 + t129 * t89) * t127 + (t9 - t202 + ((t121 * t163 - t215 + t27) * t121 + (t120 * t163 + t166 + t26) * t120) * qJD(4)) * t171 + (t15 + t18) * t186 / 0.2e1 - (t10 + t14 + t19) * t185 / 0.2e1 + ((t35 + t47) * t120 + (t36 + t48) * t121) * qJD(4) * t219;
t31 = t178 + t226;
t131 = (t172 * t203 + (t30 * (-rSges(5,3) - pkin(6)) + t31 * t172) * t120) * t127;
t119 = pkin(2) * t123;
t91 = t159 * qJD(4);
t76 = t108 * t121;
t75 = t108 * t120;
t66 = t127 * t81 + t178;
t55 = -t127 * t68 - t181;
t54 = -t127 * t206 - t182;
t32 = qJD(1) + (t120 * t63 + t121 * t64) * qJD(4);
t23 = -t181 - t91 * t185 + (t127 * t228 - t45 + t78) * t127;
t22 = -t182 - t91 * t186 + (-pkin(3) * t195 - t176 + t44 + t96) * t127;
t11 = (t164 * t120 + t165 * t121) * qJD(4);
t1 = [m(5) * t11; m(4) * (t55 * (-t80 - t218) + t54 * (t119 + t81) + (-t68 - t178 + t66) * t65) + t138 + (t23 * (t143 - t218) + t30 * (-t178 + t183) + t22 * (t119 + t145) + t131 + (-t146 + t30 - t179 + t227) * t31) * m(5); t138 + (t23 * t143 + t22 * t145 + t131 + (t176 + t227) * t31 + (t183 + t226) * t30) * m(5) + (-(-t65 * t81 - t66 * t80) * t127 + t54 * t81 - t55 * t80 - t65 * t68 - t66 * t206) * m(4); ((t127 * t36 - t14) * t121 + (t127 * t35 + t15) * t120) * t219 + ((-t186 * t197 + t196) * t120 + (t142 + (-t222 * t121 + (t198 + t144) * t120) * qJD(4)) * t121) * t171 + ((-t185 * t198 - t196) * t121 + (t142 + (t144 * t120 + (-t222 + t197) * t121) * qJD(4)) * t120) * t168 - t127 * ((t188 * t128 + t187 * t129) * t127 + ((t211 * t120 - t212 * t121) * t129 + (t120 * t213 + t121 * t214) * t128) * qJD(4)) / 0.2e1 + (t18 * t127 + ((-t120 * t223 + t127 * t27) * t121 + (t120 * t224 + t127 * t26) * t120) * t230) * t120 / 0.2e1 - (t19 * t127 + ((t121 * t223 + t127 * t25) * t121 + (-t121 * t224 + t127 * t24) * t120) * t230) * t121 / 0.2e1 + (t9 + t151) * t195 / 0.2e1 + (t10 + t150) * t192 / 0.2e1 + ((t11 * t64 + t32 * t165 - t30 * t91) * t121 + (t11 * t63 + t32 * t164 - t31 * t91) * t120 + ((-t127 * t31 - t23) * t121 + (t127 * t30 - t22) * t120) * t108 - (t30 * t75 - t31 * t76) * t127 - (t32 * (-t120 * t75 - t121 * t76) + (-t31 * t120 - t203) * t159) * qJD(4)) * m(5);];
tauc = t1(:);
