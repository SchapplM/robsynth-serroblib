% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:47
% DurationCPUTime: 3.52s
% Computational Cost: add. (2808->263), mult. (6623->357), div. (0->0), fcn. (6808->6), ass. (0->153)
t126 = sin(qJ(4));
t127 = cos(qJ(4));
t198 = Icges(5,4) * t127;
t159 = -Icges(5,2) * t126 + t198;
t200 = sin(pkin(6));
t201 = cos(pkin(6));
t227 = sin(qJ(1));
t228 = cos(qJ(1));
t88 = -t200 * t227 - t201 * t228;
t89 = t228 * t200 - t227 * t201;
t51 = Icges(5,6) * t88 + t159 * t89;
t199 = Icges(5,4) * t126;
t161 = Icges(5,1) * t127 - t199;
t54 = Icges(5,5) * t88 + t161 * t89;
t247 = t126 * t51 - t127 * t54;
t157 = Icges(5,5) * t127 - Icges(5,6) * t126;
t48 = Icges(5,3) * t88 + t157 * t89;
t14 = -t247 * t89 + t48 * t88;
t50 = Icges(5,3) * t89 - t157 * t88;
t252 = t89 * t50;
t167 = rSges(5,1) * t126 + rSges(5,2) * t127;
t243 = qJD(4) * t167;
t251 = t88 * t243;
t250 = t89 * t243;
t186 = t227 * pkin(2);
t215 = rSges(5,2) * t126;
t168 = rSges(5,1) * t127 - t215;
t57 = t88 * rSges(5,3) + t168 * t89;
t249 = t89 * pkin(3) + t88 * pkin(5) - t186 + t57;
t100 = Icges(5,1) * t126 + t198;
t158 = Icges(5,2) * t127 + t199;
t162 = -t127 * t100 + t126 * t158;
t156 = Icges(5,5) * t126 + Icges(5,6) * t127;
t63 = t156 * t88;
t235 = t162 * t89 - t63;
t248 = t235 * qJD(1);
t23 = -t126 * t54 - t127 * t51;
t62 = t156 * t89;
t124 = t228 * pkin(2);
t193 = t228 * pkin(1) + t227 * qJ(2);
t238 = t124 + t193;
t246 = -t88 * pkin(3) + t238;
t120 = t228 * qJ(2);
t187 = t227 * pkin(1);
t103 = t187 - t120;
t117 = qJD(2) * t227;
t176 = qJD(1) * t228;
t194 = qJ(2) * t176 + t117;
t245 = qJD(1) * t103 - t117 + t194;
t244 = -t88 * rSges(4,1) - t89 * rSges(4,2) + t238;
t82 = t88 * qJD(1);
t83 = qJD(1) * t89;
t242 = t83 * pkin(3) + t82 * pkin(5);
t118 = qJD(2) * t228;
t241 = -qJD(1) * t238 + t118;
t240 = 0.2e1 * qJD(4);
t139 = t89 * rSges(4,1) - t88 * rSges(4,2) - t186;
t239 = -t103 + t139;
t144 = t228 * rSges(3,1) + t227 * rSges(3,3);
t237 = t193 + t144;
t56 = Icges(5,5) * t89 - t161 * t88;
t234 = t54 * t88 + t56 * t89;
t128 = qJD(1) ^ 2;
t233 = t82 / 0.2e1;
t232 = t83 / 0.2e1;
t229 = t83 * pkin(5);
t226 = t83 * rSges(5,3);
t208 = t127 * t56;
t223 = -t89 * t208 - t88 * t50;
t206 = t127 * t88;
t222 = -t56 * t206 + t252;
t207 = t127 * t83;
t219 = rSges(5,1) * t207 + t82 * rSges(5,3);
t218 = t83 * rSges(4,1) - t82 * rSges(4,2);
t217 = -rSges(5,1) * t206 + t89 * rSges(5,3);
t53 = Icges(5,6) * t89 - t159 * t88;
t211 = t126 * t53;
t210 = t126 * t83;
t209 = t126 * t88;
t205 = -t100 - t159;
t204 = t161 - t158;
t175 = qJD(1) * t227;
t203 = (-pkin(1) * t175 + t117 + t194) * qJD(1);
t197 = qJD(4) * t89;
t196 = t88 * qJD(4);
t195 = t157 * qJD(1);
t191 = qJD(4) * t127;
t190 = t126 * qJD(4);
t189 = -t228 / 0.2e1;
t188 = t227 / 0.2e1;
t185 = t227 * rSges(3,1);
t184 = t88 * t190;
t180 = -t197 / 0.2e1;
t177 = t196 / 0.2e1;
t171 = t82 * rSges(4,1) + t83 * rSges(4,2);
t21 = t251 + t117 + (-t103 + t249) * qJD(1);
t59 = rSges(5,2) * t209 + t217;
t22 = t250 - t118 + (pkin(5) * t89 + t246 + t59) * qJD(1);
t170 = -t21 * t88 - t22 * t89;
t169 = -t57 * t89 + t59 * t88;
t163 = -t208 + t211;
t16 = t206 * t54 - t209 * t51 - t48 * t89;
t114 = qJD(1) * t118;
t155 = -t124 * t128 + t114;
t154 = pkin(3) + t168;
t153 = t126 * t82 + t191 * t89;
t152 = -t127 * t82 + t190 * t89;
t151 = t191 * t88 - t210;
t150 = t184 + t207;
t15 = t211 * t89 + t223;
t149 = qJD(4) * (-t14 * t88 + t15 * t89);
t17 = t209 * t53 + t222;
t148 = (-t16 * t88 + t17 * t89) * qJD(4);
t147 = -t128 * t186 + t203;
t146 = -t187 - t186;
t145 = -t185 - t187;
t141 = t51 * t88 + t53 * t89;
t31 = rSges(5,1) * t152 + rSges(5,2) * t153 + t226;
t32 = rSges(5,1) * t184 + rSges(5,2) * t151 + t219;
t137 = t31 * t89 + t32 * t88 - t57 * t82 - t59 * t83;
t136 = (t126 * t205 + t127 * t204) * qJD(1);
t91 = t159 * qJD(4);
t92 = t161 * qJD(4);
t130 = -t126 * t91 + t127 * t92 + (-t100 * t126 - t127 * t158) * qJD(4);
t129 = t126 * t234 + t141 * t127;
t122 = t228 * rSges(3,3);
t116 = rSges(3,3) * t176;
t93 = t168 * qJD(4);
t90 = t157 * qJD(4);
t81 = qJD(1) * t193 - t118;
t69 = t167 * t88;
t68 = t167 * t89;
t61 = -qJD(1) * t81 - t128 * t144 + t114;
t60 = qJD(1) * (-rSges(3,1) * t175 + t116) + t203;
t46 = qJD(1) * t239 + t117;
t37 = (t171 - t81) * qJD(1) + t155;
t36 = qJD(1) * t218 + t147;
t35 = -t162 * t88 - t62;
t33 = t35 * qJD(1);
t26 = Icges(5,5) * t150 + Icges(5,6) * t151 + Icges(5,3) * t82;
t25 = Icges(5,5) * t152 + Icges(5,6) * t153 + Icges(5,3) * t83;
t24 = t126 * t56 + t127 * t53;
t20 = qJD(4) * t169 - qJD(3);
t13 = (-t167 * t83 + t88 * t93) * qJD(4) + (t82 * pkin(3) - t229 - t31 - t81) * qJD(1) + t155;
t12 = (t167 * t82 + t89 * t93) * qJD(4) + (t242 + t32) * qJD(1) + t147;
t11 = t130 * t88 - t156 * t82 + t162 * t83 - t89 * t90;
t10 = t130 * t89 - t156 * t83 - t162 * t82 + t88 * t90;
t9 = qJD(4) * t163 - t126 * (Icges(5,1) * t150 + Icges(5,4) * t151 + Icges(5,5) * t82) - t127 * (Icges(5,4) * t150 + Icges(5,2) * t151 + Icges(5,6) * t82);
t8 = -qJD(4) * t247 - t126 * (Icges(5,1) * t152 + Icges(5,4) * t153 + Icges(5,5) * t83) - t127 * (Icges(5,4) * t152 + Icges(5,2) * t153 + Icges(5,6) * t83);
t7 = t137 * qJD(4);
t6 = t33 + t148;
t5 = t149 - t248;
t1 = [(t126 * t92 + t127 * t91) * qJD(1) + t33 * t177 + (t9 + t11) * t197 / 0.2e1 + (t13 * (t120 + t146) + t12 * (t217 + t246) + (t13 * (rSges(5,3) + pkin(5)) + t12 * t215) * t88 + (t12 * pkin(5) + t13 * t154) * t89 + (t154 * t82 - t226 - t229 + t241 - t250) * t21 + (-rSges(5,2) * t210 - t167 * t196 + t251 + t21 + t219 + (t146 - t249) * qJD(1) + t242 + t245) * t22) * m(5) + (t37 * t239 + t36 * t244 + (t171 + t241) * t46 + (t218 + t46 + (-t139 + t146) * qJD(1) + t245) * (t244 * qJD(1) - t118)) * m(4) + (t61 * (t120 + t122 + t145) + t60 * t237 + (t116 + (t185 - t122 + t145) * qJD(1) + t245) * (t237 * qJD(1) - t118)) * m(3) + (t5 + t248) * t180 - (t8 + t10 + t6) * t196 / 0.2e1 + (-t162 * qJD(1) + ((t15 - t16 - t223) * t88 + t222 * t89) * t177 + (-t24 + t35) * t233 + (-t23 - t235) * t232 + ((t16 + (-t163 + t48) * t89) * t89 + (t17 + (t48 - t211) * t88 + t252 - t222) * t88) * t180) * qJD(4); 0.2e1 * (t12 * t189 + t13 * t188) * m(5) + 0.2e1 * (t188 * t37 + t189 * t36) * m(4) + 0.2e1 * (t188 * t61 + t189 * t60) * m(3); -m(5) * t7; qJD(1) * (-t23 * t83 - t24 * t82 - t8 * t88 + t89 * t9) / 0.2e1 + ((t197 * t63 - t195) * t89 + (t136 + (-t89 * t62 + t129) * qJD(4)) * t88) * t180 + ((t196 * t62 + t195) * t88 + (t136 + (-t88 * t63 + t129) * qJD(4)) * t89) * t177 - qJD(1) * ((t204 * t126 - t205 * t127) * qJD(1) + (t141 * t126 - t127 * t234) * qJD(4)) / 0.2e1 - (qJD(1) * t10 + (-(-t247 * t82 - t25 * t88 - t48 * t83) * t88 + t14 * t83 + t15 * t82 + (t163 * t82 - t26 * t88 + t50 * t83) * t89) * t240) * t88 / 0.2e1 + (qJD(1) * t11 + (t16 * t83 + t17 * t82 - (t247 * t83 + t25 * t89 - t48 * t82) * t88 + (-t163 * t83 + t26 * t89 + t50 * t82) * t89) * t240) * t89 / 0.2e1 + (t6 + t148) * t233 + (t5 + t149) * t232 + (t7 * t169 + t20 * t137 - t170 * t93 - (-t12 * t89 - t13 * t88 + t21 * t83 - t22 * t82) * t167 - (-t21 * t68 + t22 * t69) * qJD(1) - (t20 * (t68 * t89 + t69 * t88) - t170 * t168) * qJD(4)) * m(5);];
tauc = t1(:);
