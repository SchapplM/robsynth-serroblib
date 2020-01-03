% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR7
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:39
% DurationCPUTime: 2.80s
% Computational Cost: add. (2757->284), mult. (4089->393), div. (0->0), fcn. (3018->6), ass. (0->161)
t247 = 2 * qJD(4);
t132 = pkin(6) + qJ(4);
t117 = sin(t132);
t118 = cos(t132);
t164 = rSges(5,1) * t117 + rSges(5,2) * t118;
t136 = sin(qJ(1));
t202 = t118 * t136;
t100 = Icges(5,4) * t202;
t137 = cos(qJ(1));
t203 = t117 * t136;
t209 = Icges(5,5) * t137;
t60 = Icges(5,1) * t203 + t100 + t209;
t210 = Icges(5,4) * t118;
t158 = Icges(5,1) * t117 + t210;
t61 = -Icges(5,5) * t136 + t137 * t158;
t85 = -Icges(5,2) * t117 + t210;
t68 = t85 * t137;
t144 = t136 * (t61 + t68) - t137 * (-Icges(5,2) * t203 + t100 + t60);
t211 = Icges(5,4) * t117;
t157 = Icges(5,2) * t118 + t211;
t58 = Icges(5,6) * t137 + t136 * t157;
t59 = -Icges(5,6) * t136 + t137 * t157;
t87 = Icges(5,1) * t118 - t211;
t69 = t87 * t136;
t70 = t87 * t137;
t145 = t136 * (t59 - t70) - t137 * (t58 - t69);
t246 = -t145 * t117 + t144 * t118;
t221 = t85 + t158;
t222 = -t157 + t87;
t245 = (t117 * t221 - t118 * t222) * qJD(1);
t129 = t137 * rSges(5,3);
t62 = rSges(5,1) * t203 + rSges(5,2) * t202 + t129;
t127 = t136 * rSges(5,3);
t63 = t137 * t164 - t127;
t29 = (-t136 * t62 - t137 * t63) * qJD(4);
t244 = 0.2e1 * t29;
t133 = sin(pkin(6));
t201 = t133 * t136;
t217 = rSges(4,2) * cos(pkin(6));
t237 = -rSges(4,1) * t201 - t137 * rSges(4,3) - t136 * t217;
t98 = t137 * pkin(1) + t136 * qJ(2);
t243 = qJ(3) * t137 - t237 + t98;
t171 = -rSges(3,2) * t137 + t136 * rSges(3,3);
t242 = t171 + t98;
t161 = t117 * t61 + t118 * t59;
t239 = t161 * t137;
t109 = pkin(3) * t201;
t238 = t109 + t62 + t98;
t156 = Icges(5,5) * t117 + Icges(5,6) * t118;
t57 = -Icges(5,3) * t136 + t137 * t156;
t205 = qJD(1) * t57;
t28 = t117 * t59 - t118 * t61;
t35 = qJD(1) * t58 - qJD(4) * t68;
t37 = -qJD(4) * t70 + (t136 * t158 + t209) * qJD(1);
t236 = t28 * qJD(4) + t117 * t37 + t118 * t35 + t205;
t78 = t157 * qJD(4);
t79 = t158 * qJD(4);
t83 = Icges(5,5) * t118 - Icges(5,6) * t117;
t235 = qJD(1) * t83 + (t117 * t85 - t118 * t87) * qJD(4) + t117 * t79 + t118 * t78;
t162 = t117 * t58 - t118 * t60;
t56 = Icges(5,3) * t137 + t136 * t156;
t206 = qJD(1) * t56;
t189 = qJD(4) * t136;
t36 = qJD(1) * t59 + t189 * t85;
t38 = qJD(1) * t61 + qJD(4) * t69;
t234 = qJD(4) * t162 - t117 * t38 - t118 * t36 + t206;
t232 = t136 / 0.2e1;
t231 = -t137 / 0.2e1;
t229 = rSges(3,2) - pkin(1);
t228 = -rSges(5,3) - pkin(1);
t227 = -qJD(1) / 0.2e1;
t219 = rSges(5,1) * t118;
t216 = rSges(5,2) * t117;
t214 = rSges(3,3) * t137;
t65 = t136 * t83;
t213 = t137 * t83;
t187 = qJD(1) * qJD(2);
t191 = qJD(1) * t136;
t190 = qJD(1) * t137;
t114 = qJ(2) * t190;
t120 = qJD(2) * t136;
t194 = t114 + t120;
t212 = t136 * t187 + qJD(1) * (-pkin(1) * t191 + t194);
t208 = qJ(3) * t136;
t207 = qJ(3) * qJD(1) ^ 2;
t89 = -t216 + t219;
t204 = qJD(4) * t89;
t200 = t133 * t137;
t160 = t117 * t87 + t118 * t85;
t26 = t137 * t160 - t65;
t199 = t26 * qJD(1);
t198 = t156 * qJD(1);
t135 = -pkin(5) - qJ(3);
t197 = qJ(3) + t135;
t179 = t133 * t190;
t183 = t137 * t217;
t196 = rSges(4,1) * t179 + qJD(1) * t183;
t195 = pkin(3) * t179 + t135 * t191;
t193 = rSges(3,2) * t191 + rSges(3,3) * t190;
t192 = qJD(3) * t137 + t120;
t188 = qJD(4) * t137;
t186 = qJD(1) * qJD(3);
t185 = -rSges(4,3) - pkin(1) - qJ(3);
t16 = t137 * t56 + t58 * t202 + t60 * t203;
t17 = -t137 * t57 - t59 * t202 - t61 * t203;
t184 = t164 * t190 + t189 * t219;
t123 = t137 * qJ(2);
t95 = pkin(1) * t136 - t123;
t91 = qJD(1) * t95;
t182 = -t91 + t192;
t181 = qJD(4) * t216;
t180 = t114 + t192;
t176 = -t189 / 0.2e1;
t174 = -t188 / 0.2e1;
t165 = rSges(4,1) * t133 + t217;
t170 = -rSges(4,3) * t136 + t137 * t165 - t208;
t121 = qJD(2) * t137;
t169 = qJD(3) * t136 - t121;
t168 = pkin(3) * t200 + t136 * t197 - t208 + t63;
t163 = t117 * t60 + t118 * t58;
t113 = t137 * t187;
t155 = -0.2e1 * t136 * t186 - t137 * t207 + t113;
t72 = t89 * t137;
t153 = -t136 * t207 + 0.2e1 * t137 * t186 + t212;
t151 = (t136 * t17 + t137 * t16) * qJD(4);
t53 = t136 * t56;
t18 = -t163 * t137 + t53;
t19 = -t136 * t57 + t239;
t150 = (t136 * t19 + t137 * t18) * qJD(4);
t149 = pkin(3) * t133 + t164;
t143 = -qJD(1) * t161 - qJD(4) * t213 + t206;
t142 = qJD(1) * t163 + qJD(4) * t65 + t205;
t141 = t160 * qJD(1) - t156 * qJD(4);
t96 = rSges(3,2) * t136 + t214;
t81 = t164 * qJD(4);
t80 = qJD(1) * t98 - t121;
t74 = t89 * t189;
t71 = t89 * t136;
t52 = qJD(1) * t242 - t121;
t51 = t120 + (-t95 + t96) * qJD(1);
t44 = t113 + (-qJD(1) * t171 - t80) * qJD(1);
t43 = qJD(1) * t193 + t212;
t42 = (-rSges(5,3) * qJD(1) - t181) * t136 + t184;
t41 = -qJD(4) * t72 + (t136 * t164 + t129) * qJD(1);
t40 = qJD(1) * t243 + t169;
t39 = (t170 - t95) * qJD(1) + t192;
t25 = t136 * t160 + t213;
t24 = (qJD(1) * t237 - t80) * qJD(1) + t155;
t23 = qJD(1) * (-rSges(4,3) * t191 + t196) + t153;
t22 = t25 * qJD(1);
t21 = -t89 * t188 + (-t135 * t137 + t238) * qJD(1) + t169;
t20 = t74 + (t168 - t95) * qJD(1) + t192;
t12 = -t81 * t189 + (-qJD(1) * t109 - t41 - t80 + (qJD(1) * t197 + t204) * t137) * qJD(1) + t155;
t11 = t81 * t188 + (t42 + (qJ(3) * qJD(1) + t204) * t136 + t195) * qJD(1) + t153;
t10 = t161 * qJD(4) - t117 * t35 + t118 * t37;
t9 = -qJD(4) * t163 - t117 * t36 + t118 * t38;
t8 = -t136 * t235 + t141 * t137;
t7 = t141 * t136 + t137 * t235;
t6 = t150 - t199;
t5 = t22 + t151;
t1 = [(-qJD(4) * t160 + t117 * t78 - t118 * t79) * qJD(1) + (t22 + ((-t18 + t53 + t17) * t136 + (t19 - t239 + (-t163 + t57) * t136 + t16) * t137) * qJD(4)) * t176 + (t6 + t199 + (t136 ^ 2 * t57 + (-t53 + t17 + (t163 + t57) * t137) * t137) * qJD(4)) * t174 + (t12 * (t135 * t136 - t127 - t95) - t20 * t169 + t11 * t238 + t21 * (-t136 * t181 + t180 + t184 + t195) + (-t11 * t135 + t12 * t149 + t20 * t204) * t137 + (t20 * (t135 + t228) * t137 + (t20 * (-qJ(2) - t149) + t21 * t228) * t136) * qJD(1) - (qJD(1) * t168 + t182 - t20 + t74) * t21) * m(5) + (t24 * (rSges(4,1) * t200 + t123 + t183) + t39 * t121 + t23 * t243 + t40 * (t180 + t196) + (-t39 * qJD(3) + t24 * t185) * t136 + (t39 * t185 * t137 + (t39 * (-qJ(2) - t165) + t40 * t185) * t136) * qJD(1) - (qJD(1) * t170 + t182 - t39) * t40) * m(4) + (t44 * (t136 * t229 + t123 + t214) + t51 * t121 + t43 * t242 + t52 * (t193 + t194) + (t51 * t229 * t137 + (t51 * (-rSges(3,3) - qJ(2)) - t52 * pkin(1)) * t136) * qJD(1) - (qJD(1) * t96 + t120 - t51 - t91) * t52) * m(3) + (t10 + t7 + t5) * t189 / 0.2e1 + (qJD(1) * t28 + t8 + t9) * t188 / 0.2e1 + (t137 * t26 + (-t162 + t25) * t136) * qJD(4) * t227; 0.2e1 * (t11 * t231 + t12 * t232) * m(5) + 0.2e1 * (t23 * t231 + t232 * t24) * m(4) + 0.2e1 * (t43 * t231 + t232 * t44) * m(3); m(4) * (t136 * t23 + t137 * t24) + m(5) * (t11 * t136 + t12 * t137); qJD(1) * (t10 * t136 + t9 * t137 + (t136 * t162 + t28 * t137) * qJD(1)) / 0.2e1 + ((t65 * t188 - t198) * t137 + (-t245 + (-t137 * t213 - t246) * qJD(4)) * t136) * t174 + ((-t189 * t213 - t198) * t136 + (t245 + (t136 * t65 + t246) * qJD(4)) * t137) * t176 + ((-t117 * t222 - t118 * t221) * qJD(1) + (t117 * t144 + t118 * t145) * qJD(4)) * t227 + (qJD(1) * t7 + ((t142 * t136 + t137 * t234) * t137 + t136 * (t143 * t136 - t137 * t236) + (-t18 * t136 + t19 * t137) * qJD(1)) * t247) * t232 + (t8 * qJD(1) + (t136 * (t136 * t236 + t143 * t137) + t137 * (-t136 * t234 + t142 * t137) + (-t16 * t136 + t17 * t137) * qJD(1)) * t247) * t137 / 0.2e1 - (t5 + t151) * t191 / 0.2e1 + (t6 + t150) * t190 / 0.2e1 + ((t21 * t81 + (-qJD(1) * t62 + t41) * t244 + (qJD(1) * t20 - t11) * t89) * t137 + (-t20 * t81 + (qJD(1) * t63 - t42) * t244 + (qJD(1) * t21 + t12) * t89) * t136 - (t20 * t72 + t21 * t71) * qJD(1) - (t29 * (-t136 * t71 - t137 * t72) - (t136 * t20 - t137 * t21) * t164) * qJD(4)) * m(5);];
tauc = t1(:);
