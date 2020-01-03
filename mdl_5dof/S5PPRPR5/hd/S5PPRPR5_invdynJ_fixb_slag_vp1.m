% Calculate vector of inverse dynamics joint torques for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 3.56s
% Computational Cost: add. (3311->341), mult. (7626->449), div. (0->0), fcn. (8074->6), ass. (0->171)
t147 = sin(qJ(5));
t145 = sin(pkin(7));
t146 = cos(pkin(7));
t238 = sin(qJ(3));
t239 = cos(qJ(3));
t124 = -t145 * t239 + t146 * t238;
t152 = t145 * t238 + t146 * t239;
t144 = Icges(6,4) * t147;
t148 = cos(qJ(5));
t162 = Icges(6,2) * t148 + t144;
t56 = Icges(6,6) * t152 + t124 * t162;
t256 = t147 * t56;
t255 = t148 * t56;
t130 = Icges(6,5) * t147 + Icges(6,6) * t148;
t53 = Icges(6,3) * t152 + t124 * t130;
t229 = t152 * t53;
t114 = t124 * qJ(4);
t181 = -pkin(3) * t152 - t114;
t118 = t124 * rSges(5,3);
t208 = rSges(5,2) * t152 - t118;
t231 = t181 + t208;
t254 = -t152 * pkin(6) + t181;
t204 = qJD(4) * t124;
t85 = -t124 * pkin(3) + qJ(4) * t152;
t253 = qJD(3) * t85 + t204;
t161 = Icges(6,5) * t148 - Icges(6,6) * t147;
t67 = t161 * t124;
t68 = t161 * t152;
t83 = -rSges(4,1) * t152 + rSges(4,2) * t124;
t252 = t83 * qJD(3);
t107 = t124 * qJD(3);
t198 = qJD(5) * t148;
t157 = -t107 * t147 + t152 * t198;
t223 = Icges(6,4) * t148;
t133 = Icges(6,2) * t147 - t223;
t134 = Icges(6,1) * t147 + t223;
t206 = t133 - t134;
t135 = -Icges(6,1) * t148 + t144;
t207 = t135 + t162;
t251 = (t147 * t206 - t148 * t207) * qJD(3);
t212 = t152 * t148;
t213 = t152 * t147;
t64 = rSges(6,1) * t213 + rSges(6,2) * t212 - t124 * rSges(6,3);
t187 = pkin(6) * t124 - t64;
t136 = t147 * rSges(6,1) + rSges(6,2) * t148;
t250 = t124 * t136;
t84 = t124 * rSges(5,2) + rSges(5,3) * t152;
t249 = qJD(3) * t84 + t253;
t202 = qJD(5) * t124;
t173 = rSges(6,1) * t148 - t147 * rSges(6,2);
t200 = qJD(5) * t173;
t108 = t152 * qJD(3);
t248 = -t108 * pkin(3) - qJ(4) * t107;
t211 = t124 * t147;
t222 = Icges(6,5) * t152;
t210 = t124 * t148;
t97 = Icges(6,4) * t210;
t59 = Icges(6,1) * t211 + t222 + t97;
t23 = -t148 * t59 + t256;
t199 = qJD(5) * t147;
t215 = t108 * t148;
t154 = -t124 * t199 + t215;
t216 = t108 * t147;
t155 = t124 * t198 + t216;
t28 = Icges(6,4) * t155 + Icges(6,2) * t154 - Icges(6,6) * t107;
t30 = Icges(6,1) * t155 + Icges(6,4) * t154 - Icges(6,5) * t107;
t247 = qJD(5) * t23 - t147 * t30 - t148 * t28;
t58 = -Icges(6,6) * t124 + t152 * t162;
t98 = Icges(6,4) * t212;
t61 = Icges(6,1) * t213 - Icges(6,5) * t124 + t98;
t24 = t147 * t58 - t148 * t61;
t189 = t152 * t199;
t217 = t107 * t148;
t156 = -t189 - t217;
t27 = Icges(6,4) * t157 + Icges(6,2) * t156 - Icges(6,6) * t108;
t29 = Icges(6,1) * t157 + Icges(6,4) * t156 - Icges(6,5) * t108;
t246 = qJD(5) * t24 - t147 * t29 - t148 * t27;
t159 = t148 * t133 + t147 * t135;
t38 = t124 * t159 - t68;
t126 = t162 * qJD(5);
t127 = t134 * qJD(5);
t158 = t147 * t133 - t135 * t148;
t245 = qJD(5) * t158 - t126 * t148 - t127 * t147;
t244 = -qJD(3) * t187 + t124 * t200 + t253;
t240 = -rSges(5,2) + pkin(3);
t235 = t124 * t135 + t56;
t234 = t135 * t152 + t58;
t233 = -Icges(6,2) * t211 + t59 + t97;
t232 = -Icges(6,2) * t213 + t61 + t98;
t230 = rSges(5,3) * t108;
t115 = t152 * rSges(6,3);
t16 = -t124 * t53 + t56 * t212 + t59 * t213;
t228 = t16 * t152;
t227 = t108 * rSges(5,2) - t107 * rSges(5,3);
t226 = t108 * qJ(4) + t204;
t113 = t152 * qJD(4);
t225 = qJD(3) * t181 + t113;
t39 = t152 * t159 + t67;
t209 = t39 * qJD(3);
t205 = qJD(2) * t146;
t203 = qJD(5) * t152;
t128 = t136 * qJD(5);
t201 = qJD(5) * t128;
t197 = t130 * qJD(3);
t196 = m(3) + m(4) + m(5);
t195 = qJDD(2) * t146;
t194 = rSges(6,3) + pkin(3) + pkin(6);
t55 = -Icges(6,3) * t124 + t130 * t152;
t15 = t152 * t55 + t58 * t210 + t61 * t211;
t193 = t124 * t55 - t58 * t212 - t61 * t213;
t192 = -m(6) - t196;
t63 = -rSges(6,1) * t211 - rSges(6,2) * t210 - t115;
t186 = t203 / 0.2e1;
t185 = -t203 / 0.2e1;
t184 = t202 / 0.2e1;
t182 = rSges(6,1) * t157 - rSges(6,2) * t217 - t108 * rSges(6,3);
t180 = t113 - t205;
t62 = t115 + t250;
t179 = -t62 + t254;
t142 = qJDD(2) * t145;
t153 = t113 + t248;
t178 = qJD(3) * t153 + qJD(4) * t108 + qJDD(3) * t85 + qJDD(4) * t124 + t142;
t164 = t147 * t61 + t148 * t58;
t165 = t147 * t59 + t255;
t25 = Icges(6,5) * t157 + Icges(6,6) * t156 - Icges(6,3) * t108;
t26 = Icges(6,5) * t155 + Icges(6,6) * t154 - Icges(6,3) * t107;
t176 = (-t165 * t107 - t108 * t53 - t124 * t26 - t152 * t247) * t152 - t124 * (-t164 * t107 - t108 * t55 - t124 * t25 - t152 * t246);
t175 = t152 * (-t107 * t53 + t165 * t108 - t247 * t124 + t152 * t26) - t124 * (-t107 * t55 + t164 * t108 - t246 * t124 + t152 * t25);
t174 = -qJD(4) * t107 + qJDD(4) * t152 - t195;
t14 = t124 * t165 + t229;
t172 = -t124 * t15 + t14 * t152;
t171 = t124 * t193 + t228;
t143 = qJD(2) * t145;
t20 = t143 + t244;
t21 = qJD(3) * t179 + t152 * t200 + t180;
t170 = -t124 * t20 - t152 * t21;
t31 = -rSges(6,2) * t189 + t182;
t32 = rSges(6,1) * t155 + rSges(6,2) * t154 - t107 * rSges(6,3);
t169 = t124 * t32 + t152 * t31;
t168 = t124 * t62 + t152 * t64;
t151 = t147 * t235 - t148 * t233;
t150 = -t147 * t234 + t148 * t232;
t77 = -qJD(5) * t108 - qJDD(5) * t124;
t8 = -t124 * t201 + qJD(3) * t31 + qJDD(3) * t64 - t173 * t77 + (-qJD(3) * t108 - qJDD(3) * t124) * pkin(6) + t178;
t52 = -t107 * pkin(3) + t226;
t78 = -qJD(5) * t107 + qJDD(5) * t152;
t9 = -t152 * t201 + t173 * t78 + t179 * qJDD(3) + (t107 * pkin(6) - t32 - t52) * qJD(3) + t174;
t149 = t107 * t21 - t108 * t20 - t124 * t8 - t152 * t9;
t125 = t130 * qJD(5);
t86 = -rSges(4,1) * t124 - rSges(4,2) * t152;
t76 = t173 * t152;
t75 = t173 * t124;
t74 = -rSges(4,1) * t108 + rSges(4,2) * t107;
t73 = -rSges(4,1) * t107 - rSges(4,2) * t108;
t66 = -t205 + t252;
t60 = -t124 * t134 - t222;
t37 = -qJD(3) * t73 + qJDD(3) * t83 - t195;
t36 = qJD(3) * t74 + qJDD(3) * t86 + t142;
t35 = qJD(3) * t231 + t180;
t34 = t143 + t249;
t22 = qJD(5) * t168 + qJD(1);
t19 = t231 * qJDD(3) + (-rSges(5,2) * t107 - t230 - t52) * qJD(3) + t174;
t18 = qJD(3) * t227 + qJDD(3) * t84 + t178;
t13 = t107 * t161 + t159 * t108 - t245 * t124 + t125 * t152;
t12 = -t159 * t107 + t108 * t161 - t124 * t125 - t152 * t245;
t11 = qJD(5) * t164 + t147 * t27 - t148 * t29;
t10 = qJD(5) * t165 + t147 * t28 - t148 * t30;
t7 = qJD(5) * t169 - t62 * t77 + t64 * t78 + qJDD(1);
t6 = qJD(5) * t171 - t209;
t5 = -qJD(3) * t38 + qJD(5) * t172;
t1 = [m(6) * t7 + (m(2) + t196) * qJDD(1) + (-m(2) + t192) * g(3); t192 * (g(1) * t145 - g(2) * t146) + m(4) * (t145 * t36 - t146 * t37) + m(5) * (t145 * t18 - t146 * t19) + m(6) * (t145 * t8 - t146 * t9) + m(3) * (t145 ^ 2 + t146 ^ 2) * qJDD(2); t6 * t186 + (t159 * qJD(5) + t147 * t126 - t127 * t148) * qJD(3) - (t39 + t24) * t77 / 0.2e1 - (t38 + t23) * t78 / 0.2e1 + (-t209 + (t228 + (t229 - t14 + (-t147 * t60 + t255) * t124 + t193) * t124) * qJD(5) + t10 + t13) * t185 + (t158 + Icges(4,3) + Icges(5,1)) * qJDD(3) + (-g(2) * (t63 + t254) + (-t114 - t250) * t9 + (-g(1) + t8) * (t85 - t187) + (-rSges(6,1) * t216 - rSges(6,2) * t215 + t194 * t107 - t173 * t202 - t226 + t244) * t21 + (-pkin(6) * t108 - qJD(3) * t63 + t182 - t225 + t248) * t20 - (t9 * t194 - t22 * (-t62 - t63) * qJD(5) + (rSges(6,2) * t199 - pkin(6) * qJD(3) - qJD(4) + t200) * t20) * t152) * m(6) + (t19 * (-t152 * t240 - t114 - t118) - g(2) * t231 + (t18 - g(1)) * (t85 + t84) + (t107 * t240 - t226 - t230 + t249) * t35 + (-qJD(3) * t208 + t153 - t225 + t227) * t34) * m(5) + (-t66 * t73 + (t74 - t252) * (qJD(3) * t86 + t143) + (t66 * qJD(3) - g(1) + t36) * t86 + (-g(2) + t37) * t83) * m(4) + ((-t229 * t152 + (-t16 + (t164 - t53) * t124 - (-t55 + (t59 + t60) * t147) * t152) * t124) * qJD(5) + t11 + t12 + t5 + (t148 * t60 - t23 + t256 + t38) * qJD(3)) * t184; (-m(5) - m(6)) * (g(1) * t124 + g(2) * t152) + m(5) * (-t107 * t35 + t108 * t34 + t124 * t18 + t152 * t19) - m(6) * t149 + 0.2e1 * (-m(5) * (-t124 * t35 + t152 * t34) / 0.2e1 - m(6) * (-t124 * t21 + t152 * t20) / 0.2e1) * qJD(3); -t107 * t5 / 0.2e1 + t152 * (-t13 * qJD(3) + t175 * qJD(5) - t38 * qJDD(3) + t14 * t78 + t15 * t77) / 0.2e1 + t78 * t172 / 0.2e1 + (-t107 * t14 - t108 * t15 + t175) * t186 - t108 * t6 / 0.2e1 - t124 * (-t12 * qJD(3) + t176 * qJD(5) - t39 * qJDD(3) + t16 * t78 - t193 * t77) / 0.2e1 + t77 * t171 / 0.2e1 - (-t107 * t16 + t108 * t193 + t176) * t202 / 0.2e1 - qJDD(3) * (-t124 * t24 + t152 * t23) / 0.2e1 - qJD(3) * (t10 * t152 - t23 * t107 - t24 * t108 - t11 * t124) / 0.2e1 + (-(-t67 * t203 + t197) * t152 + (t251 + (-t150 * t124 - (t68 + t151) * t152) * qJD(5)) * t124) * t185 + ((t68 * t202 + t197) * t124 - (-t251 + (t151 * t152 + (t67 + t150) * t124) * qJD(5)) * t152) * t184 + qJD(3) * (-(t207 * t147 + t206 * t148) * qJD(3) + ((-t124 * t234 + t152 * t235) * t148 + (-t124 * t232 + t152 * t233) * t147) * qJD(5)) / 0.2e1 + (t7 * t168 + t22 * (-t107 * t64 + t108 * t62 + t169) + t170 * t128 - t149 * t173 - (t20 * t76 - t21 * t75) * qJD(3) - (t22 * (t124 * t75 + t152 * t76) + t170 * t136) * qJD(5) - g(1) * t75 - g(2) * t76 - g(3) * t136) * m(6);];
tau = t1;
