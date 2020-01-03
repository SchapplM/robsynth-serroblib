% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP4
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:42
% DurationCPUTime: 2.32s
% Computational Cost: add. (5001->192), mult. (5388->263), div. (0->0), fcn. (4879->6), ass. (0->126)
t154 = qJ(1) + pkin(6);
t150 = sin(t154);
t151 = cos(t154);
t155 = sin(qJ(3));
t152 = Icges(5,5) * t155;
t157 = cos(qJ(3));
t200 = Icges(5,1) * t157;
t170 = t152 + t200;
t86 = Icges(5,4) * t150 + t170 * t151;
t199 = Icges(4,4) * t155;
t121 = Icges(4,1) * t157 - t199;
t88 = Icges(4,5) * t150 + t121 * t151;
t257 = t86 + t88;
t222 = rSges(5,1) + pkin(3);
t181 = t222 * t155;
t206 = rSges(5,3) + qJ(4);
t255 = -t206 * t157 + t181;
t256 = t255 * t150;
t148 = t150 ^ 2;
t149 = t151 ^ 2;
t183 = t148 + t149;
t236 = t206 * t155 + t222 * t157;
t192 = t151 * t157;
t193 = t151 * t155;
t133 = Icges(5,5) * t192;
t78 = Icges(5,6) * t150 + Icges(5,3) * t193 + t133;
t114 = Icges(4,5) * t157 - Icges(4,6) * t155;
t80 = Icges(4,3) * t150 + t114 * t151;
t115 = Icges(5,4) * t157 + Icges(5,6) * t155;
t82 = Icges(5,2) * t150 + t115 * t151;
t254 = t78 * t193 + t257 * t192 + (t80 + t82) * t150;
t252 = (-Icges(4,6) + Icges(5,6)) * t157 + (-Icges(5,4) - Icges(4,5)) * t155;
t153 = Icges(4,4) * t157;
t201 = Icges(4,1) * t155;
t171 = -t153 - t201;
t197 = Icges(4,2) * t155;
t84 = Icges(4,6) * t150 + (t153 - t197) * t151;
t251 = -Icges(5,1) * t193 + t171 * t151 + t133 + t78 - t84;
t198 = Icges(5,5) * t157;
t118 = Icges(5,1) * t155 - t198;
t113 = Icges(5,3) * t155 + t198;
t77 = -Icges(5,6) * t151 + t113 * t150;
t194 = t150 * t157;
t195 = t150 * t155;
t83 = Icges(4,4) * t194 - Icges(4,2) * t195 - Icges(4,6) * t151;
t250 = t77 - t83 + (-t118 + t171) * t150;
t116 = Icges(4,2) * t157 + t199;
t196 = Icges(5,3) * t157;
t167 = t196 - t152;
t249 = (-t116 - t167) * t151 + t257;
t134 = Icges(4,4) * t195;
t85 = -Icges(5,4) * t151 + t170 * t150;
t87 = Icges(4,1) * t194 - Icges(4,5) * t151 - t134;
t248 = -Icges(4,2) * t194 - t167 * t150 - t134 + t85 + t87;
t59 = t88 * t194;
t177 = t151 * t80 - t59;
t79 = Icges(4,5) * t194 - Icges(4,6) * t195 - Icges(4,3) * t151;
t218 = -t150 * t79 - t87 * t192;
t247 = -t83 * t193 - t84 * t195 - t177 - t218;
t246 = -t84 * t193 + t254;
t208 = t151 * (-Icges(5,2) * t151 + t115 * t150);
t245 = t208 + t254;
t243 = -t150 / 0.2e1;
t224 = t150 / 0.2e1;
t223 = -t151 / 0.2e1;
t241 = t252 * t150;
t240 = t252 * t151;
t124 = rSges(4,1) * t155 + rSges(4,2) * t157;
t104 = t124 * t150;
t106 = t124 * t151;
t178 = -sin(qJ(1)) * pkin(1) + t151 * pkin(5);
t233 = pkin(2) + t236;
t45 = t151 * rSges(5,2) - t150 * t233 + t178;
t221 = cos(qJ(1)) * pkin(1);
t46 = t221 + (rSges(5,2) + pkin(5)) * t150 + t233 * t151;
t210 = rSges(4,1) * t157;
t179 = pkin(2) + t210;
t184 = rSges(4,2) * t195 + t151 * rSges(4,3);
t50 = -t179 * t150 + t178 + t184;
t136 = rSges(4,2) * t193;
t51 = t221 - t136 + t179 * t151 + (rSges(4,3) + pkin(5)) * t150;
t56 = -t151 * t181 + t206 * t192;
t9 = (t121 / 0.2e1 - t116 / 0.2e1 + t152 + t200 / 0.2e1 - t196 / 0.2e1) * t155 + (t153 + t201 / 0.2e1 - t197 / 0.2e1 + t118 / 0.2e1 - t113 / 0.2e1) * t157 + m(5) * (t256 * t45 + t46 * t56) + m(4) * (t104 * t50 - t106 * t51);
t239 = t9 * qJD(1);
t238 = (t155 * t77 + t157 * t85) * t150;
t235 = -t249 * t155 + t251 * t157;
t234 = t248 * t155 - t250 * t157;
t232 = m(5) / 0.2e1;
t220 = t45 * t192 + t46 * t194;
t229 = m(5) * ((t150 * t56 + t151 * t256) * t155 + t220);
t72 = t255 * t151;
t228 = m(5) * (-t193 * t256 + t72 * t195 + t220);
t109 = t183 * t155;
t225 = t109 / 0.2e1;
t219 = -t72 * t192 - t194 * t256;
t212 = m(5) * qJD(3);
t211 = m(5) * qJD(4);
t207 = t155 * t83;
t191 = t155 * t157;
t76 = (t225 - t155 / 0.2e1) * m(5);
t188 = t76 * qJD(2);
t187 = t183 * t191;
t25 = t46 * t193 - t45 * t195;
t182 = m(5) * t25 * qJD(1);
t180 = t114 / 0.2e1 + t115 / 0.2e1;
t176 = t155 * t84 - t79;
t175 = t151 * t82 - t86 * t194 - t78 * t195;
t30 = -t208 + t238;
t166 = (t30 - t238 + t245) * t243 + t246 * t224 + (-t59 + (t80 + t207) * t151 + t218 + t247) * t223;
t165 = t175 * t243 + (t176 * t151 - t245 + t246) * t151 / 0.2e1 + (-(-t157 * t87 + t207) * t150 - t151 * t79 + t30) * t223 + (t176 * t150 + t85 * t192 + t77 * t193 + t175 + t177 + t247) * t224;
t128 = -rSges(4,2) * t155 + t210;
t75 = m(5) * t225 + t155 * t232;
t74 = t187 - t191;
t73 = t236 * t151;
t71 = t236 * t150;
t49 = -t104 * t150 - t106 * t151;
t38 = -t148 * t255 + t56 * t151;
t28 = t183 * t236;
t19 = t109 * t28 + t219;
t13 = t228 / 0.2e1;
t11 = t229 / 0.2e1;
t4 = t13 - t229 / 0.2e1;
t3 = t13 + t11;
t2 = t11 - t228 / 0.2e1;
t1 = t165 * t150 + t166 * t151;
t5 = [t9 * qJD(3) + t25 * t211, 0, t239 + t3 * qJD(4) + (m(5) * (-t45 * t73 - t46 * t71 + (-t56 - t72) * t256) + (m(4) * (-t104 * t124 - t128 * t50) + t180 * t151 - t166) * t151 + (m(4) * (t106 * t124 - t128 * t51) + t180 * t150 - t165) * t150 + (t251 * t155 + t249 * t157) * t224 + (t250 * t155 + t248 * t157) * t223) * qJD(3), t3 * qJD(3) + t182; 0, 0, 0.2e1 * (m(4) * t49 / 0.2e1 + t38 * t232) * qJD(3) + t75 * qJD(4), t75 * qJD(3); t1 * qJD(3) + t4 * qJD(4) - t239, qJD(4) * t76, t1 * qJD(1) + (m(4) * (t183 * t128 * t124 + (t150 * (rSges(4,1) * t194 - t184) + t151 * (rSges(4,1) * t192 + t150 * rSges(4,3) - t136)) * t49) + m(5) * (t256 * t71 + t28 * t38 + t72 * t73) + (t240 * t148 + (t234 * t151 + (t235 - t241) * t150) * t151) * t224 + (t241 * t149 + (t235 * t150 + (t234 - t240) * t151) * t150) * t223) * qJD(3) + t19 * t211, t4 * qJD(1) + t188 + t19 * t212 + (-t109 * t157 + t187 - t74) * t211; t2 * qJD(3) - t182, -t76 * qJD(3), t2 * qJD(1) - t188 + (-t157 * t38 + (-t150 * t71 - t151 * t73 + t28) * t155 - t19 + t219) * t212 + t74 * t211, t74 * t212;];
Cq = t5;
