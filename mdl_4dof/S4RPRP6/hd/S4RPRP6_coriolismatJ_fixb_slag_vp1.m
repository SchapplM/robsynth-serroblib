% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:04
% DurationCPUTime: 2.66s
% Computational Cost: add. (2076->170), mult. (4806->251), div. (0->0), fcn. (4264->4), ass. (0->109)
t243 = Icges(4,5) + Icges(5,5);
t242 = Icges(4,6) + Icges(5,6);
t152 = cos(qJ(3));
t150 = sin(qJ(3));
t184 = Icges(5,4) * t150;
t164 = Icges(5,2) * t152 + t184;
t186 = Icges(4,4) * t150;
t165 = Icges(4,2) * t152 + t186;
t244 = t164 + t165;
t151 = sin(qJ(1));
t181 = t151 * t152;
t241 = (-Icges(4,4) - Icges(5,4)) * t181;
t153 = cos(qJ(1));
t237 = t244 * t151 + t242 * t153;
t182 = t150 * t151;
t235 = -t241 + (Icges(4,1) + Icges(5,1)) * t182 + t243 * t153;
t195 = rSges(5,1) + pkin(3);
t233 = t152 * rSges(5,2) + t150 * t195;
t82 = t233 * t153;
t240 = t151 * rSges(5,3) - t82;
t239 = Icges(5,3) + Icges(4,3);
t160 = Icges(5,5) * t150 + Icges(5,6) * t152;
t162 = Icges(4,5) * t150 + Icges(4,6) * t152;
t238 = t160 + t162;
t236 = t242 * t151 - t244 * t153;
t183 = Icges(5,4) * t152;
t166 = Icges(5,1) * t150 + t183;
t185 = Icges(4,4) * t152;
t167 = Icges(4,1) * t150 + t185;
t234 = (-t166 - t167) * t153 + t243 * t151;
t225 = -t239 * t151 + t238 * t153;
t224 = t235 * t150 + t237 * t152;
t232 = -t242 * t150 + t243 * t152;
t231 = (t234 * t150 + t236 * t152) * t153;
t148 = t151 ^ 2;
t230 = m(5) / 0.2e1;
t138 = t153 * qJ(2);
t192 = -qJ(4) - pkin(5);
t53 = t138 + (-pkin(1) + t192) * t151 - t240;
t137 = t153 * t192;
t54 = -t137 + (rSges(5,3) + pkin(1)) * t153 + (qJ(2) + t233) * t151;
t229 = m(5) * (-t151 * t53 + t153 * t54);
t228 = t235 * t182 + t237 * t181 + (t238 * t151 + t239 * t153) * t153;
t227 = -t225 * t153 + t236 * t181 + t234 * t182;
t226 = -t225 * t151 - t231;
t223 = t224 * t153;
t222 = t148 / 0.2e1;
t197 = t151 / 0.2e1;
t196 = t153 / 0.2e1;
t218 = qJD(1) * t229;
t80 = t233 * t151;
t215 = t232 * t151;
t214 = t232 * t153;
t124 = Icges(5,1) * t152 - t184;
t126 = Icges(4,1) * t152 - t186;
t213 = t124 + t126;
t131 = rSges(4,1) * t152 - rSges(4,2) * t150;
t113 = t131 * t151;
t114 = t131 * t153;
t136 = pkin(3) * t181;
t174 = rSges(5,1) * t181 - rSges(5,2) * t182;
t77 = t136 + t174;
t188 = rSges(5,2) * t150;
t78 = (t152 * t195 - t188) * t153;
t191 = (t151 * t78 - t153 * t77) * t230 + m(4) * (-t113 * t153 + t151 * t114) / 0.2e1;
t120 = -Icges(5,2) * t150 + t183;
t122 = -Icges(4,2) * t150 + t185;
t212 = (t126 / 0.2e1 - t165 / 0.2e1 + t124 / 0.2e1 - t164 / 0.2e1) * t150 - (-t167 / 0.2e1 - t122 / 0.2e1 - t166 / 0.2e1 - t120 / 0.2e1) * t152;
t149 = t153 ^ 2;
t211 = 0.4e1 * qJD(1);
t210 = pkin(1) + pkin(5);
t209 = m(3) * (t153 * (rSges(3,3) * t153 + t138) + (rSges(3,3) + qJ(2)) * t148);
t173 = rSges(4,1) * t150 + rSges(4,2) * t152;
t154 = -t151 * rSges(4,3) + t153 * t173;
t57 = -t151 * t210 + t138 + t154;
t58 = (rSges(4,3) + t210) * t153 + (qJ(2) + t173) * t151;
t208 = m(4) * (t113 * t58 + t114 * t57);
t207 = m(4) * (t58 * t151 + t153 * t57);
t205 = m(5) * (t53 * t78 + t54 * t77);
t204 = m(5) * (t151 * t54 + t153 * t53);
t202 = m(5) * (t151 * t77 + t153 * t78);
t130 = rSges(5,1) * t152 - t188;
t81 = t130 * t151 + t136;
t193 = t152 * pkin(3);
t83 = (t130 + t193) * t153;
t201 = m(5) * (-t83 * t151 + t153 * t81);
t199 = m(5) * (-t151 * t81 - t153 * t83);
t189 = m(5) * qJD(3);
t180 = t151 * t153;
t67 = 0.2e1 * (t222 + t149 / 0.2e1) * m(5);
t179 = t67 * qJD(1);
t178 = t148 + t149;
t177 = -t162 / 0.2e1 - t160 / 0.2e1;
t175 = t178 * t173;
t159 = -t227 * t151 / 0.2e1 - t228 * t153 / 0.2e1 + (t223 + t227) * t197 + ((-t224 + t225) * t151 + t226 + t228 + t231) * t196;
t158 = t225 * t222 + t226 * t197 + ((t224 + t225) * t153 - t223 + t227) * t196;
t37 = t199 / 0.2e1;
t35 = t201 / 0.2e1;
t34 = t202 / 0.2e1;
t14 = t37 - t202 / 0.2e1;
t13 = t37 + t34;
t12 = t34 - t199 / 0.2e1;
t11 = -t201 / 0.2e1 + t191;
t10 = t35 + t191;
t9 = t35 - t191;
t8 = t204 + t207 + t209;
t6 = t205 + t208 - t212;
t1 = t151 * t159 + t153 * t158;
t2 = [t8 * qJD(2) + t6 * qJD(3) + qJD(4) * t229, qJD(1) * t8 + qJD(3) * t10, t6 * qJD(1) + t10 * qJD(2) + (-t53 * t80 + t54 * t82 - t77 * t83 + t78 * t81) * t189 + t13 * qJD(4) + ((m(4) * (-t113 * t131 + t173 * t58) + t177 * t153 - t158) * t153 + (m(4) * (t114 * t131 - t173 * t57) + t177 * t151 - t159) * t151 + (((t120 + t122) * t153 - t234) * t197 + ((Icges(4,2) + Icges(5,2)) * t182 - t235 + t241) * t196) * t150 + ((-t213 * t153 - t236) * t197 + (t213 * t151 - t237) * t196) * t152) * qJD(3), t13 * qJD(3) + t218; t11 * qJD(3) - t67 * qJD(4) + (-t209 / 0.4e1 - t207 / 0.4e1 - t204 / 0.4e1) * t211, 0, t11 * qJD(1) + 0.2e1 * ((-t80 * t151 - t153 * t82) * t230 - m(4) * t175 / 0.2e1) * qJD(3), -t179; t9 * qJD(2) + t1 * qJD(3) + t14 * qJD(4) + (-t208 / 0.4e1 - t205 / 0.4e1) * t211 + t212 * qJD(1), t9 * qJD(1), t1 * qJD(1) + ((-t214 * t148 + t215 * t180) * t197 + (t215 * t149 - t214 * t180) * t196 + m(4) * ((-t153 * t154 + (-rSges(4,3) * t153 - t151 * t173) * t151) * (-t151 * t113 - t114 * t153) - t131 * t175) + m(5) * (-t81 * t80 - t83 * t82 + ((t137 - t80) * t151 + ((-rSges(5,3) - t192) * t151 + t240) * t153) * (-t130 * t149 - t151 * t174 - t178 * t193))) * qJD(3), t14 * qJD(1); t67 * qJD(2) + t12 * qJD(3) - t218, t179, t12 * qJD(1) + (t151 * t82 - t153 * t80) * t189, 0;];
Cq = t2;
