% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP4
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:49
% DurationCPUTime: 2.30s
% Computational Cost: add. (4929->189), mult. (5308->261), div. (0->0), fcn. (4805->4), ass. (0->125)
t152 = pkin(6) + qJ(2);
t148 = sin(t152);
t149 = cos(t152);
t153 = sin(qJ(3));
t150 = Icges(5,5) * t153;
t154 = cos(qJ(3));
t195 = Icges(5,1) * t154;
t166 = t150 + t195;
t86 = Icges(5,4) * t148 + t166 * t149;
t194 = Icges(4,4) * t153;
t121 = Icges(4,1) * t154 - t194;
t88 = Icges(4,5) * t148 + t121 * t149;
t251 = t86 + t88;
t216 = rSges(5,1) + pkin(3);
t176 = t216 * t153;
t201 = rSges(5,3) + qJ(4);
t249 = -t201 * t154 + t176;
t250 = t249 * t148;
t146 = t148 ^ 2;
t147 = t149 ^ 2;
t178 = t146 + t147;
t230 = t201 * t153 + t216 * t154;
t187 = t149 * t154;
t188 = t149 * t153;
t131 = Icges(5,5) * t187;
t78 = Icges(5,6) * t148 + Icges(5,3) * t188 + t131;
t114 = Icges(4,5) * t154 - Icges(4,6) * t153;
t80 = Icges(4,3) * t148 + t149 * t114;
t115 = Icges(5,4) * t154 + Icges(5,6) * t153;
t82 = Icges(5,2) * t148 + t149 * t115;
t248 = t78 * t188 + t251 * t187 + (t80 + t82) * t148;
t246 = (-Icges(4,6) + Icges(5,6)) * t154 + (-Icges(5,4) - Icges(4,5)) * t153;
t151 = Icges(4,4) * t154;
t196 = Icges(4,1) * t153;
t167 = -t151 - t196;
t192 = Icges(4,2) * t153;
t84 = Icges(4,6) * t148 + (t151 - t192) * t149;
t245 = -Icges(5,1) * t188 + t167 * t149 + t131 + t78 - t84;
t193 = Icges(5,5) * t154;
t118 = Icges(5,1) * t153 - t193;
t113 = Icges(5,3) * t153 + t193;
t77 = -Icges(5,6) * t149 + t113 * t148;
t189 = t148 * t154;
t190 = t148 * t153;
t83 = Icges(4,4) * t189 - Icges(4,2) * t190 - Icges(4,6) * t149;
t244 = t77 - t83 + (-t118 + t167) * t148;
t116 = Icges(4,2) * t154 + t194;
t191 = Icges(5,3) * t154;
t163 = t191 - t150;
t243 = (-t116 - t163) * t149 + t251;
t132 = Icges(4,4) * t190;
t85 = -Icges(5,4) * t149 + t166 * t148;
t87 = Icges(4,1) * t189 - Icges(4,5) * t149 - t132;
t242 = -Icges(4,2) * t189 - t163 * t148 - t132 + t85 + t87;
t59 = t88 * t189;
t173 = t149 * t80 - t59;
t79 = Icges(4,5) * t189 - Icges(4,6) * t190 - Icges(4,3) * t149;
t213 = -t148 * t79 - t87 * t187;
t241 = -t83 * t188 - t84 * t190 - t173 - t213;
t240 = -t84 * t188 + t248;
t203 = t149 * (-Icges(5,2) * t149 + t148 * t115);
t239 = t203 + t248;
t237 = -t148 / 0.2e1;
t218 = t148 / 0.2e1;
t217 = -t149 / 0.2e1;
t235 = t246 * t148;
t234 = t246 * t149;
t124 = t153 * rSges(4,1) + rSges(4,2) * t154;
t104 = t124 * t148;
t106 = t124 * t149;
t144 = t149 * pkin(5);
t227 = pkin(2) + t230;
t45 = t149 * rSges(5,2) - t148 * t227 + t144;
t46 = (rSges(5,2) + pkin(5)) * t148 + t227 * t149;
t205 = rSges(4,1) * t154;
t174 = pkin(2) + t205;
t179 = rSges(4,2) * t190 + t149 * rSges(4,3);
t53 = -t174 * t148 + t144 + t179;
t134 = rSges(4,2) * t188;
t54 = -t134 + t174 * t149 + (rSges(4,3) + pkin(5)) * t148;
t56 = -t149 * t176 + t201 * t187;
t9 = (t121 / 0.2e1 - t116 / 0.2e1 + t150 + t195 / 0.2e1 - t191 / 0.2e1) * t153 + (t151 + t196 / 0.2e1 - t192 / 0.2e1 + t118 / 0.2e1 - t113 / 0.2e1) * t154 + m(5) * (t250 * t45 + t46 * t56) + m(4) * (t104 * t53 - t106 * t54);
t233 = t9 * qJD(2);
t232 = (t153 * t77 + t154 * t85) * t148;
t229 = -t243 * t153 + t245 * t154;
t228 = t242 * t153 - t244 * t154;
t226 = m(5) / 0.2e1;
t215 = t45 * t187 + t46 * t189;
t223 = m(5) * ((t148 * t56 + t149 * t250) * t153 + t215);
t72 = t249 * t149;
t222 = m(5) * (-t188 * t250 + t72 * t190 + t215);
t107 = t178 * t153;
t219 = t107 / 0.2e1;
t214 = -t72 * t187 - t189 * t250;
t207 = m(5) * qJD(3);
t206 = m(5) * qJD(4);
t202 = t153 * t83;
t186 = t153 * t154;
t76 = (t219 - t153 / 0.2e1) * m(5);
t183 = t76 * qJD(1);
t182 = t178 * t186;
t25 = t46 * t188 - t45 * t190;
t177 = m(5) * t25 * qJD(2);
t175 = t114 / 0.2e1 + t115 / 0.2e1;
t172 = t153 * t84 - t79;
t171 = t149 * t82 - t86 * t189 - t78 * t190;
t30 = -t203 + t232;
t162 = (t30 - t232 + t239) * t237 + t240 * t218 + (-t59 + (t80 + t202) * t149 + t213 + t241) * t217;
t161 = t171 * t237 + (t172 * t149 - t239 + t240) * t149 / 0.2e1 + (-(-t154 * t87 + t202) * t148 - t149 * t79 + t30) * t217 + (t172 * t148 + t85 * t187 + t77 * t188 + t171 + t173 + t241) * t218;
t127 = -t153 * rSges(4,2) + t205;
t75 = m(5) * t219 + t153 * t226;
t74 = t182 - t186;
t73 = t230 * t149;
t71 = t230 * t148;
t49 = -t104 * t148 - t106 * t149;
t38 = -t146 * t249 + t56 * t149;
t28 = t178 * t230;
t19 = t107 * t28 + t214;
t17 = t222 / 0.2e1;
t11 = t223 / 0.2e1;
t4 = t17 - t223 / 0.2e1;
t3 = t17 + t11;
t2 = t11 - t222 / 0.2e1;
t1 = t161 * t148 + t162 * t149;
t5 = [0, 0, 0.2e1 * (m(4) * t49 / 0.2e1 + t38 * t226) * qJD(3) + t75 * qJD(4), t75 * qJD(3); 0, t9 * qJD(3) + t25 * t206, t233 + t3 * qJD(4) + (m(5) * (-t45 * t73 - t46 * t71 + (-t56 - t72) * t250) + (m(4) * (-t104 * t124 - t127 * t53) + t175 * t149 - t162) * t149 + (m(4) * (t106 * t124 - t127 * t54) + t175 * t148 - t161) * t148 + (t245 * t153 + t243 * t154) * t218 + (t244 * t153 + t242 * t154) * t217) * qJD(3), t3 * qJD(3) + t177; qJD(4) * t76, t1 * qJD(3) + t4 * qJD(4) - t233, t1 * qJD(2) + (m(4) * (t178 * t127 * t124 + (t148 * (rSges(4,1) * t189 - t179) + t149 * (rSges(4,1) * t187 + t148 * rSges(4,3) - t134)) * t49) + m(5) * (t250 * t71 + t28 * t38 + t72 * t73) + (t234 * t146 + (t228 * t149 + (t229 - t235) * t148) * t149) * t218 + (t235 * t147 + (t229 * t148 + (t228 - t234) * t149) * t148) * t217) * qJD(3) + t19 * t206, t183 + t4 * qJD(2) + t19 * t207 + (-t107 * t154 + t182 - t74) * t206; -t76 * qJD(3), t2 * qJD(3) - t177, -t183 + t2 * qJD(2) + (-t154 * t38 + (-t148 * t71 - t149 * t73 + t28) * t153 - t19 + t214) * t207 + t74 * t206, t74 * t207;];
Cq = t5;
