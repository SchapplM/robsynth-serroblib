% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:30
% DurationCPUTime: 1.42s
% Computational Cost: add. (7050->172), mult. (6714->241), div. (0->0), fcn. (5798->6), ass. (0->119)
t262 = m(5) / 0.2e1;
t263 = m(4) / 0.2e1;
t167 = qJ(1) + qJ(2);
t165 = sin(t167);
t166 = cos(t167);
t220 = sin(qJ(1)) * pkin(1);
t154 = t166 * qJ(3);
t168 = sin(qJ(4));
t170 = cos(qJ(4));
t183 = t168 * rSges(5,1) + t170 * rSges(5,2);
t172 = -t165 * rSges(5,3) + t166 * t183;
t244 = -pkin(2) - pkin(6);
t94 = t165 * t244 + t154 + t172;
t92 = t94 - t220;
t219 = cos(qJ(1)) * pkin(1);
t95 = (rSges(5,3) - t244) * t166 + (qJ(3) + t183) * t165;
t93 = t95 + t219;
t45 = t93 * t165 + t92 * t166;
t49 = t95 * t165 + t94 * t166;
t260 = pkin(2) - rSges(4,2);
t115 = t166 * rSges(4,3) - t260 * t165 + t154;
t111 = t115 - t220;
t116 = (rSges(4,3) + qJ(3)) * t165 + t260 * t166;
t112 = t116 + t219;
t60 = t111 * t166 + t112 * t165;
t68 = t115 * t166 + t116 * t165;
t216 = (t49 + t45) * t262 + (t68 + t60) * t263;
t217 = (t45 - t49) * t262 + (t60 - t68) * t263;
t4 = t217 - t216;
t265 = t4 * qJD(1);
t210 = Icges(5,4) * t170;
t146 = -Icges(5,2) * t168 + t210;
t182 = Icges(5,1) * t168 + t210;
t264 = t146 + t182;
t243 = m(3) * (t219 * (-t165 * rSges(3,1) - t166 * rSges(3,2)) + (t166 * rSges(3,1) - t165 * rSges(3,2)) * t220);
t239 = m(4) * (-t116 * t111 + t112 * t115);
t231 = m(5) * (-t95 * t92 + t93 * t94);
t179 = Icges(5,5) * t168 + Icges(5,6) * t170;
t259 = t165 * t179;
t258 = t166 * t179;
t211 = Icges(5,4) * t168;
t181 = Icges(5,2) * t170 + t211;
t120 = -Icges(5,6) * t165 + t166 * t181;
t122 = -Icges(5,5) * t165 + t166 * t182;
t257 = t170 * t120 + t168 * t122;
t119 = Icges(5,6) * t166 + t165 * t181;
t153 = t165 * t210;
t205 = t165 * t168;
t121 = Icges(5,1) * t205 + Icges(5,5) * t166 + t153;
t256 = t170 * t119 + t168 * t121;
t151 = t170 * rSges(5,1) - t168 * rSges(5,2);
t130 = t151 * t165;
t131 = t151 * t166;
t91 = -t166 * t130 + t165 * t131;
t213 = t95 * t130;
t39 = t94 * t131 + t213;
t255 = m(5) * t39;
t254 = t257 * t166;
t148 = Icges(5,1) * t170 - t211;
t253 = t168 * t264 - t170 * (-t181 + t148);
t192 = t146 * t166 + t122;
t194 = -t148 * t166 + t120;
t252 = t168 * t194 - t170 * t192;
t193 = -Icges(5,2) * t205 + t121 + t153;
t195 = -t148 * t165 + t119;
t251 = t168 * t195 - t170 * t193;
t186 = -t264 * t170 / 0.2e1 + (t181 / 0.2e1 - t148 / 0.2e1) * t168;
t163 = t165 ^ 2;
t164 = t166 ^ 2;
t250 = 0.4e1 * qJD(1);
t237 = m(4) * t60;
t236 = m(4) * t68;
t214 = t93 * t130;
t36 = t92 * t131 + t214;
t235 = m(5) * (t39 + t36);
t234 = m(5) * (t214 - t213 + (t92 - t94) * t131);
t229 = m(5) * t36;
t227 = m(5) * t45;
t226 = m(5) * t49;
t225 = m(5) * t91;
t224 = t165 / 0.2e1;
t223 = -t166 / 0.2e1;
t222 = t166 / 0.2e1;
t117 = Icges(5,3) * t166 + t259;
t108 = t165 * t117;
t118 = -Icges(5,3) * t165 + t258;
t50 = t166 * t117 + t256 * t165;
t51 = -t166 * t118 - t257 * t165;
t52 = -t166 * t256 + t108;
t53 = -t165 * t118 + t254;
t11 = (-t52 + t108 + t51) * t165 + (t53 - t254 + (t118 - t256) * t165 + t50) * t166;
t12 = t163 * t118 + (-t108 + t51 + (t118 + t256) * t166) * t166;
t24 = t51 * t165 + t50 * t166;
t25 = t53 * t165 + t52 * t166;
t2 = (t12 / 0.2e1 + t25 / 0.2e1) * t166 + (-t24 / 0.2e1 + t11 / 0.2e1) * t165;
t221 = -qJD(3) * t225 / 0.2e1 + t2 * qJD(4);
t189 = t225 / 0.2e1;
t188 = qJD(1) / 0.4e1 + qJD(2) / 0.4e1;
t187 = (t163 + t164) * t183;
t185 = t91 * t151;
t184 = t235 / 0.2e1 + t186;
t180 = Icges(5,5) * t170 - Icges(5,6) * t168;
t124 = t165 * t180;
t174 = -t165 * t11 / 0.2e1 + (t12 + t25) * t223 + (-t253 * t165 - t168 * t193 - t170 * t195 - t258) * t222 + (t253 * t166 + t168 * t192 + t170 * t194 + t24 - t259) * t224;
t173 = -t186 + (t222 + t223) * (t168 * t120 - t170 * t122);
t125 = t180 * t166;
t81 = qJD(3) * t189;
t80 = qJD(4) * t189;
t33 = t226 + t236;
t30 = t186 + t255;
t29 = t227 + t237;
t28 = t186 + t229;
t18 = t234 / 0.2e1;
t13 = t231 + t239 + t243;
t8 = -t234 / 0.2e1 + t184;
t7 = t18 + t184;
t5 = t216 + t217;
t3 = t18 - t235 / 0.2e1 + t173;
t1 = [t13 * qJD(2) + t29 * qJD(3) + t28 * qJD(4), t13 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + 0.2e1 * (t231 / 0.2e1 + t239 / 0.2e1 + t243 / 0.2e1) * qJD(2), t29 * qJD(1) + t5 * qJD(2) + t80, t28 * qJD(1) + t7 * qJD(2) + t81 + ((-(t165 * t92 - t166 * t93) * t183 + t185) * m(5) + t174) * qJD(4); -t4 * qJD(3) + t8 * qJD(4) + (-t243 / 0.4e1 - t239 / 0.4e1 - t231 / 0.4e1) * t250, t33 * qJD(3) + t30 * qJD(4), t33 * qJD(2) - t265 + t80, t8 * qJD(1) + t30 * qJD(2) + t81 + ((-(t165 * t94 - t166 * t95) * t183 + t185) * m(5) + t174) * qJD(4); t4 * qJD(2) + t80 + (-t237 / 0.4e1 - t227 / 0.4e1) * t250, t265 + t80 + 0.4e1 * (-t236 / 0.4e1 - t226 / 0.4e1) * qJD(2), 0, 0.2e1 * (t188 * t91 - qJD(4) * t187 / 0.2e1) * m(5); (t173 - t229) * qJD(1) + t3 * qJD(2) + t221, t3 * qJD(1) + (t173 - t255) * qJD(2) + t221, -0.2e1 * t188 * t225, (m(5) * ((-t166 * t172 + (-t166 * rSges(5,3) - t165 * t183) * t165) * (-t165 * t130 - t166 * t131) - t151 * t187) + (t164 * t124 + (t252 * t165 + (-t125 - t251) * t166) * t165) * t222 + (-t163 * t125 + (t251 * t166 + (t124 - t252) * t165) * t166) * t224) * qJD(4) + (qJD(1) + qJD(2)) * t2;];
Cq = t1;
