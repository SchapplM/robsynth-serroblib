% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:39
% EndTime: 2019-12-31 17:35:42
% DurationCPUTime: 2.09s
% Computational Cost: add. (11436->188), mult. (11995->270), div. (0->0), fcn. (14086->8), ass. (0->138)
t251 = m(6) / 0.2e1;
t164 = sin(qJ(5));
t166 = cos(qJ(5));
t224 = Icges(6,4) * t166;
t183 = -Icges(6,2) * t164 + t224;
t184 = Icges(6,1) * t164 + t224;
t272 = t183 + t184;
t225 = Icges(6,4) * t164;
t182 = Icges(6,2) * t166 + t225;
t185 = Icges(6,1) * t166 - t225;
t190 = t272 * t166 / 0.2e1 + (-t182 / 0.2e1 + t185 / 0.2e1) * t164;
t207 = qJ(3) + qJ(4);
t196 = sin(t207);
t197 = cos(t207);
t226 = sin(pkin(8));
t227 = cos(pkin(8));
t134 = -t196 * t226 - t197 * t227;
t186 = rSges(6,1) * t164 + rSges(6,2) * t166;
t263 = t186 * t134;
t135 = t196 * t227 - t197 * t226;
t118 = t186 * t135;
t167 = cos(qJ(3));
t270 = pkin(3) * t167;
t258 = t270 * t226;
t165 = sin(qJ(3));
t195 = t227 * t165;
t160 = pkin(3) * t195;
t230 = rSges(6,2) * t164;
t192 = t134 * rSges(6,3) - t135 * t230;
t231 = rSges(6,1) * t166;
t198 = pkin(4) + t231;
t83 = -t134 * pkin(7) - t135 * t198 - t192;
t266 = t83 - t160;
t73 = t258 + t266;
t67 = t73 * t118;
t194 = t226 * t165;
t170 = -pkin(3) * t194 - t270 * t227;
t133 = t135 * pkin(7);
t191 = -t135 * rSges(6,3) - t134 * t230;
t84 = t134 * t198 - t133 + t191;
t74 = t84 + t170;
t29 = -t263 * t74 + t67;
t268 = t263 * t84;
t32 = t83 * t118 - t268;
t187 = (t32 + t29) * t251 + t190;
t181 = Icges(6,5) * t166 - Icges(6,6) * t164;
t264 = t135 * t181;
t99 = -Icges(6,3) * t134 - t264;
t269 = t135 * t99;
t102 = -Icges(6,6) * t134 - t135 * t183;
t216 = t164 * t102;
t265 = t134 * t181;
t98 = -Icges(6,3) * t135 + t265;
t267 = t98 + t216;
t101 = -Icges(6,6) * t135 + t134 * t183;
t104 = -Icges(6,5) * t135 + t134 * t185;
t260 = -t101 * t164 + t104 * t166;
t178 = t134 * t231 + t191;
t257 = qJD(3) + qJD(4);
t259 = t134 * rSges(5,1) + t135 * rSges(5,2);
t82 = t134 * pkin(4) - t133 + t178;
t26 = (t82 - t84) * t83;
t14 = m(6) * t26;
t255 = t134 ^ 2;
t254 = t135 ^ 2;
t253 = 2 * qJD(3);
t252 = m(5) / 0.2e1;
t121 = -t135 * rSges(5,1) + t134 * rSges(5,2);
t177 = -t121 + t160;
t108 = -t177 + t258;
t109 = t259 + t170;
t43 = t108 * t259 - t109 * t121;
t250 = m(5) * t43;
t71 = t82 * t263;
t247 = m(6) * (t67 + t71 - (t134 * t74 + t135 * t83) * t186);
t22 = t71 - t268;
t246 = m(6) * t22;
t23 = t73 * t82 - t74 * t83;
t245 = m(6) * t23;
t241 = -t134 / 0.2e1;
t238 = t135 / 0.2e1;
t31 = m(5) * (t121 * t227 + t226 * t259) + m(6) * (t226 * t82 + t227 * t83);
t168 = (-t134 * t227 - t135 * t226) * t186 * t251;
t172 = m(6) * (t118 * t226 + t227 * t263);
t57 = t172 / 0.2e1 - t168;
t237 = t31 * qJD(4) + t57 * qJD(5);
t105 = -Icges(6,5) * t134 - t135 * t185;
t212 = t166 * t105;
t234 = t134 * t99 + t135 * t212;
t233 = -t134 * t212 + t269;
t232 = m(6) * qJD(5);
t220 = t263 * t186;
t155 = t230 - t231;
t218 = t134 * t155;
t56 = -t172 / 0.2e1 - t168;
t208 = t56 * qJD(2);
t206 = -t184 * t135 + t102;
t205 = -t184 * t134 - t101;
t204 = -t182 * t135 - t105;
t203 = t182 * t134 - t104;
t33 = t135 * t216 - t234;
t34 = t134 * t98 + t260 * t135;
t19 = -t134 * t33 + t135 * t34;
t35 = t134 * t216 + t233;
t20 = -t134 * t35 + t135 * (t260 * t134 - t135 * t98);
t8 = (t267 * t135 - t234 - t33) * t135 + (-t34 + (-t212 + t267) * t134 + t269) * t134;
t9 = (t34 - t35 + t233) * t135 + t234 * t134;
t2 = (t9 / 0.2e1 - t19 / 0.2e1) * t135 + (-t20 / 0.2e1 - t8 / 0.2e1) * t134;
t202 = t2 * qJD(5) + t208;
t193 = t226 * t167;
t189 = pkin(3) * t193;
t180 = Icges(6,5) * t164 + Icges(6,6) * t166;
t173 = (t182 - t185) * t166 + t272 * t164;
t176 = (t135 * t173 - t265) * t241 - (t164 * t205 - t166 * t203 + t9) * t135 / 0.2e1 + (t134 * t173 + t19 + t264) * t238 + (t164 * t206 + t166 * t204 + t20 + t8) * t134 / 0.2e1;
t145 = -t167 * t227 - t194;
t175 = t164 * t204 - t166 * t206;
t174 = t164 * t203 + t166 * t205;
t171 = t145 * pkin(3);
t146 = -t193 + t195;
t113 = t180 * t134;
t112 = t180 * t135;
t111 = -t189 + t177;
t110 = t171 + t259;
t96 = t118 * t263;
t76 = t171 + t82;
t75 = -t189 - t266;
t72 = t118 * t135 + t134 * t263;
t54 = t56 * qJD(5);
t28 = m(6) * t32 + t190;
t27 = m(6) * t29 + t190;
t21 = t246 / 0.2e1;
t15 = t247 / 0.2e1;
t13 = t14 * qJD(4);
t10 = t245 + t250;
t5 = t21 - t247 / 0.2e1 + t187;
t4 = t15 - t246 / 0.2e1 + t187;
t3 = t15 + t21 - t187;
t1 = [0, 0, 0, 0, t232 * t72; 0, 0, (m(4) * ((rSges(4,1) * t145 + rSges(4,2) * t146) * t226 + (-rSges(4,1) * t146 + rSges(4,2) * t145) * t227) / 0.2e1 + (t110 * t226 - t111 * t227) * t252 + (t226 * t76 - t227 * t75) * t251) * t253 + t237, qJD(3) * t31 + t237, (-t134 * t226 + t135 * t227) * t155 * t232 + t257 * t57; 0, -t54, (m(6) * (t73 * t76 + t74 * t75) + m(5) * (t108 * t110 + t109 * t111)) * qJD(3) + t10 * qJD(4) + t27 * qJD(5), t10 * qJD(3) + t4 * qJD(5) + (0.2e1 * (t26 + t23) * t251 + 0.2e1 * t43 * t252 - t14) * qJD(4), t27 * qJD(3) + t4 * qJD(4) + (m(6) * (-t73 * t218 + t96 + (-t155 * t74 - t220) * t135) + t176) * qJD(5) - t208; 0, -t54, t5 * qJD(5) + t13 + 0.4e1 * (-t245 / 0.4e1 - t250 / 0.4e1) * qJD(3) + ((t75 * t84 + t76 * t83 + t23) * t251 + (t110 * t121 + t111 * t259 + t43) * t252) * t253, qJD(3) * t14 + qJD(5) * t28 + t13, t5 * qJD(3) + t28 * qJD(4) + (m(6) * (-t83 * t218 + t96 + (-t155 * t84 - t220) * t135) + t176) * qJD(5) - t208; 0, t257 * t56, (-t190 + (t67 - (-t135 * t75 + (t74 - t76) * t134) * t186 - t29) * m(6)) * qJD(3) + t3 * qJD(4) + t202, t3 * qJD(3) + (-t190 + (t22 - t32) * m(6)) * qJD(4) + t202, (m(6) * ((t135 * (-t135 * t231 - t192) - t134 * t178) * t72 - (t254 + t255) * t155 * t186) + (t254 * t113 + (t175 * t134 + (-t112 + t174) * t135) * t134) * t238 + (t255 * t112 + (t174 * t135 + (-t113 + t175) * t134) * t135) * t241) * qJD(5) + t257 * t2;];
Cq = t1;
