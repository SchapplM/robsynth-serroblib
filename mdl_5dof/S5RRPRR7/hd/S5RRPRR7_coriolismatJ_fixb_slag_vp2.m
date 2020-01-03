% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:30
% EndTime: 2019-12-31 20:15:34
% DurationCPUTime: 1.91s
% Computational Cost: add. (4871->209), mult. (8262->261), div. (0->0), fcn. (7393->6), ass. (0->136)
t141 = cos(qJ(4));
t138 = t141 ^ 2;
t227 = pkin(4) * t141;
t263 = m(6) * t227;
t139 = sin(qJ(4));
t229 = sin(qJ(5));
t230 = cos(qJ(5));
t178 = -t230 * t139 - t229 * t141;
t179 = t229 * t139 - t230 * t141;
t166 = -t178 * mrSges(6,1) - t179 * mrSges(6,2);
t267 = t166 * qJD(5);
t134 = t139 * pkin(8);
t143 = -pkin(2) - pkin(7);
t207 = t139 * t143;
t118 = -t134 + t207;
t119 = (-pkin(8) + t143) * t141;
t84 = -t229 * t118 + t230 * t119;
t181 = t84 * mrSges(6,2);
t85 = -t230 * t118 - t229 * t119;
t261 = t85 * mrSges(6,1);
t266 = t261 - t181;
t142 = cos(qJ(2));
t228 = pkin(1) * t142;
t191 = -pkin(2) - t228;
t125 = -pkin(7) + t191;
t210 = t125 * t139;
t102 = -t134 + t210;
t103 = (-pkin(8) + t125) * t141;
t71 = -t229 * t102 + t230 * t103;
t182 = t71 * mrSges(6,2);
t72 = -t230 * t102 - t229 * t103;
t262 = t72 * mrSges(6,1);
t265 = t262 - t182;
t264 = t178 / 0.2e1;
t250 = mrSges(6,3) * t179;
t50 = t72 * t250;
t61 = t85 * t250;
t137 = t139 ^ 2;
t203 = t137 + t138;
t260 = t203 * mrSges(5,3) + mrSges(3,1);
t133 = t141 * mrSges(5,2);
t161 = -t139 * mrSges(5,1) - t133 - t166;
t160 = -mrSges(4,3) + t161;
t140 = sin(qJ(2));
t136 = t140 * pkin(1);
t127 = t136 + qJ(3);
t214 = t141 * mrSges(5,1);
t215 = t139 * mrSges(5,2);
t188 = t214 - t215;
t101 = t127 * t188;
t116 = qJ(3) * t188;
t135 = t139 * pkin(4);
t120 = t135 + t127;
t128 = qJ(3) + t135;
t169 = t178 * mrSges(6,3);
t173 = Ifges(6,2) * t178;
t83 = -t179 * mrSges(6,1) + t178 * mrSges(6,2);
t65 = t120 * t83;
t74 = t128 * t83;
t180 = -t50 / 0.2e1 - t61 / 0.2e1 + 0.2e1 * Ifges(6,4) * t178 * t264 + t65 / 0.2e1 + t74 / 0.2e1 + (t173 / 0.2e1 - Ifges(6,4) * t179 + Ifges(6,2) * t264 + (-t264 - t178 / 0.2e1) * Ifges(6,1)) * t179;
t257 = (Ifges(5,1) - Ifges(5,2)) * t141;
t51 = t71 * t169;
t62 = t84 * t169;
t77 = t166 * t227;
t258 = t180 + t101 / 0.2e1 + t116 / 0.2e1 + t51 / 0.2e1 + t62 / 0.2e1 + t77 + (t120 + t128) * t263 / 0.2e1 - t138 * Ifges(5,4) + (t72 + t85) * t250 / 0.2e1 - (t71 + t84) * t169 / 0.2e1 + (-t257 / 0.2e1 + Ifges(5,4) * t139 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t141) * t139;
t194 = t136 / 0.2e1;
t256 = -t181 / 0.2e1;
t255 = -t182 / 0.2e1;
t254 = t50 - t65;
t253 = t61 - t74;
t100 = t178 * pkin(4) * t230;
t246 = m(5) / 0.4e1 + m(4) / 0.4e1;
t157 = t179 * t178;
t236 = mrSges(6,3) / 0.2e1;
t147 = -t157 * t236 + t250 * t264;
t245 = qJD(3) * t147;
t244 = qJD(5) * t147;
t92 = t179 * t136;
t93 = t178 * t136;
t218 = -t92 * mrSges(6,1) / 0.2e1 + t93 * mrSges(6,2) / 0.2e1;
t243 = m(4) / 0.2e1;
t241 = m(5) / 0.2e1;
t238 = m(6) / 0.2e1;
t237 = m(6) * pkin(4);
t235 = t72 / 0.2e1;
t234 = t85 / 0.2e1;
t213 = t142 * mrSges(3,2);
t170 = t179 * t92;
t158 = mrSges(4,2) * t136 - mrSges(6,3) * t170 - t160 * t228 - t93 * t169;
t190 = t203 * t140;
t209 = t127 * t142;
t9 = -m(6) * (t120 * t228 - t71 * t92 + t72 * t93) - t158 + (t213 - m(5) * (t125 * t190 + t209) - m(4) * (t191 * t140 + t209)) * pkin(1) + t260 * t136;
t212 = t9 * qJD(1);
t36 = m(6) * t120 + 0.4e1 * t246 * t127 - t160;
t211 = qJD(1) * t36;
t206 = Ifges(6,5) * t178 + Ifges(6,6) * t179;
t201 = t194 + qJ(3);
t197 = pkin(4) * t229;
t189 = pkin(1) * t190;
t146 = (-t178 ^ 2 + t179 ^ 2) * Ifges(6,4) + Ifges(6,1) * t157 - t179 * t173;
t7 = t146 + t254 - t50;
t187 = -t7 * qJD(1) - t245;
t145 = (t138 - t137) * Ifges(5,4) - t77 + t146 + t139 * t257;
t5 = -t51 - t120 * t263 - t101 + (t71 * t178 - t179 * t72) * mrSges(6,3) + t145 + t254;
t186 = t5 * qJD(1);
t185 = (qJD(1) + qJD(2)) * t147;
t172 = t178 * t93;
t29 = (t189 / 0.2e1 - t201) * m(5) + (t194 - t201) * m(4) + (t172 / 0.2e1 + t170 / 0.2e1 - t135 - t201) * m(6) + t160;
t42 = m(6) * t128 + (m(5) + m(4)) * qJ(3) - t160;
t183 = -qJD(1) * t29 + qJD(2) * t42;
t153 = (-t229 * t93 - t230 * t92) * t237 / 0.2e1 - t215 * t136 / 0.2e1 + t194 * t214 + t218;
t2 = t153 - t258;
t6 = -t62 - t116 - t128 * t263 + (t84 * t178 - t179 * t85) * mrSges(6,3) + t145 + t253;
t177 = -t2 * qJD(1) - t6 * qJD(2);
t148 = t180 + (t234 + t235) * t250;
t4 = t148 - t218;
t8 = t146 + t253 - t61;
t176 = t4 * qJD(1) - t8 * qJD(2) - t245;
t165 = t179 * t197;
t117 = (t229 * mrSges(6,1) + t230 * mrSges(6,2)) * pkin(4);
t98 = t250 * t197;
t159 = t98 / 0.2e1 - t165 * t236;
t156 = t255 + t159;
t12 = t182 / 0.2e1 + (t235 - t72 / 0.2e1) * mrSges(6,1) + t156;
t155 = t256 + t159;
t19 = t181 / 0.2e1 + (t234 - t85 / 0.2e1) * mrSges(6,1) + t155;
t162 = -t12 * qJD(1) - t19 * qJD(2) + t117 * qJD(4);
t154 = -mrSges(6,3) * t100 - Ifges(5,5) * t139 - Ifges(5,6) * t141 + t206 + t98;
t129 = qJ(3) * t228;
t114 = t117 * qJD(5);
t30 = m(6) * t135 + (t172 + t170) * t238 + (t203 * t241 + t243) * t136 - t160 + (m(6) / 0.4e1 + t246) * (0.2e1 * t136 + 0.4e1 * qJ(3));
t18 = t261 + t256 + t155 + t206;
t10 = t262 + t255 + t156 + t206;
t3 = t148 + t218;
t1 = t153 + t258;
t11 = [-qJD(2) * t9 + qJD(3) * t36 - qJD(4) * t5 - qJD(5) * t7, -t212 + t158 * qJD(2) + t30 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t260 * t140 - t213) * qJD(2) * pkin(1) + 0.2e1 * ((t128 * t228 - t84 * t92 + t85 * t93) * t238 + (t143 * t189 + t129) * t241 + (-pkin(2) * t136 + t129) * t243) * qJD(2), qJD(2) * t30 + t211 - t244, t1 * qJD(2) + (-mrSges(5,1) * t210 - t125 * t133 + (t229 * t71 + t230 * t72) * t237 + t154 + t265) * qJD(4) + t10 * qJD(5) - t186, t3 * qJD(2) + t10 * qJD(4) + (t206 + t265) * qJD(5) + t187; -qJD(3) * t29 - qJD(4) * t2 + qJD(5) * t4 + t212, qJD(3) * t42 - qJD(4) * t6 - qJD(5) * t8, t183 - t244, (-mrSges(5,1) * t207 - t143 * t133 + (t229 * t84 + t230 * t85) * t237 + t154 + t266) * qJD(4) + t18 * qJD(5) + t177, t18 * qJD(4) + (t206 + t266) * qJD(5) + t176; qJD(2) * t29 - t211 - t244, -t183 - t244, 0, (m(6) * (-t165 + t100) + t161) * qJD(4) - t267, -qJD(4) * t166 - t185 - t267; qJD(2) * t2 + qJD(5) * t12 + t186, qJD(5) * t19 - t177, 0, -t114, -t114 - t162; -qJD(2) * t4 - qJD(4) * t12 - t187, -qJD(4) * t19 - t176, t185, t162, 0;];
Cq = t11;
