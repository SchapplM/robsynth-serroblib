% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:38
% DurationCPUTime: 0.94s
% Computational Cost: add. (3463->143), mult. (4236->158), div. (0->0), fcn. (2312->6), ass. (0->85)
t180 = sin(qJ(4));
t181 = sin(qJ(2));
t183 = cos(qJ(4));
t184 = cos(qJ(2));
t179 = (qJD(1) + qJD(2));
t173 = qJD(4) - t179;
t167 = t173 ^ 2;
t178 = qJDD(1) + qJDD(2);
t168 = -qJDD(4) + t178;
t190 = t167 * t180 + t168 * t183;
t194 = -t167 * t183 + t168 * t180;
t210 = t181 * t194 + t184 * t190;
t105 = pkin(5) * t210 - (t180 * t184 - t181 * t183) * g(3);
t182 = sin(qJ(1));
t185 = cos(qJ(1));
t117 = t181 * t190 - t184 * t194;
t175 = t184 * g(3);
t207 = t181 * g(3);
t224 = -pkin(5) * t117 + t175 * t183 + t180 * t207;
t93 = t117 * t185 + t182 * t210;
t232 = -pkin(4) * t93 - t105 * t182 + t185 * t224;
t177 = t179 ^ 2;
t165 = t178 * qJ(3);
t159 = g(1) * t185 + g(2) * t182;
t187 = qJD(1) ^ 2;
t155 = -pkin(1) * t187 - t159;
t158 = g(1) * t182 - g(2) * t185;
t188 = qJDD(1) * pkin(1) + t158;
t130 = t155 * t184 + t181 * t188;
t200 = (2 * qJD(3) * t179) + t130;
t193 = t165 + t200;
t208 = pkin(2) + pkin(3);
t109 = -t177 * t208 + t193;
t171 = t178 * pkin(2);
t129 = t155 * t181 - t184 * t188;
t191 = -qJDD(3) - t129;
t115 = -qJ(3) * t177 - t171 - t191;
t113 = -pkin(3) * t178 + t115;
t87 = t180 * t109 - t113 * t183;
t88 = t109 * t183 + t113 * t180;
t78 = t180 * t88 - t183 * t87;
t80 = t180 * t87 + t183 * t88;
t219 = t181 * t80 - t184 * t78;
t76 = t181 * t78 + t184 * t80;
t69 = t182 * t76 + t185 * t219;
t70 = -t182 * t219 + t185 * t76;
t223 = -t117 * t182 + t185 * t210;
t230 = pkin(4) * t223 + t105 * t185 + t182 * t224;
t150 = t177 * t184 + t178 * t181;
t153 = t177 * t181 - t178 * t184;
t122 = t150 * t185 - t153 * t182;
t134 = pkin(5) * t153 - t207;
t137 = pkin(5) * t150 - t175;
t229 = pkin(4) * t122 - t134 * t182 + t137 * t185;
t126 = t150 * t182 + t153 * t185;
t228 = pkin(4) * t126 + t134 * t185 + t137 * t182;
t197 = t129 * t181 + t130 * t184;
t102 = t129 * t184 - t130 * t181;
t204 = t185 * t102;
t225 = -t182 * t197 + t204;
t206 = t182 * t102;
t84 = t185 * t197 + t206;
t222 = pkin(1) * t150;
t114 = -pkin(2) * t177 + t193;
t198 = t114 * t184 + t115 * t181;
t90 = t114 * t181 - t115 * t184;
t82 = -t182 * t90 + t185 * t198;
t81 = t182 * t198 + t185 * t90;
t186 = pkin(1) * g(3);
t201 = qJ(3) * t207 + t186;
t132 = -t158 * t182 - t159 * t185;
t157 = qJDD(1) * t185 - t182 * t187;
t192 = -pkin(4) * t157 - g(3) * t182;
t131 = t158 * t185 - t159 * t182;
t170 = qJ(3) * t175;
t166 = t208 * g(3);
t156 = qJDD(1) * t182 + t185 * t187;
t146 = pkin(1) * t153;
t145 = -pkin(4) * t156 + g(3) * t185;
t99 = pkin(5) * t197 + t186;
t86 = -pkin(2) * t207 - pkin(5) * t90 + t170;
t85 = pkin(2) * t175 + pkin(5) * t198 + t201;
t72 = -pkin(5) * t219 - t166 * t181 + t170;
t71 = pkin(5) * t76 + t166 * t184 + t201;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t156, -t157, 0, t132, 0, 0, 0, 0, 0, 0, -t122, t126, 0, t84, 0, 0, 0, 0, 0, 0, -t122, 0, -t126, t82, 0, 0, 0, 0, 0, 0, -t93, t223, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t157, -t156, 0, t131, 0, 0, 0, 0, 0, 0, -t126, -t122, 0, -t225, 0, 0, 0, 0, 0, 0, -t126, 0, t122, t81, 0, 0, 0, 0, 0, 0, t223, t93, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t157, 0, -t156, 0, t192, -t145, -t131, -pkin(4) * t131, 0, 0, -t126, 0, -t122, 0, t228, t229, t225, pkin(4) * t225 + pkin(5) * t204 - t182 * t99, 0, -t126, 0, 0, t122, 0, t228, -t81, -t229, -pkin(4) * t81 - t182 * t85 + t185 * t86, 0, 0, -t223, 0, -t93, 0, -t230, t232, t69, -pkin(4) * t69 - t182 * t71 + t185 * t72; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t156, 0, t157, 0, t145, t192, t132, pkin(4) * t132, 0, 0, t122, 0, -t126, 0, -t229, t228, t84, pkin(4) * t84 + pkin(5) * t206 + t185 * t99, 0, t122, 0, 0, t126, 0, -t229, t82, -t228, pkin(4) * t82 + t182 * t86 + t185 * t85, 0, 0, -t93, 0, t223, 0, t232, t230, -t70, pkin(4) * t70 + t182 * t72 + t185 * t71; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t158, t159, 0, 0, 0, 0, 0, 0, 0, t178, -t129 - t146, -t130 - t222, 0, -pkin(1) * t102, 0, 0, 0, t178, 0, 0, -t146 + 0.2e1 * t171 + t191, 0, 0.2e1 * t165 + t200 + t222, pkin(1) * t90 - pkin(2) * t115 + qJ(3) * t114, 0, 0, 0, 0, 0, t168, pkin(1) * t210 + qJ(3) * t194 + t190 * t208 + t87, pkin(1) * t117 + qJ(3) * t190 - t194 * t208 + t88, 0, pkin(1) * t219 + qJ(3) * t80 - t208 * t78;];
tauB_reg  = t1;
