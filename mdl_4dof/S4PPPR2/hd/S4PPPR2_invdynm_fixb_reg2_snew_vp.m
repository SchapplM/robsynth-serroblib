% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPPR2
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
%   pkin=[a2,a3,a4,d4,theta2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:38:08
% EndTime: 2019-05-04 18:38:09
% DurationCPUTime: 0.80s
% Computational Cost: add. (1238->128), mult. (1500->89), div. (0->0), fcn. (1258->4), ass. (0->58)
t145 = sin(qJ(4));
t146 = cos(qJ(4));
t147 = qJD(4) ^ 2;
t131 = t145 * qJDD(4) + t146 * t147;
t132 = -t146 * qJDD(4) + t145 * t147;
t143 = sin(pkin(5));
t144 = cos(pkin(5));
t105 = t144 * t131 + t143 * t132;
t142 = g(3) - qJDD(2);
t113 = pkin(4) * t131 + t146 * t142;
t157 = pkin(4) * t132 + t145 * t142;
t183 = -qJ(2) * t105 + t144 * t113 + t143 * t157;
t166 = g(2) - qJDD(1);
t127 = t143 * g(1) - t144 * t166;
t124 = -qJDD(3) + t127;
t128 = t144 * g(1) + t143 * t166;
t100 = -t145 * t124 - t146 * t128;
t99 = t146 * t124 - t145 * t128;
t156 = t146 * t100 + t145 * t99;
t90 = t145 * t100 - t146 * t99;
t79 = t143 * t90 + t144 * t156;
t182 = t143 * t156 - t144 * t90;
t162 = -t143 * t131 + t144 * t132;
t84 = qJ(2) * t162 + t143 * t113 - t144 * t157;
t176 = pkin(3) * t90;
t175 = pkin(4) * t90;
t174 = pkin(4) * t156;
t117 = t143 * t128;
t102 = t144 * t127 - t117;
t173 = pkin(1) * t102;
t172 = qJ(2) * t102;
t171 = qJ(3) * t142;
t168 = t143 * t142;
t167 = t144 * t142;
t165 = pkin(2) * t124 - qJ(3) * t128;
t96 = t144 * t124 - t117;
t118 = t144 * t128;
t164 = -t143 * t124 - t118;
t163 = -t143 * t127 - t118;
t160 = -pkin(2) * t90 + qJ(3) * t156 - t176;
t159 = -pkin(1) * t96 - t165;
t158 = pkin(3) * t131 + t100;
t80 = -t174 + (pkin(2) + pkin(3)) * t142;
t82 = t171 - t175;
t155 = -qJ(2) * t182 - t143 * t80 + t144 * t82;
t154 = pkin(3) * t132 + t99;
t153 = -pkin(1) * t182 - t160;
t152 = pkin(2) * t131 + qJ(3) * t132 + t158;
t151 = -pkin(2) * t168 - qJ(2) * t96 + qJ(3) * t167;
t150 = pkin(2) * t132 - qJ(3) * t131 + t154;
t149 = -pkin(1) * t105 - t152;
t148 = -pkin(1) * t162 - t150;
t141 = pkin(1) * t142;
t139 = qJ(1) * t142;
t94 = qJ(2) * t163 + t141;
t93 = qJ(2) * t164 + t141 + (pkin(2) * t144 + qJ(3) * t143) * t142;
t75 = qJ(2) * t79 + t143 * t82 + t144 * t80 + t141;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t166, -g(3), -qJ(1) * g(3), 0, 0, 0, 0, 0, 0, -t168, -t167, -t102, -t139 - t172, 0, 0, 0, 0, 0, 0, -t168, -t96, t167, -t139 + t151, 0, 0, -t162, 0, -t105, 0, -t84, t183, t182, -t139 + t155; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t167, -t168, t163, t94, 0, 0, 0, 0, 0, 0, t167, t164, t168, t93, 0, 0, -t105, 0, t162, 0, t183, t84, -t79, t75; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t166, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t127, t128, 0, -qJ(1) * t163 + t173, 0, 0, 0, 0, 0, 0, t124, 0, -t128, -qJ(1) * t164 - t159, 0, 0, 0, 0, 0, -qJDD(4), qJ(1) * t105 - t148, -qJ(1) * t162 - t149, 0, -qJ(1) * t79 - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -g(3), 0, 0, 0, 0, 0, 0, 0, -t168, -t167, -t102, -t172, 0, 0, 0, 0, 0, 0, -t168, -t96, t167, t151, 0, 0, -t162, 0, -t105, 0, -t84, t183, t182, t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, -t127, -t128, 0, -t173, 0, 0, 0, 0, 0, 0, -t124, 0, t128, t159, 0, 0, 0, 0, 0, qJDD(4), t148, t149, 0, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t167, -t168, t163, t94, 0, 0, 0, 0, 0, 0, t167, t164, t168, t93, 0, 0, -t105, 0, t162, 0, t183, t84, -t79, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t127, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t142, t171, 0, 0, -t132, 0, -t131, 0, t157, t113, -t90, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, -t128, 0, 0, 0, 0, 0, 0, 0, t142, -t128, 0, pkin(2) * t142, 0, 0, -t131, 0, t132, 0, t113, -t157, -t156, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t128, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t128, t165, 0, 0, 0, 0, 0, -qJDD(4), t150, t152, 0, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t142, 0, 0, 0, -t132, 0, -t131, 0, t157, t113, -t90, -t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, 0, -t128, 0, 0, 0, 0, 0, 0, -qJDD(4), t154, t158, 0, -t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t128, 0, 0, 0, 0, t131, 0, -t132, 0, -t113, t157, t156, -pkin(3) * t142 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), 0, -t147, 0, 0, t142, t99, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, qJDD(4), 0, -t142, 0, t100, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t99, -t100, 0, 0;];
m_new_reg  = t1;
