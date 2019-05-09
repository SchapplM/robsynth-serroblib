% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:40:48
% EndTime: 2019-05-04 18:40:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (1663->114), mult. (1863->85), div. (0->0), fcn. (1584->4), ass. (0->65)
t146 = sin(pkin(5));
t144 = g(2) - qJDD(1);
t145 = g(1) - qJDD(2);
t147 = cos(pkin(5));
t128 = -t146 * t144 + t147 * t145;
t129 = -t147 * t144 - t146 * t145;
t148 = sin(qJ(4));
t149 = cos(qJ(4));
t105 = t149 * t128 + t148 * t129;
t106 = -t148 * t128 + t149 * t129;
t162 = t148 * t105 + t149 * t106;
t97 = t149 * t105 - t148 * t106;
t171 = t147 * t97;
t83 = -t146 * t162 + t171;
t172 = t146 * t97;
t186 = t147 * t162 + t172;
t150 = qJD(4) ^ 2;
t132 = t148 * qJDD(4) + t149 * t150;
t133 = t149 * qJDD(4) - t148 * t150;
t110 = t147 * t132 + t146 * t133;
t143 = g(3) + qJDD(3);
t123 = pkin(4) * t132 + t149 * t143;
t157 = -pkin(4) * t133 + t148 * t143;
t160 = t147 * t123 - t146 * t157;
t185 = qJ(3) * t110 + t160;
t113 = t146 * t132 - t147 * t133;
t161 = t146 * t123 + t147 * t157;
t165 = qJ(3) * t113 + t161;
t180 = pkin(2) * t83;
t179 = pkin(1) * t144;
t103 = t147 * t128 - t146 * t129;
t178 = pkin(2) * t103;
t177 = pkin(2) * t110;
t176 = pkin(2) * t113;
t175 = qJ(1) * g(3);
t174 = -pkin(1) - qJ(3);
t173 = -pkin(2) - qJ(1);
t156 = t146 * t128 + t147 * t129;
t170 = qJ(2) * t156;
t169 = qJ(2) * t144;
t135 = t146 * t143;
t168 = t147 * t143;
t142 = pkin(2) * t143;
t167 = qJ(3) * t156 - t142;
t164 = pkin(3) * t132 + t106;
t94 = pkin(3) * t97;
t163 = qJ(2) * t186 + t94;
t159 = -pkin(1) * t156 - t167;
t158 = -pkin(3) * t133 + t105;
t92 = -pkin(3) * t143 + pkin(4) * t162;
t155 = pkin(4) * t171 - t146 * t92;
t154 = qJ(2) * t113 + t164;
t153 = pkin(4) * t172 + qJ(3) * t186 + t147 * t92 - t142;
t152 = -qJ(2) * t110 + t158;
t151 = -pkin(1) * t186 - t153;
t140 = qJ(1) * t143;
t139 = qJ(2) * t143;
t134 = pkin(1) * t145 + qJ(2) * g(3);
t93 = -t103 * t174 + t139;
t90 = -t113 * t174 + t161;
t89 = -pkin(1) * t113 - t165;
t88 = pkin(1) * t110 + t185;
t87 = -t110 * t174 + t160;
t81 = -t174 * t83 + t139 + t155;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t144, 0, -g(3), -t175, 0, 0, 0, 0, 0, 0, -g(3), t144, 0, -t175 - t179, 0, 0, 0, 0, 0, 0, -t168, t135, t156, -t140 - t159, 0, 0, t110, 0, -t113, 0, -t88, -t89, t186, -t140 - t151; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, g(3), t134, 0, 0, 0, 0, 0, 0, t135, t168, t103, t93, 0, 0, -t113, 0, -t110, 0, t90, t87, t83, t81; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t145, 0, -t144, qJ(1) * t145 - t169, 0, 0, 0, 0, 0, 0, t128, t129, 0, -t103 * t173 + t170, 0, 0, 0, 0, 0, -qJDD(4), -t113 * t173 + t152, -t110 * t173 + t154, 0, -t173 * t83 + t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -g(1), 0, 0, 0, 0, 0, 0, 0, -t145, 0, t144, t169, 0, 0, 0, 0, 0, 0, -t128, -t129, 0, -t170 - t178, 0, 0, 0, 0, 0, qJDD(4), -t152 - t176, -t154 - t177, 0, -t163 - t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, g(3), 0, 0, 0, 0, 0, 0, 0, g(3), -t144, 0, t179, 0, 0, 0, 0, 0, 0, t168, -t135, -t156, t159, 0, 0, -t110, 0, t113, 0, t88, t89, -t186, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, g(3), t134, 0, 0, 0, 0, 0, 0, t135, t168, t103, t93, 0, 0, -t113, 0, -t110, 0, t90, t87, t83, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, g(3), 0, 0, 0, 0, 0, 0, 0, t135, t168, t103, qJ(3) * t103, 0, 0, -t113, 0, -t110, 0, t165, t185, t83, qJ(3) * t83 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, -t144, 0, 0, 0, 0, 0, 0, 0, t128, t129, 0, t178, 0, 0, 0, 0, 0, -qJDD(4), t158 + t176, t164 + t177, 0, t94 + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), t144, 0, 0, 0, 0, 0, 0, 0, 0, -t168, t135, t156, t167, 0, 0, t110, 0, -t113, 0, -t185, t165, t186, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t128, 0, 0, 0, t133, 0, -t132, 0, t157, t123, t97, pkin(4) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, 0, t129, 0, 0, 0, t132, 0, t133, 0, -t123, t157, t162, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t129, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t158, -t164, 0, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), 0, -t150, 0, 0, t143, t105, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, qJDD(4), 0, -t143, 0, t106, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -t105, -t106, 0, 0;];
m_new_reg  = t1;
