% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3PRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3PRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:25:02
% EndTime: 2019-05-04 18:25:02
% DurationCPUTime: 0.55s
% Computational Cost: add. (1085->92), mult. (1416->79), div. (0->0), fcn. (970->4), ass. (0->52)
t141 = qJD(2) + qJD(3);
t139 = t141 ^ 2;
t140 = qJDD(2) + qJDD(3);
t143 = sin(qJ(3));
t145 = cos(qJ(3));
t126 = t139 * t145 + t140 * t143;
t129 = t139 * t143 - t140 * t145;
t144 = sin(qJ(2));
t146 = cos(qJ(2));
t102 = t126 * t146 - t129 * t144;
t118 = pkin(4) * t126 - g(3) * t145;
t170 = pkin(4) * t129 - g(3) * t143;
t94 = pkin(3) * t102 + t118 * t146 - t144 * t170;
t106 = t126 * t144 + t129 * t146;
t174 = pkin(3) * t106 + t118 * t144 + t146 * t170;
t142 = g(2) - qJDD(1);
t133 = g(1) * t144 - t142 * t146;
t131 = qJDD(2) * pkin(2) + t133;
t134 = g(1) * t146 + t144 * t142;
t148 = qJD(2) ^ 2;
t132 = -pkin(2) * t148 - t134;
t108 = -t145 * t131 + t132 * t143;
t109 = t131 * t143 + t132 * t145;
t159 = t108 * t143 + t145 * t109;
t99 = t108 * t145 - t109 * t143;
t161 = t146 * t99;
t89 = -t144 * t159 + t161;
t162 = t144 * t99;
t171 = t146 * t159 + t162;
t112 = t133 * t146 - t134 * t144;
t167 = pkin(1) * t112;
t166 = pkin(3) * t112;
t163 = qJ(1) * g(3);
t96 = pkin(2) * t99;
t160 = pkin(1) * t89 + t96;
t157 = -t133 * t144 - t146 * t134;
t136 = qJDD(2) * t146 - t144 * t148;
t156 = -pkin(3) * t136 - g(3) * t144;
t155 = -pkin(2) * t129 - t108;
t135 = qJDD(2) * t144 + t146 * t148;
t154 = pkin(1) * t135 - t134;
t153 = -pkin(2) * t126 - t109;
t152 = pkin(1) * t106 - t155;
t151 = -pkin(1) * t136 - t133;
t95 = pkin(2) * g(3) + pkin(4) * t159;
t150 = pkin(3) * t89 + pkin(4) * t161 - t144 * t95;
t149 = pkin(1) * t102 - t153;
t147 = pkin(1) * g(3);
t123 = -pkin(3) * t135 + g(3) * t146;
t110 = pkin(3) * t157 + t147;
t86 = pkin(3) * t171 + pkin(4) * t162 + t146 * t95 + t147;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t142, -g(3), -t163, 0, 0, t136, 0, -t135, 0, t156, -t123, -t112, -t163 - t166, 0, 0, -t106, 0, -t102, 0, t174, t94, t89, t150 - t163; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t135, 0, t136, 0, t123, t156, t157, t110, 0, 0, t102, 0, -t106, 0, -t94, t174, t171, t86; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t142, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), qJ(1) * t135 - t151, qJ(1) * t136 - t154, 0, -qJ(1) * t157 + t167, 0, 0, 0, 0, 0, t140, qJ(1) * t102 - t152, -qJ(1) * t106 - t149, 0, -qJ(1) * t171 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -g(3), 0, 0, 0, t136, 0, -t135, 0, t156, -t123, -t112, -t166, 0, 0, -t106, 0, -t102, 0, t174, t94, t89, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), t151, t154, 0, -t167, 0, 0, 0, 0, 0, -t140, t152, t149, 0, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t135, 0, t136, 0, t123, t156, t157, t110, 0, 0, t102, 0, -t106, 0, -t94, t174, t171, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t148, 0, 0, -g(3), -t133, 0, 0, 0, -t129, 0, -t126, 0, t170, t118, t99, pkin(4) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, qJDD(2), 0, g(3), 0, -t134, 0, 0, 0, t126, 0, -t129, 0, -t118, t170, t159, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t133, t134, 0, 0, 0, 0, 0, 0, 0, t140, t155, t153, 0, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, -t139, 0, 0, -g(3), t108, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t140, 0, g(3), 0, t109, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t108, -t109, 0, 0;];
m_new_reg  = t1;
