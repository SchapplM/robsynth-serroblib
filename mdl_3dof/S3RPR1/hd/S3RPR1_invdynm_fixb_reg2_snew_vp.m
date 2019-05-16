% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3RPR1
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
%   pkin=[a2,a3,d1,d3]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3RPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:56
% EndTime: 2019-05-04 18:28:57
% DurationCPUTime: 0.56s
% Computational Cost: add. (1118->102), mult. (1679->95), div. (0->0), fcn. (730->4), ass. (0->57)
t153 = sin(qJ(1));
t155 = cos(qJ(1));
t152 = sin(qJ(3));
t149 = -qJD(1) + qJD(3);
t147 = t149 ^ 2;
t148 = qJDD(1) - qJDD(3);
t154 = cos(qJ(3));
t160 = t152 * t147 + t154 * t148;
t175 = pkin(4) * t160 + t152 * g(3);
t134 = t154 * t147 - t152 * t148;
t180 = t155 * t134 + t153 * t160;
t183 = pkin(4) * t134 + t154 * g(3);
t189 = -pkin(3) * t180 + t153 * t175 + t155 * t183;
t162 = t153 * t134 - t155 * t160;
t187 = -pkin(3) * t162 + t153 * t183 - t155 * t175;
t169 = (qJD(2) * qJD(1));
t145 = 2 * t169;
t157 = qJD(1) ^ 2;
t142 = t155 * g(1) + t153 * g(2);
t150 = qJDD(1) * qJ(2);
t159 = -t142 + t150;
t174 = pkin(2) + pkin(1);
t121 = -t174 * t157 + t145 + t159;
t141 = t153 * g(1) - t155 * g(2);
t165 = -qJDD(2) + t141;
t158 = -t157 * qJ(2) - t165;
t126 = -t174 * qJDD(1) + t158;
t112 = t152 * t121 - t154 * t126;
t113 = t154 * t121 + t152 * t126;
t108 = t154 * t112 - t152 * t113;
t164 = t152 * t112 + t154 * t113;
t182 = t153 * t108 - t155 * t164;
t181 = t155 * t108 + t153 * t164;
t173 = pkin(4) * t108;
t170 = qJDD(1) * pkin(1);
t127 = t157 * pkin(1) - t159 - (2 * t169);
t129 = -t158 + t170;
t168 = -t155 * t127 - t153 * t129;
t167 = -t153 * t141 - t155 * t142;
t166 = -pkin(2) * g(3) + pkin(4) * t164;
t139 = t153 * qJDD(1) + t155 * t157;
t131 = -pkin(3) * t139 + t155 * g(3);
t140 = t155 * qJDD(1) - t153 * t157;
t130 = pkin(3) * t140 + t153 * g(3);
t163 = t153 * t127 - t155 * t129;
t161 = t155 * t141 - t153 * t142;
t156 = pkin(1) * g(3);
t151 = qJ(2) * g(3);
t137 = t165 + 0.2e1 * t170;
t132 = -t142 + t145 + 0.2e1 * t150;
t114 = pkin(1) * t129 - qJ(2) * t127;
t111 = qJ(2) * t160 + t134 * t174 + t113;
t110 = -qJ(2) * t134 + t160 * t174 + t112;
t105 = t151 + t173;
t104 = t156 - t166;
t103 = qJ(2) * t164 + t108 * t174;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t140, 0, -t139, 0, -t130, -t131, -t161, -pkin(3) * t161, 0, t140, 0, 0, t139, 0, -t130, t163, t131, pkin(3) * t163 + (-t153 * pkin(1) + t155 * qJ(2)) * g(3), 0, 0, t162, 0, -t180, 0, -t187, t189, t181, -pkin(3) * t181 - t153 * t104 + t155 * t105; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t139, 0, t140, 0, t131, -t130, t167, pkin(3) * t167, 0, t139, 0, 0, -t140, 0, t131, t168, t130, pkin(3) * t168 + (t155 * pkin(1) + t153 * qJ(2)) * g(3), 0, 0, -t180, 0, -t162, 0, t189, t187, t182, -pkin(3) * t182 + t155 * t104 + t153 * t105; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t141, t142, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t137, 0, t132, t114, 0, 0, 0, 0, 0, t148, t110, t111, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t157, 0, 0, -g(3), -t141, 0, 0, qJDD(1), 0, 0, t157, 0, 0, -t129, g(3), t151, 0, 0, -t160, 0, -t134, 0, t175, t183, t108, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, qJDD(1), 0, g(3), 0, -t142, 0, 0, t157, 0, 0, -qJDD(1), 0, g(3), -t127, 0, t156, 0, 0, -t134, 0, t160, 0, t183, -t175, -t164, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t141, t142, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t137, 0, t132, t114, 0, 0, 0, 0, 0, t148, t110, t111, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t157, 0, 0, -t129, g(3), 0, 0, 0, -t160, 0, -t134, 0, t175, t183, t108, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t129, 0, -t127, 0, 0, 0, 0, 0, 0, t148, pkin(2) * t160 + t112, pkin(2) * t134 + t113, 0, pkin(2) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, 0, 0, qJDD(1), 0, -g(3), t127, 0, 0, 0, 0, t134, 0, -t160, 0, -t183, t175, t164, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, 0, -t147, 0, 0, g(3), t112, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, -t148, 0, -g(3), 0, t113, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t112, -t113, 0, 0;];
m_new_reg  = t1;
