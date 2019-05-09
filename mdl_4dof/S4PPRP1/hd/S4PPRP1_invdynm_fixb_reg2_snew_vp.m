% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PPRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:43:16
% EndTime: 2019-05-04 18:43:17
% DurationCPUTime: 0.94s
% Computational Cost: add. (1279->150), mult. (1956->117), div. (0->0), fcn. (1384->4), ass. (0->71)
t175 = sin(qJ(3));
t176 = cos(qJ(3));
t177 = qJD(3) ^ 2;
t158 = t175 * qJDD(3) + t176 * t177;
t159 = -t176 * qJDD(3) + t175 * t177;
t173 = sin(pkin(5));
t174 = cos(pkin(5));
t133 = -t173 * t158 + t174 * t159;
t171 = g(3) - qJDD(1);
t202 = t175 * t171;
t141 = pkin(4) * t159 + t202;
t201 = t176 * t171;
t143 = pkin(4) * t158 + t201;
t229 = qJ(1) * t133 - t174 * t141 + t173 * t143;
t193 = t174 * t158 + t173 * t159;
t119 = -qJ(1) * t193 + t173 * t141 + t174 * t143;
t199 = (qJD(4) * qJD(3));
t160 = t173 * g(1) - t174 * g(2);
t157 = -qJDD(2) + t160;
t161 = t174 * g(1) + t173 * g(2);
t200 = t175 * t157 + t176 * t161;
t183 = (2 * t199) - t200;
t198 = qJDD(3) * qJ(4);
t131 = -t177 * pkin(3) + t183 + t198;
t172 = qJDD(3) * pkin(3);
t194 = -t176 * t157 + t175 * t161;
t182 = -qJDD(4) + t194;
t132 = -t177 * qJ(4) - t172 - t182;
t116 = t175 * t131 - t176 * t132;
t181 = t176 * t131 + t175 * t132;
t225 = t173 * t116 + t174 * t181;
t224 = -t174 * t116 + t173 * t181;
t129 = t175 * t200 - t176 * t194;
t180 = -t175 * t194 - t176 * t200;
t223 = t173 * t129 - t174 * t180;
t222 = t174 * t129 + t173 * t180;
t216 = pkin(1) * t159 - qJ(2) * t158;
t215 = -0.2e1 * t198 - (2 * t199);
t213 = pkin(1) + pkin(2);
t125 = qJ(2) * t159 + t213 * t158 - t200;
t212 = pkin(1) * t171;
t211 = pkin(4) * t129;
t209 = qJ(2) * t171;
t205 = t173 * t171;
t163 = t174 * t171;
t197 = pkin(3) * t202 - pkin(4) * t116;
t165 = pkin(2) * t171;
t196 = -pkin(4) * t181 + t165;
t195 = pkin(4) * t180 - t165;
t152 = t173 * t161;
t191 = t174 * t157 - t152;
t190 = t174 * t160 - t152;
t153 = t174 * t161;
t189 = -t173 * t157 - t153;
t188 = -t173 * t160 - t153;
t187 = pkin(2) * t158 - t200;
t185 = pkin(3) * t132 - qJ(4) * t131;
t184 = -pkin(3) * t176 - qJ(4) * t175;
t179 = pkin(2) * t159 - t194;
t178 = qJDD(4) - 0.2e1 * t172 + t179;
t139 = pkin(1) * t157 - qJ(2) * t161;
t126 = t179 + t216;
t124 = t178 + t216;
t123 = -t125 + t215;
t122 = t209 + t211;
t121 = -t195 + t212;
t114 = (-qJ(4) * t176 + qJ(2)) * t171 + t197;
t113 = (pkin(1) - t184) * t171 + t196;
t112 = qJ(2) * t180 + t129 * t213;
t111 = qJ(2) * t181 - t213 * t116 + t185;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t205, -t163, -t190, -qJ(1) * t190, 0, 0, 0, 0, 0, 0, -t205, -t191, t163, -qJ(1) * t191 + (-t173 * pkin(1) + t174 * qJ(2)) * t171, 0, 0, -t133, 0, -t193, 0, -t229, t119, t222, -qJ(1) * t222 - t173 * t121 + t174 * t122, 0, -t133, 0, 0, t193, 0, -t229, t224, -t119, -qJ(1) * t224 - t173 * t113 + t174 * t114; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t163, -t205, t188, qJ(1) * t188, 0, 0, 0, 0, 0, 0, t163, t189, t205, qJ(1) * t189 + (t174 * pkin(1) + t173 * qJ(2)) * t171, 0, 0, -t193, 0, t133, 0, t119, t229, t223, -qJ(1) * t223 + t174 * t121 + t173 * t122, 0, -t193, 0, 0, -t133, 0, t119, -t225, -t229, qJ(1) * t225 + t174 * t113 + t173 * t114; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t160, t161, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, -t161, t139, 0, 0, 0, 0, 0, -qJDD(3), t126, t125, 0, t112, 0, 0, 0, -qJDD(3), 0, 0, t124, 0, t123, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t160, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t171, t209, 0, 0, -t159, 0, -t158, 0, t141, t143, t129, t122, 0, -t159, 0, 0, t158, 0, t141, -t116, -t143, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, -t161, 0, 0, 0, 0, 0, 0, 0, t171, -t161, 0, t212, 0, 0, -t158, 0, t159, 0, t143, -t141, -t180, t121, 0, -t158, 0, 0, -t159, 0, t143, -t181, t141, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t161, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, -t161, t139, 0, 0, 0, 0, 0, -qJDD(3), t126, t125, 0, t112, 0, 0, 0, -qJDD(3), 0, 0, t124, 0, t123, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t171, 0, 0, 0, -t159, 0, -t158, 0, t141, t143, t129, t211, 0, -t159, 0, 0, t158, 0, t141, -t116, -t143, -qJ(4) * t201 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, -t161, 0, 0, 0, 0, 0, 0, -qJDD(3), t179, t187, 0, pkin(2) * t129, 0, 0, 0, -qJDD(3), 0, 0, t178, 0, -t187 + t215, -pkin(2) * t116 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t161, 0, 0, 0, 0, t158, 0, -t159, 0, -t143, t141, t180, t195, 0, t158, 0, 0, t159, 0, -t143, t181, -t141, t184 * t171 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t177, 0, 0, t171, -t194, 0, 0, qJDD(3), 0, 0, t177, 0, 0, t132, -t171, -qJ(4) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, 0, qJDD(3), 0, -t171, 0, -t200, 0, 0, t177, 0, 0, -qJDD(3), 0, -t171, t131, 0, -pkin(3) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t194, t200, 0, 0, 0, 0, 0, qJDD(3), 0, 0, 0.2e1 * t172 + t182, 0, t183 + 0.2e1 * t198, -t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, 0, t177, 0, 0, t132, -t171, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, 0, -t132, 0, t131, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, 0, 0, qJDD(3), 0, t171, -t131, 0, 0;];
m_new_reg  = t1;
