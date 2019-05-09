% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3RRP1
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
%   pkin=[a2,a3,d1,d2]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3RRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:31:26
% EndTime: 2019-05-04 18:31:27
% DurationCPUTime: 0.61s
% Computational Cost: add. (1190->100), mult. (1652->99), div. (0->0), fcn. (968->4), ass. (0->60)
t190 = (qJD(1) + qJD(2));
t188 = t190 ^ 2;
t189 = qJDD(1) + qJDD(2);
t191 = sin(qJ(2));
t193 = cos(qJ(2));
t171 = t193 * t188 + t191 * t189;
t174 = t191 * t188 - t193 * t189;
t192 = sin(qJ(1));
t194 = cos(qJ(1));
t152 = t194 * t171 - t192 * t174;
t160 = pkin(4) * t174 - t191 * g(3);
t163 = pkin(4) * t171 - t193 * g(3);
t230 = pkin(3) * t152 - t192 * t160 + t194 * t163;
t205 = t192 * t171 + t194 * t174;
t229 = pkin(3) * t205 + t194 * t160 + t192 * t163;
t180 = t194 * g(1) + t192 * g(2);
t196 = qJD(1) ^ 2;
t176 = -t196 * pkin(1) - t180;
t179 = t192 * g(1) - t194 * g(2);
t198 = qJDD(1) * pkin(1) + t179;
t157 = t191 * t176 - t193 * t198;
t158 = t193 * t176 + t191 * t198;
t207 = t191 * t157 + t193 * t158;
t141 = t193 * t157 - t191 * t158;
t213 = t194 * t141;
t228 = -t192 * t207 + t213;
t215 = t192 * t141;
t227 = t194 * t207 + t215;
t226 = pkin(1) * t171;
t210 = (2 * qJD(3) * t190) + t158;
t216 = t189 * qJ(3);
t147 = -t188 * pkin(2) + t210 + t216;
t184 = t189 * pkin(2);
t201 = -qJDD(3) - t157;
t151 = -t188 * qJ(3) - t184 - t201;
t132 = t191 * t147 - t193 * t151;
t208 = t193 * t147 + t191 * t151;
t225 = -t192 * t132 + t194 * t208;
t224 = t194 * t132 + t192 * t208;
t211 = -pkin(2) * t151 + qJ(3) * t147;
t204 = -t192 * t179 - t194 * t180;
t203 = 0.2e1 * t216 + t210;
t178 = t194 * qJDD(1) - t192 * t196;
t202 = -pkin(3) * t178 - t192 * g(3);
t199 = t194 * t179 - t192 * t180;
t197 = 0.2e1 * t184 + t201;
t195 = pkin(1) * g(3);
t177 = t192 * qJDD(1) + t194 * t196;
t167 = pkin(1) * t174;
t166 = -pkin(3) * t177 + t194 * g(3);
t149 = -t158 - t226;
t148 = -t157 - t167;
t144 = -t167 + t197;
t143 = t203 + t226;
t138 = pkin(1) * t141;
t137 = pkin(4) * t207 + t195;
t130 = -pkin(4) * t132 + (-pkin(2) * t191 + qJ(3) * t193) * g(3);
t129 = pkin(4) * t208 + t195 + (pkin(2) * t193 + qJ(3) * t191) * g(3);
t128 = pkin(1) * t132 + t211;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t178, 0, -t177, 0, t202, -t166, -t199, -pkin(3) * t199, 0, 0, -t205, 0, -t152, 0, t229, t230, t228, pkin(3) * t228 + pkin(4) * t213 - t192 * t137, 0, -t205, 0, 0, t152, 0, t229, -t224, -t230, -pkin(3) * t224 - t192 * t129 + t194 * t130; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t177, 0, t178, 0, t166, t202, t204, pkin(3) * t204, 0, 0, t152, 0, -t205, 0, -t230, t229, t227, pkin(3) * t227 + pkin(4) * t215 + t194 * t137, 0, t152, 0, 0, t205, 0, -t230, t225, -t229, pkin(3) * t225 + t194 * t129 + t192 * t130; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t179, t180, 0, 0, 0, 0, 0, 0, 0, t189, t148, t149, 0, -t138, 0, 0, 0, t189, 0, 0, t144, 0, t143, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t196, 0, 0, -g(3), -t179, 0, 0, 0, -t174, 0, -t171, 0, t160, t163, t141, pkin(4) * t141, 0, -t174, 0, 0, t171, 0, t160, -t132, -t163, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, qJDD(1), 0, g(3), 0, -t180, 0, 0, 0, t171, 0, -t174, 0, -t163, t160, t207, t137, 0, t171, 0, 0, t174, 0, -t163, t208, -t160, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t179, t180, 0, 0, 0, 0, 0, 0, 0, t189, t148, t149, 0, -t138, 0, 0, 0, t189, 0, 0, t144, 0, t143, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, -t188, 0, 0, -g(3), t157, 0, 0, t189, 0, 0, t188, 0, 0, t151, g(3), qJ(3) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, t189, 0, g(3), 0, t158, 0, 0, t188, 0, 0, -t189, 0, g(3), t147, 0, pkin(2) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, -t157, -t158, 0, 0, 0, 0, 0, t189, 0, 0, t197, 0, t203, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, 0, t188, 0, 0, t151, g(3), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0, 0, -t151, 0, t147, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, 0, 0, t189, 0, -g(3), -t147, 0, 0;];
m_new_reg  = t1;
