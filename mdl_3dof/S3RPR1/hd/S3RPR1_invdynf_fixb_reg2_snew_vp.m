% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3RPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:57
% EndTime: 2019-05-04 18:28:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (334->59), mult. (482->46), div. (0->0), fcn. (220->4), ass. (0->28)
t214 = sin(qJ(1));
t216 = cos(qJ(1));
t211 = -qJD(1) + qJD(3);
t209 = t211 ^ 2;
t210 = qJDD(1) - qJDD(3);
t213 = sin(qJ(3));
t215 = cos(qJ(3));
t220 = t213 * t209 + t215 * t210;
t221 = -t215 * t209 + t213 * t210;
t224 = t214 * t221 + t216 * t220;
t223 = t214 * t220 - t216 * t221;
t222 = -pkin(2) - pkin(1);
t204 = t214 * g(1) - t216 * g(2);
t205 = -t216 * g(1) - t214 * g(2);
t219 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t205;
t217 = qJD(1) ^ 2;
t218 = -t217 * qJ(2) + qJDD(2) - t204;
t203 = t216 * qJDD(1) - t214 * t217;
t202 = t214 * qJDD(1) + t216 * t217;
t197 = qJDD(1) * pkin(1) - t218;
t196 = -t217 * pkin(1) + t219;
t195 = t222 * qJDD(1) + t218;
t194 = t222 * t217 + t219;
t193 = t215 * t194 + t213 * t195;
t192 = -t213 * t194 + t215 * t195;
t191 = -t213 * t192 + t215 * t193;
t190 = t215 * t192 + t213 * t193;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t202, -t203, 0, -t214 * t204 + t216 * t205, 0, 0, 0, 0, 0, 0, -t202, 0, t203, t216 * t196 - t214 * t197, 0, 0, 0, 0, 0, 0, -t223, t224, 0, t214 * t190 + t216 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t203, -t202, 0, t216 * t204 + t214 * t205, 0, 0, 0, 0, 0, 0, t203, 0, t202, t214 * t196 + t216 * t197, 0, 0, 0, 0, 0, 0, t224, t223, 0, -t216 * t190 + t214 * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, -qJDD(1), 0, t205, 0, 0, 0, 0, 0, 0, -t217, 0, qJDD(1), t196, 0, 0, 0, 0, 0, 0, t221, t220, 0, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t217, 0, t204, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t217, t197, 0, 0, 0, 0, 0, 0, t220, -t221, 0, -t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, 0, qJDD(1), t196, 0, 0, 0, 0, 0, 0, t221, t220, 0, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t217, -t197, 0, 0, 0, 0, 0, 0, -t220, t221, 0, t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t210, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, -t209, 0, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
f_new_reg  = t1;
