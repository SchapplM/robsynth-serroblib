% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S3RPP1
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
%
% Output:
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3RPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:27:39
% EndTime: 2019-05-04 18:27:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (142->58), mult. (205->31), div. (0->0), fcn. (84->2), ass. (0->15)
t158 = pkin(1) + qJ(3);
t153 = sin(qJ(1));
t154 = cos(qJ(1));
t146 = t153 * g(1) - t154 * g(2);
t147 = -t154 * g(1) - t153 * g(2);
t155 = qJD(1) ^ 2;
t157 = t155 * qJ(2) - qJDD(2) + t146;
t156 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t147;
t145 = t153 * qJDD(1) + t154 * t155;
t144 = t154 * qJDD(1) - t153 * t155;
t143 = qJDD(1) * pkin(1) + t157;
t142 = t155 * pkin(1) + t156;
t141 = t158 * t155 - qJDD(3) + t156;
t140 = (2 * qJD(3) * qJD(1)) + t158 * qJDD(1) + t157;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t145, -t144, 0, -t153 * t146 + t154 * t147, 0, 0, 0, 0, 0, 0, 0, t145, t144, -t154 * t142 - t153 * t143, 0, 0, 0, 0, 0, 0, 0, t144, -t145, -t153 * t140 - t154 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t144, -t145, 0, t154 * t146 + t153 * t147, 0, 0, 0, 0, 0, 0, 0, -t144, t145, -t153 * t142 + t154 * t143, 0, 0, 0, 0, 0, 0, 0, t145, t144, t154 * t140 - t153 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -qJDD(1), 0, t147, 0, 0, 0, 0, 0, 0, 0, t155, qJDD(1), -t142, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t155, -t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t155, 0, t146, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t155, t143, 0, 0, 0, 0, 0, 0, 0, t155, qJDD(1), t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -qJDD(1), t142, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), t155, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t155, -t143, 0, 0, 0, 0, 0, 0, 0, -t155, -qJDD(1), -t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -qJDD(1), -t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t155, -t141;];
f_new_reg  = t1;
