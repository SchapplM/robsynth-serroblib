% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3RRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:31:27
% EndTime: 2019-05-04 18:31:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (378->66), mult. (490->50), div. (0->0), fcn. (300->4), ass. (0->27)
t258 = (qJD(1) + qJD(2));
t256 = t258 ^ 2;
t257 = qJDD(1) + qJDD(2);
t259 = sin(qJ(2));
t261 = cos(qJ(2));
t245 = t261 * t256 + t259 * t257;
t248 = t259 * t256 - t261 * t257;
t260 = sin(qJ(1));
t262 = cos(qJ(1));
t265 = t262 * t245 - t260 * t248;
t266 = t260 * t245 + t262 * t248;
t253 = -t262 * g(1) - t260 * g(2);
t263 = qJD(1) ^ 2;
t249 = -t263 * pkin(1) + t253;
t252 = t260 * g(1) - t262 * g(2);
t264 = qJDD(1) * pkin(1) + t252;
t240 = t261 * t249 + t259 * t264;
t239 = -t259 * t249 + t261 * t264;
t251 = -t260 * qJDD(1) - t262 * t263;
t250 = t262 * qJDD(1) - t260 * t263;
t236 = -t257 * pkin(2) - t256 * qJ(3) + qJDD(3) - t239;
t235 = -t256 * pkin(2) + t257 * qJ(3) + (2 * qJD(3) * t258) + t240;
t234 = -t259 * t239 + t261 * t240;
t233 = t261 * t239 + t259 * t240;
t232 = t261 * t235 + t259 * t236;
t231 = t259 * t235 - t261 * t236;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t251, -t250, 0, -t260 * t252 + t262 * t253, 0, 0, 0, 0, 0, 0, -t265, t266, 0, -t260 * t233 + t262 * t234, 0, 0, 0, 0, 0, 0, -t265, 0, -t266, -t260 * t231 + t262 * t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t250, t251, 0, t262 * t252 + t260 * t253, 0, 0, 0, 0, 0, 0, -t266, -t265, 0, t262 * t233 + t260 * t234, 0, 0, 0, 0, 0, 0, -t266, 0, t265, t262 * t231 + t260 * t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, -qJDD(1), 0, t253, 0, 0, 0, 0, 0, 0, -t245, t248, 0, t234, 0, 0, 0, 0, 0, 0, -t245, 0, -t248, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t263, 0, t252, 0, 0, 0, 0, 0, 0, -t248, -t245, 0, t233, 0, 0, 0, 0, 0, 0, -t248, 0, t245, t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, -t257, 0, t240, 0, 0, 0, 0, 0, 0, -t256, 0, t257, t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, -t256, 0, t239, 0, 0, 0, 0, 0, 0, t257, 0, t256, -t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, 0, t257, t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, 0, -t256, t236;];
f_new_reg  = t1;
