% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRPP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:54:39
% EndTime: 2019-05-04 18:54:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (539->70), mult. (778->49), div. (0->0), fcn. (488->4), ass. (0->32)
t327 = sin(pkin(5));
t328 = cos(pkin(5));
t331 = qJD(2) ^ 2;
t316 = t327 * qJDD(2) + t328 * t331;
t317 = t328 * qJDD(2) - t327 * t331;
t329 = sin(qJ(2));
t330 = cos(qJ(2));
t334 = -t329 * t316 + t330 * t317;
t302 = t330 * t316 + t329 * t317;
t325 = -g(2) + qJDD(1);
t315 = -t330 * g(1) + t329 * t325;
t311 = -t331 * pkin(2) + t315;
t314 = t329 * g(1) + t330 * t325;
t332 = qJDD(2) * pkin(2) + t314;
t299 = t328 * t311 + t327 * t332;
t298 = -t327 * t311 + t328 * t332;
t324 = g(3) - qJDD(3);
t319 = t330 * qJDD(2) - t329 * t331;
t318 = -t329 * qJDD(2) - t330 * t331;
t301 = -t329 * t314 + t330 * t315;
t300 = t330 * t314 + t329 * t315;
t297 = -qJDD(2) * pkin(3) - t331 * qJ(4) + qJDD(4) - t298;
t296 = -t331 * pkin(3) + qJDD(2) * qJ(4) + (2 * qJD(4) * qJD(2)) + t299;
t295 = -t327 * t298 + t328 * t299;
t294 = t328 * t298 + t327 * t299;
t293 = t328 * t296 + t327 * t297;
t292 = t327 * t296 - t328 * t297;
t291 = -t329 * t294 + t330 * t295;
t290 = t330 * t294 + t329 * t295;
t289 = -t329 * t292 + t330 * t293;
t288 = t330 * t292 + t329 * t293;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t318, -t319, 0, t301, 0, 0, 0, 0, 0, 0, -t302, -t334, 0, t291, 0, 0, 0, 0, 0, 0, -t302, 0, t334, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, 0, 0, 0, 0, 0, 0, t319, t318, 0, t300, 0, 0, 0, 0, 0, 0, t334, -t302, 0, t290, 0, 0, 0, 0, 0, 0, t334, 0, t302, t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t318, -t319, 0, t301, 0, 0, 0, 0, 0, 0, -t302, -t334, 0, t291, 0, 0, 0, 0, 0, 0, -t302, 0, t334, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, 0, 0, 0, 0, 0, 0, t319, t318, 0, t300, 0, 0, 0, 0, 0, 0, t334, -t302, 0, t290, 0, 0, 0, 0, 0, 0, t334, 0, t302, t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t331, -qJDD(2), 0, t315, 0, 0, 0, 0, 0, 0, -t316, -t317, 0, t295, 0, 0, 0, 0, 0, 0, -t316, 0, t317, t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t331, 0, t314, 0, 0, 0, 0, 0, 0, t317, -t316, 0, t294, 0, 0, 0, 0, 0, 0, t317, 0, t316, t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t331, -qJDD(2), 0, t299, 0, 0, 0, 0, 0, 0, -t331, 0, qJDD(2), t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t331, 0, t298, 0, 0, 0, 0, 0, 0, qJDD(2), 0, t331, -t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t331, 0, qJDD(2), t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t331, t297;];
f_new_reg  = t1;
