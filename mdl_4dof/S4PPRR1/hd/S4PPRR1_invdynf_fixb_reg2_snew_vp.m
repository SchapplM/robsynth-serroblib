% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PPRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PPRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:47:51
% EndTime: 2019-05-04 18:47:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (724->68), mult. (1074->62), div. (0->0), fcn. (904->6), ass. (0->39)
t328 = qJD(3) + qJD(4);
t326 = t328 ^ 2;
t327 = qJDD(3) + qJDD(4);
t332 = sin(qJ(4));
t334 = cos(qJ(4));
t312 = t332 * t326 - t334 * t327;
t333 = sin(qJ(3));
t335 = cos(qJ(3));
t337 = -t334 * t326 - t332 * t327;
t304 = t335 * t312 - t333 * t337;
t330 = sin(pkin(6));
t331 = cos(pkin(6));
t342 = t333 * t312 + t335 * t337;
t346 = t330 * t304 - t331 * t342;
t345 = t331 * t304 + t330 * t342;
t323 = t330 * g(1) - t331 * g(2);
t320 = -qJDD(2) + t323;
t324 = -t331 * g(1) - t330 * g(2);
t309 = -t333 * t320 + t335 * t324;
t308 = -t335 * t320 - t333 * t324;
t336 = qJD(3) ^ 2;
t321 = t333 * qJDD(3) + t335 * t336;
t322 = -t335 * qJDD(3) + t333 * t336;
t339 = -t330 * t321 + t331 * t322;
t338 = t331 * t321 + t330 * t322;
t329 = g(3) - qJDD(1);
t318 = t331 * t324;
t317 = t330 * t324;
t307 = -t336 * pkin(3) + t309;
t306 = qJDD(3) * pkin(3) + t308;
t301 = -t333 * t308 + t335 * t309;
t300 = t335 * t308 + t333 * t309;
t299 = t332 * t306 + t334 * t307;
t298 = t334 * t306 - t332 * t307;
t297 = -t332 * t298 + t334 * t299;
t296 = t334 * t298 + t332 * t299;
t295 = -t333 * t296 + t335 * t297;
t294 = t335 * t296 + t333 * t297;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t330 * t323 + t318, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t330 * t320 + t318, 0, 0, 0, 0, 0, 0, -t338, t339, 0, t330 * t300 + t331 * t301, 0, 0, 0, 0, 0, 0, -t346, t345, 0, t330 * t294 + t331 * t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t331 * t323 + t317, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331 * t320 + t317, 0, 0, 0, 0, 0, 0, t339, t338, 0, -t331 * t300 + t330 * t301, 0, 0, 0, 0, 0, 0, t345, t346, 0, -t331 * t294 + t330 * t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, 0, 0, 0, 0, 0, 0, -t321, t322, 0, t301, 0, 0, 0, 0, 0, 0, t342, t304, 0, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, 0, 0, 0, 0, 0, 0, t322, t321, 0, -t300, 0, 0, 0, 0, 0, 0, t304, -t342, 0, -t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, 0, 0, 0, 0, 0, 0, -t321, t322, 0, t301, 0, 0, 0, 0, 0, 0, t342, t304, 0, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t320, 0, 0, 0, 0, 0, 0, -t322, -t321, 0, t300, 0, 0, 0, 0, 0, 0, -t304, t342, 0, t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, -qJDD(3), 0, t309, 0, 0, 0, 0, 0, 0, t337, t312, 0, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t336, 0, t308, 0, 0, 0, 0, 0, 0, -t312, t337, 0, t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t326, -t327, 0, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327, -t326, 0, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329;];
f_new_reg  = t1;
