% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S3RRR1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
%
% Output:
% f_new_reg [(3*4)x(4*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S3RRR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynf_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynf_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:32:50
% EndTime: 2019-05-04 18:32:50
% DurationCPUTime: 0.34s
% Computational Cost: add. (789->65), mult. (1070->67), div. (0->0), fcn. (708->6), ass. (0->42)
t342 = qJD(1) + qJD(2);
t338 = qJD(3) + t342;
t336 = t338 ^ 2;
t341 = qJDD(1) + qJDD(2);
t337 = qJDD(3) + t341;
t343 = sin(qJ(3));
t346 = cos(qJ(3));
t321 = t336 * t343 - t337 * t346;
t344 = sin(qJ(2));
t347 = cos(qJ(2));
t351 = -t336 * t346 - t337 * t343;
t313 = t321 * t347 - t344 * t351;
t345 = sin(qJ(1));
t348 = cos(qJ(1));
t358 = t321 * t344 + t347 * t351;
t362 = t313 * t345 + t348 * t358;
t361 = t313 * t348 - t345 * t358;
t340 = t342 ^ 2;
t328 = t340 * t344 - t341 * t347;
t350 = -t340 * t347 - t341 * t344;
t357 = t328 * t345 + t348 * t350;
t356 = t328 * t348 - t345 * t350;
t334 = g(1) * t345 - g(2) * t348;
t330 = qJDD(1) * pkin(1) + t334;
t335 = -g(1) * t348 - g(2) * t345;
t349 = qJD(1) ^ 2;
t331 = -pkin(1) * t349 + t335;
t318 = t330 * t344 + t331 * t347;
t317 = t330 * t347 - t331 * t344;
t333 = -qJDD(1) * t345 - t348 * t349;
t332 = qJDD(1) * t348 - t345 * t349;
t316 = -pkin(2) * t340 + t318;
t315 = pkin(2) * t341 + t317;
t310 = -t317 * t344 + t318 * t347;
t309 = t317 * t347 + t318 * t344;
t308 = t315 * t343 + t316 * t346;
t307 = t315 * t346 - t316 * t343;
t306 = -t307 * t343 + t308 * t346;
t305 = t307 * t346 + t308 * t343;
t304 = -t305 * t344 + t306 * t347;
t303 = t305 * t347 + t306 * t344;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t333, -t332, 0, -t334 * t345 + t335 * t348, 0, 0, 0, 0, 0, 0, t357, t356, 0, -t309 * t345 + t310 * t348, 0, 0, 0, 0, 0, 0, t362, t361, 0, -t303 * t345 + t304 * t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t332, t333, 0, t334 * t348 + t335 * t345, 0, 0, 0, 0, 0, 0, -t356, t357, 0, t309 * t348 + t310 * t345, 0, 0, 0, 0, 0, 0, -t361, t362, 0, t303 * t348 + t304 * t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349, -qJDD(1), 0, t335, 0, 0, 0, 0, 0, 0, t350, t328, 0, t310, 0, 0, 0, 0, 0, 0, t358, t313, 0, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t349, 0, t334, 0, 0, 0, 0, 0, 0, -t328, t350, 0, t309, 0, 0, 0, 0, 0, 0, -t313, t358, 0, t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, -t341, 0, t318, 0, 0, 0, 0, 0, 0, t351, t321, 0, t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, -t340, 0, t317, 0, 0, 0, 0, 0, 0, -t321, t351, 0, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, -t337, 0, t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, -t336, 0, t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
f_new_reg  = t1;
