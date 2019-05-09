% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:03:32
% EndTime: 2019-05-04 19:03:32
% DurationCPUTime: 0.44s
% Computational Cost: add. (605->63), mult. (718->46), div. (0->0), fcn. (488->4), ass. (0->35)
t341 = qJD(2) + qJD(3);
t339 = t341 ^ 2;
t340 = qJDD(2) + qJDD(3);
t344 = sin(qJ(3));
t346 = cos(qJ(3));
t327 = t339 * t344 - t340 * t346;
t345 = sin(qJ(2));
t347 = cos(qJ(2));
t349 = -t339 * t346 - t340 * t344;
t352 = t327 * t345 + t347 * t349;
t316 = t327 * t347 - t345 * t349;
t343 = -g(2) + qJDD(1);
t331 = g(1) * t345 + t343 * t347;
t329 = qJDD(2) * pkin(2) + t331;
t332 = -g(1) * t347 + t343 * t345;
t348 = qJD(2) ^ 2;
t330 = -pkin(2) * t348 + t332;
t319 = t329 * t344 + t330 * t346;
t318 = t329 * t346 - t330 * t344;
t342 = g(3) - qJDD(4);
t334 = qJDD(2) * t347 - t345 * t348;
t333 = -qJDD(2) * t345 - t347 * t348;
t321 = -t331 * t345 + t332 * t347;
t320 = t331 * t347 + t332 * t345;
t313 = -pkin(3) * t339 + t319;
t312 = pkin(3) * t340 + t318;
t311 = -t318 * t344 + t319 * t346;
t310 = t318 * t346 + t319 * t344;
t309 = -t312 * t344 + t313 * t346;
t308 = t312 * t346 + t313 * t344;
t307 = -t310 * t345 + t311 * t347;
t306 = t310 * t347 + t311 * t345;
t305 = -t308 * t345 + t309 * t347;
t304 = t308 * t347 + t309 * t345;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t333, -t334, 0, t321, 0, 0, 0, 0, 0, 0, t352, t316, 0, t307, 0, 0, 0, 0, 0, 0, t352, t316, 0, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, 0, 0, 0, 0, 0, 0, t334, t333, 0, t320, 0, 0, 0, 0, 0, 0, -t316, t352, 0, t306, 0, 0, 0, 0, 0, 0, -t316, t352, 0, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t333, -t334, 0, t321, 0, 0, 0, 0, 0, 0, t352, t316, 0, t307, 0, 0, 0, 0, 0, 0, t352, t316, 0, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, 0, 0, 0, 0, 0, 0, t334, t333, 0, t320, 0, 0, 0, 0, 0, 0, -t316, t352, 0, t306, 0, 0, 0, 0, 0, 0, -t316, t352, 0, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t348, -qJDD(2), 0, t332, 0, 0, 0, 0, 0, 0, t349, t327, 0, t311, 0, 0, 0, 0, 0, 0, t349, t327, 0, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t348, 0, t331, 0, 0, 0, 0, 0, 0, -t327, t349, 0, t310, 0, 0, 0, 0, 0, 0, -t327, t349, 0, t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, -t340, 0, t319, 0, 0, 0, 0, 0, 0, -t339, -t340, 0, t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t339, 0, t318, 0, 0, 0, 0, 0, 0, t340, -t339, 0, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, -t340, 0, t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t339, 0, t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342;];
f_new_reg  = t1;
