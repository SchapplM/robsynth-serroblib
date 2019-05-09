% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4RPRP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynf_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:15:26
% EndTime: 2019-05-04 19:15:26
% DurationCPUTime: 0.46s
% Computational Cost: add. (673->79), mult. (932->56), div. (0->0), fcn. (400->4), ass. (0->33)
t342 = sin(qJ(1));
t344 = cos(qJ(1));
t338 = -qJD(1) + qJD(3);
t336 = t338 ^ 2;
t337 = qJDD(1) - qJDD(3);
t341 = sin(qJ(3));
t343 = cos(qJ(3));
t348 = t341 * t336 + t343 * t337;
t349 = -t343 * t336 + t341 * t337;
t351 = t342 * t349 + t344 * t348;
t313 = t342 * t348 - t344 * t349;
t350 = -pkin(1) - pkin(2);
t345 = qJD(1) ^ 2;
t331 = -t344 * g(1) - t342 * g(2);
t347 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t331;
t320 = t350 * t345 + t347;
t330 = t342 * g(1) - t344 * g(2);
t346 = -t345 * qJ(2) + qJDD(2) - t330;
t321 = t350 * qJDD(1) + t346;
t312 = t343 * t320 + t341 * t321;
t311 = -t341 * t320 + t343 * t321;
t340 = g(3) + qJDD(4);
t329 = t344 * qJDD(1) - t342 * t345;
t328 = t342 * qJDD(1) + t344 * t345;
t323 = qJDD(1) * pkin(1) - t346;
t322 = -t345 * pkin(1) + t347;
t310 = -t336 * pkin(3) + t312;
t309 = -t337 * pkin(3) + t311;
t308 = -t341 * t311 + t343 * t312;
t307 = t343 * t311 + t341 * t312;
t306 = -t341 * t309 + t343 * t310;
t305 = t343 * t309 + t341 * t310;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t328, -t329, 0, -t342 * t330 + t344 * t331, 0, 0, 0, 0, 0, 0, -t328, 0, t329, t344 * t322 - t342 * t323, 0, 0, 0, 0, 0, 0, -t313, t351, 0, t342 * t307 + t344 * t308, 0, 0, 0, 0, 0, 0, -t313, t351, 0, t342 * t305 + t344 * t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t329, -t328, 0, t344 * t330 + t342 * t331, 0, 0, 0, 0, 0, 0, t329, 0, t328, t342 * t322 + t344 * t323, 0, 0, 0, 0, 0, 0, t351, t313, 0, -t344 * t307 + t342 * t308, 0, 0, 0, 0, 0, 0, t351, t313, 0, -t344 * t305 + t342 * t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, -qJDD(1), 0, t331, 0, 0, 0, 0, 0, 0, -t345, 0, qJDD(1), t322, 0, 0, 0, 0, 0, 0, t349, t348, 0, t308, 0, 0, 0, 0, 0, 0, t349, t348, 0, t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t345, 0, t330, 0, 0, 0, 0, 0, 0, qJDD(1), 0, t345, t323, 0, 0, 0, 0, 0, 0, t348, -t349, 0, -t307, 0, 0, 0, 0, 0, 0, t348, -t349, 0, -t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, 0, qJDD(1), t322, 0, 0, 0, 0, 0, 0, t349, t348, 0, t308, 0, 0, 0, 0, 0, 0, t349, t348, 0, t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t345, -t323, 0, 0, 0, 0, 0, 0, -t348, t349, 0, t307, 0, 0, 0, 0, 0, 0, -t348, t349, 0, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, t337, 0, t312, 0, 0, 0, 0, 0, 0, -t336, t337, 0, t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, -t336, 0, t311, 0, 0, 0, 0, 0, 0, -t337, -t336, 0, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, t337, 0, t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, -t336, 0, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340;];
f_new_reg  = t1;
