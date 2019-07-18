% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output:
% f_new_reg [(3*5)x(5*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S4PRRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:28
% EndTime: 2019-07-18 13:27:29
% DurationCPUTime: 0.56s
% Computational Cost: add. (1227->72), mult. (1678->67), div. (0->0), fcn. (1120->6), ass. (0->49)
t352 = qJD(2) + qJD(3);
t346 = qJD(4) + t352;
t344 = t346 ^ 2;
t351 = qJDD(2) + qJDD(3);
t345 = qJDD(4) + t351;
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t329 = t354 * t344 - t357 * t345;
t355 = sin(qJ(3));
t358 = cos(qJ(3));
t362 = -t357 * t344 - t354 * t345;
t315 = t358 * t329 - t355 * t362;
t356 = sin(qJ(2));
t359 = cos(qJ(2));
t368 = t355 * t329 + t358 * t362;
t371 = t356 * t315 + t359 * t368;
t307 = t359 * t315 - t356 * t368;
t350 = t352 ^ 2;
t336 = t355 * t350 - t358 * t351;
t361 = -t358 * t350 - t355 * t351;
t367 = t356 * t336 + t359 * t361;
t321 = t359 * t336 - t356 * t361;
t342 = t359 * g(1) + t356 * g(3);
t338 = qJDD(2) * pkin(1) + t342;
t343 = t356 * g(1) - t359 * g(3);
t360 = qJD(2) ^ 2;
t339 = -t360 * pkin(1) + t343;
t324 = t355 * t338 + t358 * t339;
t323 = t358 * t338 - t355 * t339;
t353 = g(2) + qJDD(1);
t341 = -t359 * qJDD(2) + t356 * t360;
t340 = t356 * qJDD(2) + t359 * t360;
t326 = -t356 * t342 + t359 * t343;
t325 = t359 * t342 + t356 * t343;
t318 = -t350 * pkin(2) + t324;
t317 = t351 * pkin(2) + t323;
t312 = -t355 * t323 + t358 * t324;
t311 = t358 * t323 + t355 * t324;
t310 = t354 * t317 + t357 * t318;
t309 = t357 * t317 - t354 * t318;
t304 = -t356 * t311 + t359 * t312;
t303 = t359 * t311 + t356 * t312;
t302 = -t354 * t309 + t357 * t310;
t301 = t357 * t309 + t354 * t310;
t300 = -t355 * t301 + t358 * t302;
t299 = t358 * t301 + t355 * t302;
t298 = -t356 * t299 + t359 * t300;
t297 = t359 * t299 + t356 * t300;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t341, t340, 0, -t325, 0, 0, 0, 0, 0, 0, t321, -t367, 0, -t303, 0, 0, 0, 0, 0, 0, t307, -t371, 0, -t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t340, t341, 0, t326, 0, 0, 0, 0, 0, 0, t367, t321, 0, t304, 0, 0, 0, 0, 0, 0, t371, t307, 0, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t340, t341, 0, t326, 0, 0, 0, 0, 0, 0, t367, t321, 0, t304, 0, 0, 0, 0, 0, 0, t371, t307, 0, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), 0, 0, 0, 0, 0, 0, -t341, -t340, 0, t325, 0, 0, 0, 0, 0, 0, -t321, t367, 0, t303, 0, 0, 0, 0, 0, 0, -t307, t371, 0, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t360, -qJDD(2), 0, t343, 0, 0, 0, 0, 0, 0, t361, t336, 0, t312, 0, 0, 0, 0, 0, 0, t368, t315, 0, t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t360, 0, t342, 0, 0, 0, 0, 0, 0, -t336, t361, 0, t311, 0, 0, 0, 0, 0, 0, -t315, t368, 0, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, -t351, 0, t324, 0, 0, 0, 0, 0, 0, t362, t329, 0, t302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t351, -t350, 0, t323, 0, 0, 0, 0, 0, 0, -t329, t362, 0, t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t344, -t345, 0, t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345, -t344, 0, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t353;];
f_new_reg  = t1;
