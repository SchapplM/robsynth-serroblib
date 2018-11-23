% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:16:23
% EndTime: 2018-11-23 11:16:23
% DurationCPUTime: 0.45s
% Computational Cost: add. (2094->71), mult. (5637->147), div. (0->0), fcn. (4907->16), ass. (0->66)
t373 = cos(pkin(6)) * qJD(1);
t343 = qJD(2) + t373;
t354 = sin(qJ(3));
t349 = cos(pkin(7));
t359 = cos(qJ(2));
t347 = sin(pkin(6));
t374 = qJD(1) * t347;
t366 = t359 * t374;
t365 = t349 * t366;
t355 = sin(qJ(2));
t367 = t355 * t374;
t346 = sin(pkin(7));
t385 = cos(qJ(3));
t368 = t346 * t385;
t328 = -t343 * t368 + t354 * t367 - t385 * t365;
t336 = -t349 * t343 + t346 * t366 - qJD(3);
t345 = sin(pkin(8));
t348 = cos(pkin(8));
t388 = t328 * t348 + t336 * t345;
t372 = pkin(1) * t373;
t375 = pkin(10) * t366 + t355 * t372;
t326 = (t343 * t346 + t365) * pkin(11) + t375;
t342 = t359 * t372;
t327 = t343 * pkin(2) + t342 + (-pkin(11) * t349 - pkin(10)) * t367;
t335 = (-pkin(11) * t346 * t355 - pkin(2) * t359 - pkin(1)) * t374;
t364 = t349 * t385 * t327 - t354 * t326 + t335 * t368;
t377 = t349 * t354;
t378 = t346 * t354;
t329 = t343 * t378 + (t385 * t355 + t359 * t377) * t374;
t384 = pkin(12) * t329;
t310 = -t336 * pkin(3) - t348 * t384 + t364;
t319 = -t346 * t327 + t349 * t335;
t313 = t328 * pkin(3) - t345 * t384 + t319;
t387 = t310 * t348 + t313 * t345;
t369 = t385 * t326 + t327 * t377 + t335 * t378;
t309 = -pkin(12) * t388 + t369;
t353 = sin(qJ(4));
t358 = cos(qJ(4));
t386 = -t353 * t309 + t387 * t358;
t360 = qJD(1) ^ 2;
t379 = t347 ^ 2 * t360;
t320 = -t345 * t328 + t348 * t336 - qJD(4);
t370 = t358 * t309 + t387 * t353;
t297 = -t320 * pkin(13) + t370;
t300 = -t345 * t310 + t348 * t313;
t316 = t353 * t329 + t388 * t358;
t317 = t358 * t329 - t353 * t388;
t299 = t316 * pkin(4) - t317 * pkin(13) + t300;
t352 = sin(qJ(5));
t357 = cos(qJ(5));
t376 = t357 * t297 + t352 * t299;
t371 = t359 * t379;
t307 = t352 * t317 + t357 * t320;
t363 = -t352 * t297 + t357 * t299;
t296 = t320 * pkin(4) - t386;
t356 = cos(qJ(6));
t351 = sin(qJ(6));
t315 = qJD(5) + t316;
t308 = t357 * t317 - t352 * t320;
t306 = qJD(6) + t307;
t302 = t356 * t308 + t351 * t315;
t301 = t351 * t308 - t356 * t315;
t294 = t307 * pkin(5) - t308 * pkin(14) + t296;
t293 = t315 * pkin(14) + t376;
t292 = -t315 * pkin(5) - t363;
t1 = [t360 / 0.2e1, 0, 0, t355 ^ 2 * t379 / 0.2e1, t355 * t371, t343 * t367, t343 * t366, t343 ^ 2 / 0.2e1 (-pkin(10) * t367 + t342) * t343 + pkin(1) * t371, -pkin(1) * t355 * t379 - t375 * t343, t329 ^ 2 / 0.2e1, -t328 * t329, -t336 * t329, t328 * t336, t336 ^ 2 / 0.2e1, t319 * t328 - t364 * t336, t319 * t329 + t369 * t336, t317 ^ 2 / 0.2e1, -t316 * t317, -t320 * t317, t320 * t316, t320 ^ 2 / 0.2e1, t300 * t316 - t386 * t320, t300 * t317 + t370 * t320, t308 ^ 2 / 0.2e1, -t308 * t307, t308 * t315, -t307 * t315, t315 ^ 2 / 0.2e1, t296 * t307 + t363 * t315, t296 * t308 - t376 * t315, t302 ^ 2 / 0.2e1, -t302 * t301, t302 * t306, -t301 * t306, t306 ^ 2 / 0.2e1 (-t351 * t293 + t356 * t294) * t306 + t292 * t301 -(t356 * t293 + t351 * t294) * t306 + t292 * t302;];
T_reg  = t1;
