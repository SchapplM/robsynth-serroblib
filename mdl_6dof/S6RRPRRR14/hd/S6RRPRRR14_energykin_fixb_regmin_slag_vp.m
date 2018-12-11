% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:32:32
% EndTime: 2018-12-10 18:32:32
% DurationCPUTime: 0.47s
% Computational Cost: add. (2027->72), mult. (5772->148), div. (0->0), fcn. (4990->16), ass. (0->66)
t381 = cos(pkin(6)) * qJD(1);
t351 = qJD(2) + t381;
t353 = sin(pkin(14));
t358 = cos(pkin(7));
t367 = cos(qJ(2));
t356 = sin(pkin(6));
t382 = qJD(1) * t356;
t376 = t367 * t382;
t374 = t358 * t376;
t355 = sin(pkin(7));
t391 = cos(pkin(14));
t375 = t355 * t391;
t363 = sin(qJ(2));
t377 = t363 * t382;
t337 = -t351 * t375 + t353 * t377 - t391 * t374;
t354 = sin(pkin(8));
t357 = cos(pkin(8));
t371 = -t358 * t351 + t355 * t376;
t395 = t357 * t337 + t354 * t371;
t380 = pkin(1) * t381;
t383 = pkin(10) * t376 + t363 * t380;
t336 = (t351 * t355 + t374) * qJ(3) + t383;
t350 = t367 * t380;
t339 = t351 * pkin(2) + t350 + (-qJ(3) * t358 - pkin(10)) * t377;
t344 = (-qJ(3) * t355 * t363 - pkin(2) * t367 - pkin(1)) * t382;
t322 = t358 * t391 * t339 - t353 * t336 + t344 * t375;
t386 = t353 * t358;
t387 = t353 * t355;
t338 = t351 * t387 + (t391 * t363 + t367 * t386) * t382;
t392 = pkin(11) * t338;
t318 = -t371 * pkin(3) - t357 * t392 + t322;
t328 = -t355 * t339 + t358 * t344 + qJD(3);
t321 = t337 * pkin(3) - t354 * t392 + t328;
t394 = t318 * t357 + t321 * t354;
t323 = t391 * t336 + t339 * t386 + t344 * t387;
t315 = -pkin(11) * t395 + t323;
t362 = sin(qJ(4));
t366 = cos(qJ(4));
t393 = -t362 * t315 + t394 * t366;
t368 = qJD(1) ^ 2;
t388 = t356 ^ 2 * t368;
t330 = -t354 * t337 + t357 * t371 - qJD(4);
t378 = t366 * t315 + t394 * t362;
t305 = -pkin(12) * t330 + t378;
t308 = -t354 * t318 + t357 * t321;
t326 = t362 * t338 + t395 * t366;
t327 = t366 * t338 - t362 * t395;
t307 = t326 * pkin(4) - t327 * pkin(12) + t308;
t361 = sin(qJ(5));
t365 = cos(qJ(5));
t384 = t365 * t305 + t361 * t307;
t379 = t367 * t388;
t316 = t361 * t327 + t365 * t330;
t373 = -t361 * t305 + t365 * t307;
t304 = t330 * pkin(4) - t393;
t364 = cos(qJ(6));
t360 = sin(qJ(6));
t325 = qJD(5) + t326;
t317 = t365 * t327 - t361 * t330;
t314 = qJD(6) + t316;
t310 = t364 * t317 + t360 * t325;
t309 = t360 * t317 - t364 * t325;
t302 = t316 * pkin(5) - t317 * pkin(13) + t304;
t301 = pkin(13) * t325 + t384;
t300 = -t325 * pkin(5) - t373;
t1 = [t368 / 0.2e1, 0, 0, t363 ^ 2 * t388 / 0.2e1, t363 * t379, t351 * t377, t351 * t376, t351 ^ 2 / 0.2e1 (-pkin(10) * t377 + t350) * t351 + pkin(1) * t379, -pkin(1) * t363 * t388 - t383 * t351, -t322 * t371 + t328 * t337, t323 * t371 + t328 * t338, -t322 * t338 - t323 * t337, t323 ^ 2 / 0.2e1 + t322 ^ 2 / 0.2e1 + t328 ^ 2 / 0.2e1, t327 ^ 2 / 0.2e1, -t326 * t327, -t330 * t327, t330 * t326, t330 ^ 2 / 0.2e1, t308 * t326 - t330 * t393, t308 * t327 + t378 * t330, t317 ^ 2 / 0.2e1, -t316 * t317, t325 * t317, -t316 * t325, t325 ^ 2 / 0.2e1, t304 * t316 + t373 * t325, t304 * t317 - t384 * t325, t310 ^ 2 / 0.2e1, -t310 * t309, t314 * t310, -t314 * t309, t314 ^ 2 / 0.2e1 (-t360 * t301 + t364 * t302) * t314 + t300 * t309, t300 * t310 - (t364 * t301 + t360 * t302) * t314;];
T_reg  = t1;
