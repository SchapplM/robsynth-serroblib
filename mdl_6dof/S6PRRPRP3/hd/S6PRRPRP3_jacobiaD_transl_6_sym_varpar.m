% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:25
% EndTime: 2019-02-26 20:02:25
% DurationCPUTime: 0.45s
% Computational Cost: add. (642->95), mult. (1461->159), div. (0->0), fcn. (1517->12), ass. (0->67)
t367 = cos(pkin(11)) * pkin(4) + pkin(3);
t375 = sin(qJ(3));
t377 = cos(qJ(3));
t370 = pkin(11) + qJ(5);
t368 = sin(t370);
t408 = qJD(6) * t368;
t421 = r_i_i_C(2) + pkin(9) + qJ(4);
t426 = (-t367 * t375 + t421 * t377) * qJD(3) + t375 * qJD(4) + t377 * t408;
t372 = sin(pkin(10));
t376 = sin(qJ(2));
t378 = cos(qJ(2));
t418 = cos(pkin(10));
t419 = cos(pkin(6));
t393 = t419 * t418;
t356 = t372 * t378 + t376 * t393;
t373 = sin(pkin(6));
t400 = t373 * t418;
t344 = t356 * t377 - t375 * t400;
t401 = t372 * t419;
t358 = -t376 * t401 + t418 * t378;
t413 = t373 * t377;
t360 = t419 * t375 + t376 * t413;
t369 = cos(t370);
t412 = t373 * t378;
t425 = -t360 * t369 + t368 * t412;
t410 = qJD(3) * t375;
t424 = (qJD(2) * t377 - qJD(5)) * t376 + t378 * t410;
t422 = pkin(5) + r_i_i_C(1);
t420 = r_i_i_C(3) + qJ(6);
t414 = t373 * t375;
t411 = qJD(2) * t376;
t409 = qJD(5) * t377;
t407 = t369 * qJD(6);
t405 = sin(pkin(11)) * pkin(4) + pkin(8);
t404 = t373 * t411;
t403 = qJD(2) * t412;
t389 = t378 * t393;
t351 = -qJD(2) * t389 + t372 * t411;
t355 = t372 * t376 - t389;
t395 = t355 * t409 - t351;
t357 = t418 * t376 + t378 * t401;
t353 = t357 * qJD(2);
t394 = t357 * t409 - t353;
t392 = t344 * t369 + t355 * t368;
t346 = t358 * t377 + t372 * t414;
t391 = t346 * t369 + t357 * t368;
t390 = (qJD(2) - t409) * t378;
t388 = -t358 * t375 + t372 * t413;
t386 = -t377 * t367 - t421 * t375 - pkin(2);
t385 = -t356 * t375 - t377 * t400;
t384 = -t376 * t414 + t419 * t377;
t383 = -t420 * t368 - t422 * t369 - t367;
t352 = t356 * qJD(2);
t382 = qJD(5) * t356 - t352 * t377 + t355 * t410;
t354 = t358 * qJD(2);
t381 = qJD(5) * t358 - t354 * t377 + t357 * t410;
t380 = t408 + (-t422 * t368 + t420 * t369) * qJD(5);
t348 = t384 * qJD(3) + t377 * t403;
t347 = t360 * qJD(3) + t375 * t403;
t342 = t388 * qJD(3) - t353 * t377;
t341 = t346 * qJD(3) - t353 * t375;
t340 = t385 * qJD(3) - t351 * t377;
t339 = t344 * qJD(3) - t351 * t375;
t335 = -t425 * qJD(5) + t348 * t368 - t369 * t404;
t329 = t391 * qJD(5) + t342 * t368 - t354 * t369;
t327 = t392 * qJD(5) + t340 * t368 - t352 * t369;
t1 = [0, -t358 * t407 - t405 * t353 + t422 * (t394 * t368 + t381 * t369) + t420 * (t381 * t368 - t394 * t369) - t426 * t357 + t386 * t354, t346 * qJD(4) + t383 * t341 + t421 * t342 + t380 * t388, t341, t391 * qJD(6) + t420 * (t342 * t369 + t354 * t368 + (-t346 * t368 + t357 * t369) * qJD(5)) - t422 * t329, t329; 0, -t356 * t407 - t405 * t351 + t422 * (t395 * t368 + t382 * t369) + t420 * (t382 * t368 - t395 * t369) - t426 * t355 + t386 * t352, t344 * qJD(4) + t383 * t339 + t421 * t340 + t380 * t385, t339, t392 * qJD(6) + t420 * (t340 * t369 + t352 * t368 + (-t344 * t368 + t355 * t369) * qJD(5)) - t422 * t327, t327; 0 (t422 * (t368 * t390 - t424 * t369) - t420 * (t424 * t368 + t369 * t390) - t376 * t407 + t426 * t378 + (t386 * t376 + t405 * t378) * qJD(2)) * t373, t360 * qJD(4) + t383 * t347 + t421 * t348 + t380 * t384, t347, -t425 * qJD(6) + t420 * (t368 * t404 + t348 * t369 + (-t360 * t368 - t369 * t412) * qJD(5)) - t422 * t335, t335;];
JaD_transl  = t1;
