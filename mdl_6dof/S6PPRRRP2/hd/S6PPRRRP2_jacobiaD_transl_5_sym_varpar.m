% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:07
% EndTime: 2019-02-26 19:42:07
% DurationCPUTime: 0.49s
% Computational Cost: add. (599->84), mult. (1900->158), div. (0->0), fcn. (2202->14), ass. (0->66)
t461 = sin(pkin(12));
t462 = sin(pkin(11));
t447 = t462 * t461;
t465 = cos(pkin(12));
t466 = cos(pkin(11));
t453 = t466 * t465;
t468 = cos(pkin(6));
t431 = -t468 * t453 + t447;
t467 = cos(pkin(7));
t429 = t431 * t467;
t463 = sin(pkin(7));
t464 = sin(pkin(6));
t451 = t464 * t463;
t439 = t466 * t451;
t474 = t429 + t439;
t449 = t462 * t465;
t452 = t466 * t461;
t432 = t468 * t449 + t452;
t448 = t462 * t464;
t473 = t432 * t467 - t463 * t448;
t454 = t467 * t464;
t472 = t465 * t454 + t468 * t463;
t412 = t468 * t452 + t449;
t423 = sin(qJ(3));
t469 = cos(qJ(3));
t395 = t412 * t423 + t474 * t469;
t471 = t473 * t469;
t421 = sin(qJ(5));
t424 = cos(qJ(5));
t441 = qJD(5) * (r_i_i_C(1) * t421 + r_i_i_C(2) * t424);
t450 = t464 * t461;
t403 = t423 * t450 - t472 * t469;
t470 = pkin(10) + r_i_i_C(3);
t460 = qJD(3) * t423;
t459 = qJD(5) * t421;
t458 = qJD(5) * t424;
t457 = t412 * t469;
t396 = -t474 * t423 + t457;
t405 = t431 * t463 - t466 * t454;
t422 = sin(qJ(4));
t425 = cos(qJ(4));
t388 = t396 * t425 + t405 * t422;
t446 = -t396 * t422 + t405 * t425;
t413 = -t468 * t447 + t453;
t398 = t413 * t469 - t473 * t423;
t406 = t432 * t463 + t467 * t448;
t390 = t398 * t425 + t406 * t422;
t445 = -t398 * t422 + t406 * t425;
t404 = t472 * t423 + t469 * t450;
t411 = -t465 * t451 + t468 * t467;
t400 = t404 * t425 + t411 * t422;
t444 = -t404 * t422 + t411 * t425;
t443 = r_i_i_C(1) * t424 - r_i_i_C(2) * t421 + pkin(4);
t434 = -t470 * t422 - t443 * t425 - pkin(3);
t426 = t425 * t441 + (t443 * t422 - t470 * t425) * qJD(4);
t402 = t404 * qJD(3);
t401 = t403 * qJD(3);
t397 = t413 * t423 + t471;
t394 = t398 * qJD(3);
t393 = t471 * qJD(3) + t413 * t460;
t392 = -t439 * t460 + (-t423 * t429 + t457) * qJD(3);
t391 = t395 * qJD(3);
t386 = t444 * qJD(4) - t401 * t425;
t384 = t445 * qJD(4) - t393 * t425;
t382 = t446 * qJD(4) - t391 * t425;
t1 = [0, 0 (-t393 * t421 + t398 * t458) * r_i_i_C(1) + (-t393 * t424 - t398 * t459) * r_i_i_C(2) - t393 * pkin(9) + t434 * t394 + t426 * t397, t470 * t384 - t445 * t441 + t443 * (-t390 * qJD(4) + t393 * t422) (-t384 * t421 + t394 * t424) * r_i_i_C(1) + (-t384 * t424 - t394 * t421) * r_i_i_C(2) + ((-t390 * t424 - t397 * t421) * r_i_i_C(1) + (t390 * t421 - t397 * t424) * r_i_i_C(2)) * qJD(5), 0; 0, 0 (-t391 * t421 + t396 * t458) * r_i_i_C(1) + (-t391 * t424 - t396 * t459) * r_i_i_C(2) - t391 * pkin(9) + t434 * t392 + t426 * t395, t470 * t382 - t446 * t441 + t443 * (-t388 * qJD(4) + t391 * t422) (-t382 * t421 + t392 * t424) * r_i_i_C(1) + (-t382 * t424 - t392 * t421) * r_i_i_C(2) + ((-t388 * t424 - t395 * t421) * r_i_i_C(1) + (t388 * t421 - t395 * t424) * r_i_i_C(2)) * qJD(5), 0; 0, 0 (-t401 * t421 + t404 * t458) * r_i_i_C(1) + (-t401 * t424 - t404 * t459) * r_i_i_C(2) - t401 * pkin(9) + t434 * t402 + t426 * t403, t470 * t386 - t444 * t441 + t443 * (-t400 * qJD(4) + t401 * t422) (-t386 * t421 + t402 * t424) * r_i_i_C(1) + (-t386 * t424 - t402 * t421) * r_i_i_C(2) + ((-t400 * t424 - t403 * t421) * r_i_i_C(1) + (t400 * t421 - t403 * t424) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
