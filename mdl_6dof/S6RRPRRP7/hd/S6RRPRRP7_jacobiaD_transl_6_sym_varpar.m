% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:28
% EndTime: 2019-02-26 21:49:28
% DurationCPUTime: 0.53s
% Computational Cost: add. (570->71), mult. (1731->112), div. (0->0), fcn. (1664->8), ass. (0->53)
t421 = cos(qJ(2));
t465 = sin(qJ(4));
t466 = sin(qJ(2));
t467 = cos(qJ(4));
t399 = -t421 * t465 + t466 * t467;
t418 = sin(qJ(5));
t420 = cos(qJ(5));
t464 = r_i_i_C(3) + qJ(6);
t469 = -r_i_i_C(1) - pkin(5);
t472 = t469 * t418 + t464 * t420;
t424 = qJD(5) * t472 + qJD(6) * t418;
t487 = t399 * t424;
t398 = t421 * t467 + t466 * t465;
t483 = qJD(2) - qJD(4);
t392 = t483 * t398;
t419 = sin(qJ(1));
t422 = cos(qJ(1));
t475 = t399 * qJD(1);
t390 = -t392 * t419 - t475 * t422;
t393 = t483 * t399;
t429 = qJD(1) * t398;
t391 = -t393 * t419 + t422 * t429;
t435 = -t464 * t418 + t469 * t420;
t432 = -pkin(4) + t435;
t468 = -r_i_i_C(2) - pkin(9);
t486 = t432 * t390 - t468 * t391 + t419 * t487;
t388 = t393 * t422 + t419 * t429;
t397 = t398 * t422;
t389 = t483 * t397 - t475 * t419;
t485 = t468 * t388 - t432 * t389 + t422 * t487;
t484 = -t468 * t392 - t432 * t393 - t424 * t398;
t479 = pkin(2) + pkin(3);
t437 = t479 * t466;
t430 = -qJ(3) * t421 + t437;
t474 = t430 * qJD(2) - t466 * qJD(3);
t395 = t398 * t419;
t445 = t395 * t420 + t422 * t418;
t462 = qJD(1) * t419;
t384 = t445 * qJD(5) + t391 * t418 + t420 * t462;
t444 = t395 * t418 - t422 * t420;
t473 = t444 * qJD(5) - t391 * t420 + t418 * t462;
t470 = pkin(7) - pkin(8);
t463 = t419 * t418;
t461 = qJD(1) * t422;
t460 = qJD(2) * t421;
t455 = qJD(1) * t466;
t454 = qJD(2) * t466;
t443 = -t397 * t418 - t419 * t420;
t436 = -t466 * qJ(3) - t421 * t479;
t433 = -pkin(1) + t436;
t383 = t443 * qJD(5) - t388 * t420 - t418 * t461;
t382 = -t388 * t418 - qJD(5) * t463 + (qJD(5) * t397 + t461) * t420;
t1 = [-t444 * qJD(6) - t391 * pkin(4) + t468 * t390 - t469 * t473 - t464 * t384 + t474 * t419 + (-t470 * t419 + t433 * t422) * qJD(1), -t422 * qJ(3) * t454 + t437 * t462 + (-qJ(3) * t462 + (-qJD(2) * t479 + qJD(3)) * t422) * t421 - t485, -t419 * t455 + t422 * t460, t485 -(-t397 * t420 + t463) * qJD(6) + t464 * t383 + t469 * t382, t382; -t443 * qJD(6) - t388 * pkin(4) + t468 * t389 - t469 * t383 + t464 * t382 - t474 * t422 + (t433 * t419 + t470 * t422) * qJD(1), -t430 * t461 + (t436 * qJD(2) + qJD(3) * t421) * t419 - t486, t419 * t460 + t422 * t455, t486, t445 * qJD(6) + t469 * t384 - t464 * t473, t384; 0, -t474 - t484, t454, t484, t472 * t392 + (t435 * qJD(5) + t420 * qJD(6)) * t399, t399 * qJD(5) * t420 + t392 * t418;];
JaD_transl  = t1;
