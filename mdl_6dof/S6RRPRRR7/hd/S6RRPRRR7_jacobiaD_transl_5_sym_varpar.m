% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:37
% EndTime: 2019-02-26 21:57:37
% DurationCPUTime: 0.42s
% Computational Cost: add. (341->73), mult. (1050->123), div. (0->0), fcn. (979->8), ass. (0->52)
t342 = sin(qJ(4));
t346 = cos(qJ(2));
t343 = sin(qJ(2));
t385 = cos(qJ(4));
t372 = t343 * t385;
t394 = -t346 * t342 + t372;
t388 = pkin(2) + pkin(3);
t392 = -qJ(3) * t346 + t388 * t343;
t400 = t392 * qJD(2) - t343 * qJD(3);
t366 = qJD(4) * t385;
t367 = qJD(2) * t385;
t375 = qJD(4) * t342;
t399 = t346 * t375 + (-t366 + t367) * t343;
t322 = t343 * t342 + t346 * t385;
t378 = qJD(2) * t343;
t316 = t322 * qJD(4) - t342 * t378 - t346 * t367;
t344 = sin(qJ(1));
t347 = cos(qJ(1));
t395 = t394 * qJD(1);
t314 = t316 * t344 - t395 * t347;
t377 = qJD(2) * t346;
t317 = t342 * t377 - t399;
t352 = qJD(1) * t322;
t315 = t317 * t344 + t347 * t352;
t341 = sin(qJ(5));
t345 = cos(qJ(5));
t358 = t345 * r_i_i_C(1) - t341 * r_i_i_C(2) + pkin(4);
t386 = -r_i_i_C(3) - pkin(9);
t384 = t345 * r_i_i_C(2);
t393 = qJD(5) * (t341 * r_i_i_C(1) + t384);
t398 = -t394 * t344 * t393 - t358 * t314 - t386 * t315;
t376 = qJD(2) * t347;
t368 = t346 * t376;
t312 = -t342 * t368 + t344 * t352 + t399 * t347;
t321 = t322 * t347;
t381 = t346 * t347;
t313 = -t347 * t343 * t375 + qJD(2) * t321 - t395 * t344 - t366 * t381;
t397 = (t342 * t381 - t347 * t372) * t393 + t386 * t312 + t358 * t313;
t396 = t386 * t316 - t358 * t317 + t322 * t393;
t387 = pkin(7) - pkin(8);
t380 = qJD(1) * t344;
t379 = qJD(1) * t347;
t374 = qJD(5) * t394;
t319 = t322 * t344;
t360 = t319 * t345 + t341 * t347;
t359 = t319 * t341 - t345 * t347;
t357 = -qJ(3) * t343 - t388 * t346;
t353 = -pkin(1) + t357;
t340 = t341 * t380;
t311 = -t341 * t379 - t312 * t345 + (-t321 * t341 - t344 * t345) * qJD(5);
t310 = -t345 * t379 + t312 * t341 + (-t321 * t345 + t341 * t344) * qJD(5);
t1 = [t340 * r_i_i_C(1) - t358 * t315 + t386 * t314 + (t359 * r_i_i_C(1) + t360 * r_i_i_C(2)) * qJD(5) + t400 * t344 + ((t384 - t387) * t344 + t353 * t347) * qJD(1) (-qJ(3) * t376 + t388 * t380) * t343 + (-qJ(3) * t380 + (-t388 * qJD(2) + qJD(3)) * t347) * t346 - t397, -t343 * t380 + t368, t397, t310 * r_i_i_C(1) - t311 * r_i_i_C(2), 0; -t312 * pkin(4) + t311 * r_i_i_C(1) + t310 * r_i_i_C(2) + t386 * t313 - t400 * t347 + (t353 * t344 + t387 * t347) * qJD(1), -t392 * t379 + (t357 * qJD(2) + qJD(3) * t346) * t344 - t398, t343 * t379 + t344 * t377, t398 (-t315 * t341 - t345 * t380) * r_i_i_C(1) + (-t315 * t345 + t340) * r_i_i_C(2) + (-t360 * r_i_i_C(1) + t359 * r_i_i_C(2)) * qJD(5), 0; 0, -t400 - t396, t378, t396 (t316 * t345 + t341 * t374) * r_i_i_C(2) + (t316 * t341 - t345 * t374) * r_i_i_C(1), 0;];
JaD_transl  = t1;
