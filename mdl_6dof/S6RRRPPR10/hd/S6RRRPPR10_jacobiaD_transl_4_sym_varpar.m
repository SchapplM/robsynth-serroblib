% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR10_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:58
% EndTime: 2019-02-26 22:08:58
% DurationCPUTime: 0.23s
% Computational Cost: add. (278->59), mult. (828->98), div. (0->0), fcn. (802->8), ass. (0->45)
t285 = cos(pkin(6));
t287 = sin(qJ(2));
t288 = sin(qJ(1));
t311 = t288 * t287;
t305 = t285 * t311;
t290 = cos(qJ(2));
t291 = cos(qJ(1));
t307 = t291 * t290;
t273 = -qJD(1) * t305 - qJD(2) * t311 + (qJD(2) * t285 + qJD(1)) * t307;
t289 = cos(qJ(3));
t308 = t291 * t287;
t309 = t288 * t290;
t278 = t285 * t308 + t309;
t286 = sin(qJ(3));
t284 = sin(pkin(6));
t313 = t284 * t291;
t302 = t278 * t286 + t289 * t313;
t315 = t284 * t288;
t306 = t286 * t315;
t322 = -qJD(1) * t306 + t302 * qJD(3) - t273 * t289;
t312 = t286 * t291;
t301 = -t278 * t289 + t284 * t312;
t314 = t284 * t289;
t304 = qJD(1) * t314;
t321 = t301 * qJD(3) - t273 * t286 + t288 * t304;
t317 = r_i_i_C(3) + qJ(4);
t319 = pkin(3) - r_i_i_C(2);
t320 = t317 * t286 + t319 * t289 + pkin(2);
t318 = -r_i_i_C(1) - pkin(9);
t296 = t305 - t307;
t316 = t296 * t286;
t310 = t288 * t289;
t303 = qJD(2) * t284 * t290;
t300 = -t289 * t296 + t306;
t299 = t285 * t286 + t287 * t314;
t298 = t285 * t307 - t311;
t297 = t285 * t309 + t308;
t292 = qJD(4) * t286 + (-t319 * t286 + t317 * t289) * qJD(3);
t274 = t299 * qJD(3) + t286 * t303;
t272 = t297 * qJD(1) + t278 * qJD(2);
t271 = t278 * qJD(1) + t297 * qJD(2);
t270 = -t298 * qJD(1) + t296 * qJD(2);
t265 = -t271 * t289 + qJD(3) * t316 + (qJD(1) * t312 + qJD(3) * t310) * t284;
t264 = t300 * qJD(3) - t271 * t286 - t291 * t304;
t1 = [-t302 * qJD(4) - t273 * pkin(2) + t318 * t272 + t319 * t322 + t317 * t321 + (-t291 * pkin(1) - pkin(8) * t315) * qJD(1), t320 * t270 + t318 * t271 - t292 * t297, t300 * qJD(4) - t319 * t264 + t317 * t265, t264, 0, 0; -(t284 * t310 + t316) * qJD(4) - t271 * pkin(2) + t318 * t270 + t319 * t265 + t317 * t264 + (-t288 * pkin(1) + pkin(8) * t313) * qJD(1), -t272 * t320 - t318 * t273 + t292 * t298, -t301 * qJD(4) - t317 * t322 + t319 * t321, -t321, 0, 0; 0 (t292 * t290 + (-t287 * t320 - t318 * t290) * qJD(2)) * t284, t299 * qJD(4) + t317 * (t289 * t303 + (-t284 * t286 * t287 + t285 * t289) * qJD(3)) - t319 * t274, t274, 0, 0;];
JaD_transl  = t1;
