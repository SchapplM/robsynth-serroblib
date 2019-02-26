% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:44
% EndTime: 2019-02-26 20:11:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (391->49), mult. (642->86), div. (0->0), fcn. (615->10), ass. (0->45)
t321 = pkin(4) - r_i_i_C(2);
t319 = r_i_i_C(3) + qJ(5);
t288 = sin(pkin(11));
t290 = cos(pkin(11));
t291 = cos(pkin(6));
t295 = cos(qJ(2));
t311 = t291 * t295;
t307 = t290 * t311;
t293 = sin(qJ(2));
t310 = qJD(2) * t293;
t272 = -qJD(2) * t307 + t288 * t310;
t286 = qJD(3) + qJD(4);
t289 = sin(pkin(6));
t315 = t289 * t290;
t323 = -t286 * t315 - t272;
t287 = qJ(3) + qJ(4);
t284 = sin(t287);
t285 = cos(t287);
t294 = cos(qJ(3));
t322 = pkin(3) * t294 + t284 * t319 + t321 * t285 + pkin(2);
t320 = r_i_i_C(1) + pkin(9) + pkin(8);
t318 = t284 * t286;
t317 = t285 * t286;
t316 = t288 * t289;
t292 = sin(qJ(3));
t314 = t289 * t292;
t313 = t289 * t293;
t312 = t291 * t293;
t308 = t285 * t313;
t306 = qJD(2) * t289 * t295;
t304 = t288 * t311 + t290 * t293;
t274 = t304 * qJD(2);
t305 = t286 * t316 - t274;
t277 = t288 * t295 + t290 * t312;
t303 = t288 * t312 - t290 * t295;
t302 = t286 * t291 + t306;
t258 = t277 * t317 + t284 * t323;
t301 = -(-t277 * t285 + t284 * t315) * qJD(5) + t319 * (-t277 * t318 + t285 * t323) - t321 * t258;
t260 = t284 * t305 - t303 * t317;
t300 = -(-t284 * t316 + t285 * t303) * qJD(5) + t319 * (t285 * t305 + t303 * t318) - t321 * t260;
t265 = t284 * t302 + t286 * t308;
t299 = -(-t284 * t291 - t308) * qJD(5) + t319 * (t285 * t302 - t313 * t318) - t321 * t265;
t298 = qJD(2) * t322;
t297 = -pkin(3) * qJD(3) * t292 + qJD(5) * t284 + (-t284 * t321 + t285 * t319) * t286;
t1 = [0, -t274 * t320 - t297 * t304 + t298 * t303 (t274 * t292 + (-t288 * t314 + t294 * t303) * qJD(3)) * pkin(3) + t300, t300, t260, 0; 0, -t320 * t272 - t277 * t298 + t297 * (-t288 * t293 + t307) (t272 * t292 + (-t277 * t294 + t290 * t314) * qJD(3)) * pkin(3) + t301, t301, t258, 0; 0 (-t322 * t310 + (qJD(2) * t320 + t297) * t295) * t289 (-t292 * t306 + (-t291 * t292 - t294 * t313) * qJD(3)) * pkin(3) + t299, t299, t265, 0;];
JaD_transl  = t1;
