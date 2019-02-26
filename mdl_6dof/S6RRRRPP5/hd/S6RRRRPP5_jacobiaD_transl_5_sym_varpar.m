% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:44
% EndTime: 2019-02-26 22:27:44
% DurationCPUTime: 0.34s
% Computational Cost: add. (523->70), mult. (738->104), div. (0->0), fcn. (620->8), ass. (0->55)
t299 = qJ(3) + qJ(4);
t297 = cos(t299);
t346 = r_i_i_C(3) + qJ(5);
t354 = t346 * t297;
t298 = qJD(3) + qJD(4);
t296 = sin(t299);
t300 = sin(qJ(3));
t345 = pkin(3) * qJD(3);
t347 = r_i_i_C(2) + pkin(9) + pkin(8);
t309 = t347 * qJD(2) + qJD(5) * t296 - t300 * t345;
t349 = pkin(4) + r_i_i_C(1);
t332 = t349 * t296;
t307 = (t332 - t354) * t298 - t309;
t303 = cos(qJ(3));
t295 = t303 * pkin(3) + pkin(2);
t304 = cos(qJ(2));
t301 = sin(qJ(2));
t336 = qJD(2) * t301;
t353 = -t295 * t336 + t309 * t304;
t314 = -t346 * t296 - t349 * t297;
t311 = -t295 + t314;
t302 = sin(qJ(1));
t305 = cos(qJ(1));
t340 = t305 * t297;
t317 = t302 * t296 + t304 * t340;
t342 = t302 * t304;
t350 = t296 * t342 + t340;
t348 = pkin(3) * t300;
t344 = t298 * t301;
t341 = t305 * t296;
t339 = qJD(1) * t302;
t338 = qJD(1) * t304;
t337 = qJD(1) * t305;
t335 = qJD(2) * t304;
t334 = t297 * qJD(5);
t333 = t303 * t345;
t331 = t297 * t344;
t329 = t298 * t341;
t327 = t301 * t334 + t335 * t354;
t326 = pkin(7) + t348;
t323 = t302 * t336;
t322 = t305 * t336;
t321 = -qJD(3) + t338;
t320 = (-qJD(3) * t304 + qJD(1)) * t303;
t318 = t333 - t334;
t316 = -t295 * t304 - t347 * t301 - pkin(1);
t315 = t302 * t298 * t297 + t296 * t337;
t276 = qJD(1) * t350 + t296 * t322 - t298 * t317;
t277 = t304 * t329 + (t302 * t338 + t322) * t297 - t315;
t313 = t317 * qJD(5) + t349 * t276 - t346 * t277;
t278 = -t296 * t323 - t297 * t339 + t315 * t304 - t329;
t279 = t317 * qJD(1) - t297 * t323 - t298 * t350;
t312 = -(-t297 * t342 + t341) * qJD(5) + t346 * t279 - t349 * t278;
t310 = qJD(2) * t311;
t1 = [t318 * t305 - t349 * t279 - t346 * t278 - t353 * t302 + (-t326 * t302 + t316 * t305) * qJD(1) (t305 * t310 - t347 * t339) * t304 + (t307 * t305 - t311 * t339) * t301 (t305 * t320 + (t321 * t302 + t322) * t300) * pkin(3) + t313, t313, -t276, 0; t318 * t302 - t349 * t277 - t346 * t276 + t353 * t305 + (t316 * t302 + t326 * t305) * qJD(1) (t302 * t310 + t347 * t337) * t304 + (t307 * t302 + t311 * t337) * t301 (t302 * t320 + (-t321 * t305 + t323) * t300) * pkin(3) + t312, t312, t278, 0; 0, t301 * t310 - t307 * t304 (-t332 - t348) * t335 + (t314 * t298 - t333) * t301 + t327, -t349 * t331 + (-t349 * t335 - t346 * t344) * t296 + t327, t296 * t335 + t331, 0;];
JaD_transl  = t1;
