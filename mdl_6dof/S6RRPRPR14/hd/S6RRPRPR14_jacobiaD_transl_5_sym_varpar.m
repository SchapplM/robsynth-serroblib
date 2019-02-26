% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR14_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:27
% EndTime: 2019-02-26 21:45:27
% DurationCPUTime: 0.22s
% Computational Cost: add. (312->65), mult. (921->101), div. (0->0), fcn. (892->8), ass. (0->48)
t298 = sin(qJ(4));
t301 = cos(qJ(4));
t330 = r_i_i_C(3) + qJ(5);
t331 = r_i_i_C(2) - pkin(4);
t333 = t331 * t298 + t330 * t301 - qJ(3);
t332 = pkin(3) + pkin(8);
t296 = sin(pkin(6));
t300 = sin(qJ(1));
t329 = t296 * t300;
t302 = cos(qJ(2));
t328 = t296 * t302;
t303 = cos(qJ(1));
t327 = t296 * t303;
t299 = sin(qJ(2));
t326 = t300 * t299;
t325 = t300 * t302;
t324 = t303 * t299;
t323 = t303 * t302;
t322 = qJD(1) * t296;
t321 = qJD(2) * t299;
t320 = qJD(2) * t302;
t319 = -r_i_i_C(1) - pkin(9) - pkin(2);
t318 = t298 * t329;
t297 = cos(pkin(6));
t317 = t297 * t326;
t316 = t301 * t327;
t315 = t297 * t323;
t314 = t300 * t322;
t313 = t303 * t322;
t312 = t296 * t321;
t311 = qJD(2) * t297 + qJD(1);
t285 = -t315 + t326;
t310 = t285 * t301 + t298 * t327;
t287 = t297 * t325 + t324;
t309 = t287 * t298 + t301 * t329;
t308 = -t297 * t301 + t298 * t328;
t307 = t297 * t324 + t325;
t281 = t287 * qJD(1) + t307 * qJD(2);
t272 = -t281 * t301 - qJD(4) * t316 + (qJD(4) * t285 + t314) * t298;
t305 = -t301 * qJD(5) + qJD(3) + (t330 * t298 - t331 * t301) * qJD(4);
t304 = t310 * qJD(4) + t281 * t298 + t301 * t314;
t284 = -t308 * qJD(4) - t301 * t312;
t282 = -qJD(1) * t317 - t300 * t321 + t311 * t323;
t280 = t307 * qJD(1) + t287 * qJD(2);
t279 = -qJD(1) * t315 - t303 * t320 + t311 * t326;
t277 = -t279 * t298 - qJD(4) * t318 + (qJD(4) * t287 + t313) * t301;
t276 = t309 * qJD(4) + t279 * t301 + t298 * t313;
t1 = [t310 * qJD(5) - t281 * qJ(3) - t285 * qJD(3) + t331 * t304 - t330 * t272 + t319 * t282 + (-t303 * pkin(1) - t332 * t329) * qJD(1), -t319 * t279 + t333 * t280 + t305 * (-t317 + t323) -t279, t309 * qJD(5) + t331 * t276 + t330 * t277, t276, 0; -(t287 * t301 - t318) * qJD(5) - t279 * qJ(3) + t287 * qJD(3) - t331 * t277 + t330 * t276 + t319 * t280 + (-t300 * pkin(1) + t332 * t327) * qJD(1), t319 * t281 - t282 * t333 + t305 * t307, t281 -(-t285 * t298 + t316) * qJD(5) + t330 * t304 + t331 * t272, t272, 0; 0 (-t333 * t320 + (t319 * qJD(2) + t305) * t299) * t296, t312, -t308 * qJD(5) + t331 * t284 - t330 * (-t298 * t312 + (t297 * t298 + t301 * t328) * qJD(4)) t284, 0;];
JaD_transl  = t1;
