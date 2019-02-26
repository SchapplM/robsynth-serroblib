% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR11_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR11_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:52
% EndTime: 2019-02-26 22:21:52
% DurationCPUTime: 0.25s
% Computational Cost: add. (278->57), mult. (828->96), div. (0->0), fcn. (802->8), ass. (0->43)
t297 = cos(pkin(6));
t299 = sin(qJ(2));
t300 = sin(qJ(1));
t324 = t300 * t299;
t318 = t297 * t324;
t302 = cos(qJ(2));
t303 = cos(qJ(1));
t321 = t303 * t302;
t283 = -qJD(1) * t318 - qJD(2) * t324 + (qJD(2) * t297 + qJD(1)) * t321;
t296 = sin(pkin(6));
t325 = t296 * t303;
t334 = -qJD(3) * t325 + t283;
t322 = t303 * t299;
t323 = t300 * t302;
t287 = t297 * t322 + t323;
t320 = qJD(1) * t296;
t333 = -qJD(3) * t287 + t300 * t320;
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t332 = t298 * t333 + t334 * t301;
t328 = r_i_i_C(3) + qJ(4);
t330 = pkin(3) + r_i_i_C(1);
t331 = t328 * t298 + t301 * t330 + pkin(2);
t329 = pkin(9) + r_i_i_C(2);
t327 = t296 * t300;
t326 = t296 * t301;
t316 = t303 * t320;
t315 = qJD(2) * t296 * t302;
t307 = t318 - t321;
t312 = t298 * t307 + t300 * t326;
t311 = t298 * t327 - t301 * t307;
t310 = t297 * t298 + t299 * t326;
t309 = t297 * t321 - t324;
t308 = t297 * t323 + t322;
t276 = t334 * t298 - t333 * t301;
t304 = qJD(4) * t298 + (-t298 * t330 + t328 * t301) * qJD(3);
t284 = qJD(3) * t310 + t298 * t315;
t282 = qJD(1) * t308 + qJD(2) * t287;
t281 = qJD(1) * t287 + qJD(2) * t308;
t280 = -qJD(1) * t309 + qJD(2) * t307;
t275 = qJD(3) * t312 - t281 * t301 + t298 * t316;
t274 = qJD(3) * t311 - t281 * t298 - t301 * t316;
t1 = [-(t287 * t298 + t301 * t325) * qJD(4) - t283 * pkin(2) - t329 * t282 - t330 * t332 - t328 * t276 + (-t303 * pkin(1) - pkin(8) * t327) * qJD(1), t331 * t280 - t329 * t281 - t304 * t308, t311 * qJD(4) - t274 * t330 + t328 * t275, t274, 0, 0; -t312 * qJD(4) - t281 * pkin(2) - t329 * t280 + t330 * t275 + t328 * t274 + (-t300 * pkin(1) + pkin(8) * t325) * qJD(1), -t282 * t331 + t283 * t329 + t304 * t309 -(-t287 * t301 + t298 * t325) * qJD(4) + t328 * t332 - t330 * t276, t276, 0, 0; 0 (t304 * t302 + (-t299 * t331 + t302 * t329) * qJD(2)) * t296, t310 * qJD(4) + t328 * (t301 * t315 + (-t296 * t298 * t299 + t297 * t301) * qJD(3)) - t330 * t284, t284, 0, 0;];
JaD_transl  = t1;
