% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:05
% EndTime: 2019-02-26 20:58:05
% DurationCPUTime: 0.35s
% Computational Cost: add. (514->70), mult. (616->100), div. (0->0), fcn. (512->9), ass. (0->51)
t281 = cos(qJ(4));
t317 = pkin(4) * t281;
t270 = pkin(3) + t317;
t275 = pkin(9) + qJ(3);
t271 = sin(t275);
t273 = cos(t275);
t279 = sin(qJ(4));
t276 = qJ(4) + pkin(10);
t272 = sin(t276);
t316 = r_i_i_C(2) + qJ(5) + pkin(8);
t291 = t316 * qJD(3) + t272 * qJD(6);
t314 = pkin(4) * qJD(4);
t318 = pkin(4) * t279;
t327 = (-t279 * t314 + t291) * t273 + (pkin(7) + qJ(2) + t318) * qJD(1) - (qJD(3) * t270 - qJD(5)) * t271;
t274 = cos(t276);
t315 = r_i_i_C(3) + qJ(6);
t319 = r_i_i_C(1) + pkin(5);
t320 = t319 * t272 - t315 * t274 + t318;
t325 = t320 * qJD(4) - t291;
t290 = -t315 * t272 - t319 * t274;
t287 = -t270 + t290;
t324 = t287 * t271 + t316 * t273;
t282 = cos(qJ(1));
t313 = t274 * t282;
t280 = sin(qJ(1));
t312 = t280 * t272;
t311 = qJD(1) * t280;
t310 = qJD(1) * t282;
t309 = qJD(3) * t273;
t308 = qJD(3) * t280;
t307 = qJD(3) * t282;
t306 = qJD(4) * t280;
t305 = qJD(4) * t282;
t303 = t274 * qJD(6);
t301 = t271 * t308;
t300 = t271 * t307;
t299 = t272 * t306;
t298 = t274 * t305;
t296 = qJD(1) * t273 - qJD(4);
t294 = (-qJD(4) * t273 + qJD(1)) * t281;
t293 = t273 * t313 + t312;
t289 = t272 * t305 + t274 * t311;
t288 = t272 * t310 + t274 * t306;
t285 = t287 * qJD(3) + qJD(5);
t284 = t281 * t314 - t303 + qJD(2) + (-t270 * t273 - t316 * t271 - cos(pkin(9)) * pkin(2) - pkin(1)) * qJD(1);
t283 = t325 * t271 + t285 * t273;
t258 = t293 * qJD(1) - t273 * t299 - t274 * t301 - t298;
t257 = -t272 * t301 + t288 * t273 - t289;
t256 = t289 * t273 + t274 * t300 - t288;
t255 = t272 * t300 - t273 * t298 - t299 + (t273 * t312 + t313) * qJD(1);
t1 = [-t315 * t257 - t319 * t258 - t327 * t280 + t284 * t282, t310, t283 * t282 - t324 * t311, t293 * qJD(6) - t315 * t256 + t319 * t255 + (t282 * t294 + (t296 * t280 + t300) * t279) * pkin(4), -t271 * t311 + t273 * t307, -t255; -t315 * t255 - t319 * t256 + t284 * t280 + t327 * t282, t311, t283 * t280 + t324 * t310 -(-t280 * t273 * t274 + t272 * t282) * qJD(6) + t315 * t258 - t319 * t257 + (t280 * t294 + (-t296 * t282 + t301) * t279) * pkin(4), t271 * t310 + t273 * t308, t257; 0, 0, t285 * t271 - t325 * t273, -t320 * t309 + (t303 + (t290 - t317) * qJD(4)) * t271, qJD(3) * t271, qJD(4) * t271 * t274 + t272 * t309;];
JaD_transl  = t1;
