% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:42
% EndTime: 2019-02-26 20:30:43
% DurationCPUTime: 0.26s
% Computational Cost: add. (459->61), mult. (520->91), div. (0->0), fcn. (435->9), ass. (0->43)
t269 = pkin(10) + qJ(4);
t265 = sin(t269);
t267 = cos(t269);
t301 = pkin(8) + r_i_i_C(2);
t288 = t301 * t267;
t306 = (-pkin(4) * t265 + t288) * qJD(4);
t272 = sin(qJ(5));
t273 = cos(qJ(5));
t299 = r_i_i_C(3) + qJ(6);
t302 = pkin(5) + r_i_i_C(1);
t303 = t302 * t272 - t299 * t273;
t305 = -t303 * qJD(5) + qJD(6) * t272;
t280 = -t299 * t272 - t302 * t273;
t276 = -pkin(4) + t280;
t270 = qJ(1) + pkin(9);
t266 = sin(t270);
t298 = t266 * t272;
t297 = t266 * t273;
t268 = cos(t270);
t296 = t268 * t272;
t295 = t268 * t273;
t294 = qJD(1) * t266;
t293 = qJD(1) * t268;
t292 = qJD(4) * t272;
t291 = qJD(5) * t272;
t290 = qJD(5) * t273;
t287 = t265 * t292;
t286 = qJD(4) * t265 * t273;
t285 = t266 * t291;
t284 = t268 * t290;
t283 = t267 * t295 + t298;
t282 = t267 * t298 + t295;
t281 = -pkin(4) * t267 - t301 * t265 - cos(pkin(10)) * pkin(3) - pkin(2);
t278 = t266 * t290 + t272 * t293;
t277 = t268 * t291 + t273 * t294;
t275 = qJD(4) * t276;
t274 = -t301 * qJD(4) - t305;
t271 = -pkin(7) - qJ(3);
t253 = t283 * qJD(1) - t266 * t286 - t267 * t285 - t284;
t252 = -t266 * t287 + t278 * t267 - t277;
t251 = t277 * t267 + t268 * t286 - t278;
t250 = t282 * qJD(1) - t267 * t284 + t268 * t287 - t285;
t1 = [-t282 * qJD(6) + t268 * qJD(3) - t302 * t253 - t299 * t252 - t266 * t306 + (-cos(qJ(1)) * pkin(1) + t266 * t271 + t281 * t268) * qJD(1), 0, t293 (t268 * t275 - t301 * t294) * t267 + (t274 * t268 - t276 * t294) * t265, t283 * qJD(6) + t302 * t250 - t299 * t251, -t250; -(-t267 * t296 + t297) * qJD(6) + t266 * qJD(3) - t302 * t251 - t299 * t250 + t268 * t306 + (-sin(qJ(1)) * pkin(1) - t268 * t271 + t281 * t266) * qJD(1), 0, t294 (t266 * t275 + t301 * t293) * t267 + (t274 * t266 + t276 * t293) * t265 -(-t267 * t297 + t296) * qJD(6) + t299 * t253 - t302 * t252, t252; 0, 0, 0, t305 * t267 + (t276 * t265 + t288) * qJD(4), -t303 * t267 * qJD(4) + (t280 * qJD(5) + qJD(6) * t273) * t265, t265 * t290 + t267 * t292;];
JaD_transl  = t1;
