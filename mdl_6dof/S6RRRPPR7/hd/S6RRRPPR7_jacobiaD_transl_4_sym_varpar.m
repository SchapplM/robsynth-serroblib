% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:02
% EndTime: 2019-02-26 22:07:03
% DurationCPUTime: 0.22s
% Computational Cost: add. (164->54), mult. (510->88), div. (0->0), fcn. (427->6), ass. (0->42)
t248 = sin(qJ(3));
t251 = cos(qJ(3));
t281 = r_i_i_C(3) + qJ(4);
t284 = pkin(3) + r_i_i_C(1);
t285 = t284 * t248 - t281 * t251;
t287 = -t285 * qJD(3) + qJD(4) * t248;
t249 = sin(qJ(2));
t252 = cos(qJ(2));
t283 = pkin(8) + r_i_i_C(2);
t268 = t283 * t252;
t286 = -pkin(2) * t249 + t268;
t259 = -t281 * t248 - t284 * t251;
t256 = -pkin(2) + t259;
t250 = sin(qJ(1));
t280 = t250 * t248;
t279 = t250 * t252;
t253 = cos(qJ(1));
t278 = t253 * t248;
t277 = t253 * t251;
t276 = qJD(1) * t250;
t275 = qJD(1) * t253;
t274 = qJD(2) * t250;
t273 = qJD(2) * t252;
t272 = qJD(2) * t253;
t271 = qJD(3) * t251;
t270 = qJD(3) * t253;
t267 = t249 * t274;
t266 = qJD(3) * t280;
t265 = t249 * t272;
t264 = t248 * t270;
t263 = t251 * t270;
t262 = t252 * t277 + t280;
t261 = t248 * t279 + t277;
t260 = -pkin(2) * t252 - t283 * t249 - pkin(1);
t257 = t248 * t275 + t250 * t271;
t255 = qJD(2) * t256;
t254 = -t283 * qJD(2) - t287;
t237 = t262 * qJD(1) - t251 * t267 - t252 * t266 - t263;
t236 = -t248 * t267 - t251 * t276 + t257 * t252 - t264;
t235 = t252 * t264 + (t252 * t276 + t265) * t251 - t257;
t234 = t261 * qJD(1) + t248 * t265 - t252 * t263 - t266;
t1 = [-t261 * qJD(4) - t284 * t237 - t281 * t236 - t286 * t274 + (-t250 * pkin(7) + t260 * t253) * qJD(1) (t253 * t255 - t283 * t276) * t252 + (t254 * t253 - t256 * t276) * t249, t262 * qJD(4) + t284 * t234 - t281 * t235, -t234, 0, 0; -(t250 * t251 - t252 * t278) * qJD(4) - t284 * t235 - t281 * t234 + t286 * t272 + (t253 * pkin(7) + t260 * t250) * qJD(1) (t250 * t255 + t283 * t275) * t252 + (t254 * t250 + t256 * t275) * t249 -(-t251 * t279 + t278) * qJD(4) + t281 * t237 - t284 * t236, t236, 0, 0; 0, t287 * t252 + (t256 * t249 + t268) * qJD(2), -t285 * t273 + (t259 * qJD(3) + t251 * qJD(4)) * t249, t248 * t273 + t249 * t271, 0, 0;];
JaD_transl  = t1;
