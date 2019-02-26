% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:54
% EndTime: 2019-02-26 22:13:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (164->54), mult. (510->88), div. (0->0), fcn. (427->6), ass. (0->42)
t249 = sin(qJ(3));
t252 = cos(qJ(3));
t282 = r_i_i_C(3) + qJ(4);
t285 = pkin(3) + r_i_i_C(1);
t286 = t285 * t249 - t282 * t252;
t288 = -t286 * qJD(3) + qJD(4) * t249;
t250 = sin(qJ(2));
t253 = cos(qJ(2));
t284 = pkin(8) + r_i_i_C(2);
t269 = t284 * t253;
t287 = -pkin(2) * t250 + t269;
t260 = -t282 * t249 - t285 * t252;
t257 = -pkin(2) + t260;
t251 = sin(qJ(1));
t281 = t251 * t249;
t280 = t251 * t253;
t254 = cos(qJ(1));
t279 = t254 * t249;
t278 = t254 * t252;
t277 = qJD(1) * t251;
t276 = qJD(1) * t254;
t275 = qJD(2) * t251;
t274 = qJD(2) * t253;
t273 = qJD(2) * t254;
t272 = qJD(3) * t252;
t271 = qJD(3) * t254;
t268 = t250 * t275;
t267 = qJD(3) * t281;
t266 = t250 * t273;
t265 = t249 * t271;
t264 = t252 * t271;
t263 = t253 * t278 + t281;
t262 = t249 * t280 + t278;
t261 = -pkin(2) * t253 - t284 * t250 - pkin(1);
t258 = t249 * t276 + t251 * t272;
t256 = qJD(2) * t257;
t255 = -t284 * qJD(2) - t288;
t238 = t263 * qJD(1) - t252 * t268 - t253 * t267 - t264;
t237 = -t249 * t268 - t252 * t277 + t258 * t253 - t265;
t236 = t253 * t265 + (t253 * t277 + t266) * t252 - t258;
t235 = t262 * qJD(1) + t249 * t266 - t253 * t264 - t267;
t1 = [-t262 * qJD(4) - t285 * t238 - t282 * t237 - t287 * t275 + (-t251 * pkin(7) + t261 * t254) * qJD(1) (t254 * t256 - t284 * t277) * t253 + (t255 * t254 - t257 * t277) * t250, t263 * qJD(4) + t285 * t235 - t282 * t236, -t235, 0, 0; -(t251 * t252 - t253 * t279) * qJD(4) - t285 * t236 - t282 * t235 + t287 * t273 + (t254 * pkin(7) + t261 * t251) * qJD(1) (t251 * t256 + t284 * t276) * t253 + (t255 * t251 + t257 * t276) * t250 -(-t252 * t280 + t279) * qJD(4) + t282 * t238 - t285 * t237, t237, 0, 0; 0, t288 * t253 + (t257 * t250 + t269) * qJD(2), -t286 * t274 + (t260 * qJD(3) + t252 * qJD(4)) * t250, t249 * t274 + t250 * t272, 0, 0;];
JaD_transl  = t1;
