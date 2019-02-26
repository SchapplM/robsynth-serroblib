% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:28
% EndTime: 2019-02-26 20:06:29
% DurationCPUTime: 0.17s
% Computational Cost: add. (152->37), mult. (504->70), div. (0->0), fcn. (488->10), ass. (0->34)
t257 = sin(pkin(11));
t260 = cos(pkin(11));
t263 = sin(qJ(2));
t261 = cos(pkin(6));
t265 = cos(qJ(2));
t277 = t261 * t265;
t286 = -t257 * t263 + t260 * t277;
t278 = t261 * t263;
t251 = t257 * t265 + t260 * t278;
t264 = cos(qJ(3));
t258 = sin(pkin(6));
t262 = sin(qJ(3));
t280 = t258 * t262;
t285 = -t251 * t264 + t260 * t280;
t256 = sin(pkin(12));
t259 = cos(pkin(12));
t273 = t259 * r_i_i_C(1) - t256 * r_i_i_C(2) + pkin(3);
t283 = r_i_i_C(3) + qJ(4);
t284 = t283 * t262 + t273 * t264 + pkin(2);
t279 = t258 * t264;
t274 = qJD(2) * t258 * t265;
t272 = t256 * r_i_i_C(1) + t259 * r_i_i_C(2) + pkin(8);
t269 = t257 * t278 - t260 * t265;
t271 = t257 * t280 - t264 * t269;
t270 = t257 * t277 + t260 * t263;
t268 = t261 * t262 + t263 * t279;
t267 = qJD(2) * t284;
t266 = t262 * qJD(4) + (-t273 * t262 + t283 * t264) * qJD(3);
t248 = t270 * qJD(2);
t246 = t286 * qJD(2);
t244 = t268 * qJD(3) + t262 * t274;
t242 = t271 * qJD(3) - t248 * t262;
t240 = -t285 * qJD(3) + t246 * t262;
t1 = [0, -t272 * t248 - t266 * t270 + t269 * t267, t271 * qJD(4) + t283 * (-t248 * t264 + (t257 * t279 + t262 * t269) * qJD(3)) - t273 * t242, t242, 0, 0; 0, t272 * t246 - t251 * t267 + t266 * t286, -t285 * qJD(4) + t283 * (t246 * t264 + (-t251 * t262 - t260 * t279) * qJD(3)) - t273 * t240, t240, 0, 0; 0 (t266 * t265 + (-t284 * t263 + t272 * t265) * qJD(2)) * t258, t268 * qJD(4) + t283 * (t264 * t274 + (t261 * t264 - t263 * t280) * qJD(3)) - t273 * t244, t244, 0, 0;];
JaD_transl  = t1;
