% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:29
% EndTime: 2019-02-26 20:36:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (387->65), mult. (722->97), div. (0->0), fcn. (696->10), ass. (0->50)
t280 = cos(qJ(4));
t278 = sin(qJ(4));
t279 = cos(qJ(5));
t272 = pkin(5) * t279 + pkin(4);
t276 = qJ(5) + qJ(6);
t273 = sin(t276);
t274 = cos(t276);
t323 = r_i_i_C(1) * t274 - r_i_i_C(2) * t273;
t290 = t272 + t323;
t310 = r_i_i_C(3) + pkin(9) + pkin(8);
t297 = t310 * t280;
t283 = -t290 * t278 + t297;
t275 = qJD(5) + qJD(6);
t277 = sin(qJ(5));
t298 = pkin(5) * qJD(5) * t277;
t324 = r_i_i_C(1) * t273 + r_i_i_C(2) * t274;
t318 = t324 * t275 + t298;
t326 = t283 * qJD(4) - t318 * t280;
t296 = t310 * t278;
t325 = t290 * t280 + t296;
t308 = sin(pkin(10));
t309 = cos(pkin(10));
t315 = sin(qJ(1));
t316 = cos(qJ(1));
t264 = t316 * t308 - t315 * t309;
t320 = qJD(1) * t315;
t319 = qJD(1) * t316;
t317 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t263 = -t315 * t308 - t316 * t309;
t307 = t263 * t275;
t305 = t275 * t280;
t261 = t263 * qJD(1);
t302 = qJD(4) * t278;
t289 = t261 * t280 - t264 * t302;
t286 = t289 + t307;
t262 = t264 * qJD(1);
t291 = t264 * t305 + t262;
t304 = (t286 * t273 + t291 * t274) * r_i_i_C(1) + (-t291 * t273 + t286 * t274) * r_i_i_C(2);
t288 = -t262 * t280 - t263 * t302;
t285 = t264 * t275 - t288;
t292 = t263 * t305 + t261;
t259 = -t285 * t273 + t292 * t274;
t260 = t292 * t273 + t285 * t274;
t303 = t259 * r_i_i_C(1) - t260 * r_i_i_C(2);
t301 = qJD(4) * t280;
t300 = qJD(5) * t279;
t299 = qJD(5) * t280;
t287 = t323 * t275 * t278 + t324 * t301;
t282 = t325 * qJD(4) - t278 * t318;
t1 = [(-t262 * t273 + t274 * t307) * r_i_i_C(1) + (-t262 * t274 - t273 * t307) * r_i_i_C(2) - t262 * pkin(7) + (-t262 * t277 + t263 * t300) * pkin(5) - qJ(2) * t320 - (-pkin(3) - t325) * t261 + t326 * t264 + t317 * t316, t319, 0, t283 * t262 + t282 * t263 ((t263 * t299 + t261) * t279 + (-qJD(5) * t264 + t288) * t277) * pkin(5) + t303, t303; t261 * pkin(7) + t260 * r_i_i_C(1) + t259 * r_i_i_C(2) + (t261 * t277 + t264 * t300) * pkin(5) + qJ(2) * t319 + (t280 * t298 + (t272 * t278 - t297) * qJD(4)) * t263 + (t272 * t280 + pkin(3) + t296) * t262 + t317 * t315, t320, 0, -t283 * t261 + t282 * t264 ((t264 * t299 + t262) * t279 + (qJD(5) * t263 + t289) * t277) * pkin(5) + t304, t304; 0, 0, 0, -t326 (t277 * t301 + t278 * t300) * pkin(5) + t287, t287;];
JaD_transl  = t1;
