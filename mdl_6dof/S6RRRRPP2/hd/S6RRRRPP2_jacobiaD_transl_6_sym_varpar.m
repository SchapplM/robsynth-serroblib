% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP2
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
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:32
% EndTime: 2019-02-26 22:25:32
% DurationCPUTime: 0.39s
% Computational Cost: add. (616->78), mult. (887->104), div. (0->0), fcn. (714->8), ass. (0->57)
t334 = r_i_i_C(3) + qJ(6);
t275 = sin(qJ(4));
t326 = r_i_i_C(2) + qJ(5);
t336 = t275 * t326;
t278 = cos(qJ(4));
t280 = cos(qJ(1));
t315 = qJD(4) * t280;
t277 = sin(qJ(1));
t318 = qJD(1) * t277;
t291 = t275 * t315 + t278 * t318;
t314 = qJD(5) * t275;
t316 = qJD(4) * t278;
t335 = -t316 * t326 - t314;
t312 = r_i_i_C(1) + pkin(5) + pkin(4);
t274 = qJ(2) + qJ(3);
t271 = sin(t274);
t272 = cos(t274);
t273 = qJD(2) + qJD(3);
t276 = sin(qJ(2));
t324 = pkin(2) * qJD(2);
t310 = t276 * t324;
t311 = pkin(9) - t334;
t332 = (t273 * t311 + t314) * t272 - (pkin(3) * t273 + qJD(6)) * t271 - qJD(1) * (-pkin(8) - pkin(7)) - t310;
t331 = t271 * t312;
t328 = pkin(2) * t276;
t327 = pkin(9) * t272;
t323 = t272 * t273;
t322 = t273 * t277;
t321 = t273 * t280;
t320 = t277 * t275;
t319 = t280 * t278;
t317 = qJD(1) * t280;
t313 = t278 * qJD(5);
t309 = t271 * t322;
t308 = t271 * t321;
t307 = t271 * t318;
t306 = t272 * t318;
t303 = qJD(4) * t320;
t301 = t278 * t315;
t299 = t312 * t275;
t294 = t303 * t331 + t334 * t309 + t317 * t327;
t292 = t272 * t319 + t320;
t290 = t275 * t317 + t277 * t316;
t289 = -t278 * t312 - t336;
t288 = -pkin(3) + t289;
t279 = cos(qJ(2));
t287 = -t313 + (-t279 * pkin(2) - pkin(3) * t272 - t271 * t311 - pkin(1)) * qJD(1);
t286 = t334 * (t306 + t308) + t291 * t331 + (pkin(3) + t336) * t307;
t285 = t271 * t288 - t272 * t334;
t284 = (-pkin(9) * t273 + t335) * t271 + (t273 * t288 - qJD(6)) * t272;
t283 = -t279 * t324 + t284;
t282 = pkin(9) * t323 - t271 * qJD(6) + t273 * t285 + (-qJD(4) * t299 - t335) * t272;
t236 = qJD(1) * t292 - t272 * t303 - t278 * t309 - t301;
t235 = t272 * t290 - t275 * t309 - t291;
t234 = t272 * t291 + t278 * t308 - t290;
t233 = t275 * t308 - t272 * t301 - t303 + (t272 * t320 + t319) * qJD(1);
t1 = [-t326 * t235 - t312 * t236 - t332 * t277 + t287 * t280, t283 * t280 + t286 + (-t327 + t328) * t318, -pkin(9) * t306 + t280 * t284 + t286, qJD(5) * t292 + t233 * t312 - t234 * t326, -t233, -t272 * t321 + t307; -t326 * t233 - t312 * t234 + t287 * t277 + t332 * t280 (t285 - t328) * t317 + t283 * t277 + t294, t277 * t284 + t285 * t317 + t294 -(-t277 * t272 * t278 + t280 * t275) * qJD(5) + t326 * t236 - t312 * t235, t235, -t271 * t317 - t272 * t322; 0, t282 - t310, t282 (t278 * t326 - t299) * t323 + (qJD(4) * t289 + t313) * t271, t271 * t316 + t275 * t323, -t273 * t271;];
JaD_transl  = t1;
