% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:53
% EndTime: 2019-02-26 21:34:53
% DurationCPUTime: 0.36s
% Computational Cost: add. (521->73), mult. (636->102), div. (0->0), fcn. (525->10), ass. (0->53)
t283 = cos(qJ(4));
t321 = pkin(4) * t283;
t270 = pkin(3) + t321;
t277 = qJ(2) + pkin(9);
t273 = sin(t277);
t275 = cos(t277);
t320 = r_i_i_C(2) + qJ(5) + pkin(8);
t296 = t320 * t275 - sin(qJ(2)) * pkin(2);
t299 = qJD(4) * t275 - qJD(1);
t307 = t273 * qJD(5);
t276 = qJ(4) + pkin(10);
t272 = sin(t276);
t308 = t272 * qJD(6);
t280 = sin(qJ(4));
t322 = pkin(4) * t280;
t332 = (-t270 * t273 + t296) * qJD(2) - t299 * t322 - qJD(1) * (-qJ(3) - pkin(7)) + t275 * t308 + t307;
t274 = cos(t276);
t319 = r_i_i_C(3) + qJ(6);
t325 = pkin(5) + r_i_i_C(1);
t293 = -t319 * t272 - t325 * t274;
t291 = -t270 + t293;
t288 = t291 * t273 + t296;
t326 = t325 * t272 - t319 * t274 + t322;
t330 = t326 * qJD(4) - t308;
t328 = -t320 * t273 - cos(qJ(2)) * pkin(2);
t285 = cos(qJ(1));
t317 = t274 * t285;
t282 = sin(qJ(1));
t316 = t282 * t272;
t315 = qJD(1) * t282;
t314 = qJD(1) * t285;
t313 = qJD(2) * t275;
t312 = qJD(2) * t282;
t311 = qJD(2) * t285;
t310 = qJD(4) * t282;
t309 = qJD(4) * t285;
t306 = t274 * qJD(6);
t305 = t273 * t312;
t304 = t273 * t311;
t303 = t272 * t310;
t302 = t272 * t309;
t301 = t274 * t309;
t298 = qJD(1) * t275 - qJD(4);
t297 = t299 * t283;
t294 = t275 * t317 + t316;
t292 = t272 * t314 + t274 * t310;
t287 = qJD(4) * t321 - t306 + qJD(3) + (-t270 * t275 - pkin(1) + t328) * qJD(1);
t286 = qJD(5) * t275 + t330 * t273 + (t291 * t275 + t328) * qJD(2);
t259 = t294 * qJD(1) - t274 * t305 - t275 * t303 - t301;
t258 = -t272 * t305 - t274 * t315 + t292 * t275 - t302;
t257 = t275 * t302 + (t275 * t315 + t304) * t274 - t292;
t256 = t272 * t304 - t275 * t301 - t303 + (t275 * t316 + t317) * qJD(1);
t1 = [-t319 * t258 - t325 * t259 - t332 * t282 + t287 * t285, t286 * t285 - t288 * t315, t314, t294 * qJD(6) - t319 * t257 + t325 * t256 + (-t285 * t297 + (t298 * t282 + t304) * t280) * pkin(4), -t273 * t315 + t275 * t311, -t256; -t319 * t256 - t325 * t257 + t287 * t282 + t332 * t285, t286 * t282 + t288 * t314, t315 -(-t282 * t275 * t274 + t272 * t285) * qJD(6) + t319 * t259 - t325 * t258 + (-t282 * t297 + (-t298 * t285 + t305) * t280) * pkin(4), t273 * t314 + t275 * t312, t258; 0, t288 * qJD(2) - t330 * t275 + t307, 0, -t326 * t313 + (t306 + (t293 - t321) * qJD(4)) * t273, qJD(2) * t273, qJD(4) * t273 * t274 + t272 * t313;];
JaD_transl  = t1;
