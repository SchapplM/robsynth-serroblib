% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:50
% EndTime: 2019-02-26 20:50:50
% DurationCPUTime: 0.27s
% Computational Cost: add. (423->55), mult. (494->85), div. (0->0), fcn. (386->10), ass. (0->48)
t242 = cos(qJ(3));
t239 = sin(qJ(5));
t260 = pkin(5) * t239 + qJ(4);
t238 = qJ(5) + qJ(6);
t234 = sin(t238);
t235 = cos(t238);
t288 = r_i_i_C(1) * t234 + r_i_i_C(2) * t235;
t249 = t260 + t288;
t240 = sin(qJ(3));
t264 = pkin(3) + r_i_i_C(3) + pkin(9) + pkin(8);
t280 = t264 * t240;
t285 = -t249 * t242 + t280;
t290 = qJD(1) * t285;
t241 = cos(qJ(5));
t273 = t241 * pkin(5);
t253 = qJD(5) * t273 + qJD(4);
t289 = (-t260 * t242 + t280) * qJD(3) - t253 * t240;
t287 = r_i_i_C(1) * t235 - r_i_i_C(2) * t234;
t268 = qJD(1) * t240;
t284 = t241 * (qJD(5) + t268);
t236 = qJD(5) + qJD(6);
t259 = t236 * t240 + qJD(1);
t266 = qJD(3) * t242;
t279 = t259 * t234 - t235 * t266;
t278 = t234 * t266 + t259 * t235;
t272 = pkin(7) + pkin(4) + t273;
t237 = qJ(1) + pkin(10);
t232 = sin(t237);
t233 = cos(t237);
t258 = -t236 - t268;
t251 = t258 * t235;
t224 = t279 * t232 + t233 * t251;
t252 = t258 * t234;
t225 = -t278 * t232 + t233 * t252;
t270 = -t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
t226 = t232 * t251 - t279 * t233;
t227 = t232 * t252 + t278 * t233;
t269 = t226 * r_i_i_C(1) - t227 * r_i_i_C(2);
t267 = qJD(3) * t240;
t265 = qJD(5) * t239;
t263 = pkin(5) * t265;
t255 = t264 * t242;
t250 = t288 * t236 * t242 + t287 * t267;
t248 = t241 * t266 + (-qJD(5) * t240 - qJD(1)) * t239;
t247 = -t260 * t240 - pkin(2) - t255;
t246 = t287 * t236 + t253;
t244 = t246 * t242 + (-t249 * t240 - t255) * qJD(3);
t1 = [-t233 * t263 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t289 * t232 + (-cos(qJ(1)) * pkin(1) - t272 * t232 + t247 * t233) * qJD(1), 0, t232 * t290 + t244 * t233, -t232 * t268 + t233 * t266 (-t232 * t284 + t248 * t233) * pkin(5) + t269, t269; -t232 * t263 + t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t289 * t233 + (-sin(qJ(1)) * pkin(1) + t272 * t233 + t247 * t232) * qJD(1), 0, t244 * t232 - t233 * t290, t232 * t266 + t233 * t268 (t248 * t232 + t233 * t284) * pkin(5) + t270, t270; 0, 0, -qJD(3) * t285 + t246 * t240, t267 (t241 * t267 + t242 * t265) * pkin(5) + t250, t250;];
JaD_transl  = t1;
