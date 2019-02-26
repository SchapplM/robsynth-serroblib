% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:19
% EndTime: 2019-02-26 20:45:19
% DurationCPUTime: 0.28s
% Computational Cost: add. (354->54), mult. (562->83), div. (0->0), fcn. (466->8), ass. (0->40)
t241 = sin(qJ(3));
t259 = pkin(3) + pkin(8) + r_i_i_C(2);
t243 = cos(qJ(3));
t240 = sin(qJ(5));
t242 = cos(qJ(5));
t270 = r_i_i_C(3) + qJ(6);
t271 = -r_i_i_C(1) - pkin(5);
t248 = t271 * t240 + t270 * t242 - qJ(4);
t276 = t248 * t243;
t279 = qJD(1) * (t259 * t241 + t276);
t262 = qJD(3) * t243;
t273 = t259 * qJD(3) + qJD(6) * t242 - qJD(4);
t277 = -qJ(4) * t262 + t273 * t241;
t239 = qJ(1) + pkin(9);
t237 = sin(t239);
t238 = cos(t239);
t264 = qJD(1) * t241;
t274 = t237 * t262 + t238 * t264;
t272 = pkin(4) + pkin(7);
t269 = t237 * t240;
t268 = t237 * t242;
t267 = t238 * t240;
t266 = t238 * t242;
t265 = t240 * t241;
t263 = qJD(3) * t241;
t261 = qJD(5) * t243;
t260 = t240 * qJD(6);
t255 = t238 * t262;
t254 = qJD(5) * t241 + qJD(1);
t253 = qJD(5) + t264;
t251 = t238 * t265 + t268;
t250 = -qJ(4) * t241 - t259 * t243 - pkin(2);
t247 = t240 * t262 + t254 * t242;
t245 = (t270 * t240 - t271 * t242) * qJD(5) - t273;
t244 = t245 * t243 + t248 * t263;
t232 = t247 * t238 - t253 * t269;
t231 = -t242 * t255 + t251 * qJD(5) + (t241 * t268 + t267) * qJD(1);
t230 = t247 * t237 + t253 * t267;
t229 = -qJD(5) * t266 - t274 * t242 + t254 * t269;
t1 = [t238 * t260 + t271 * t230 - t270 * t229 + t277 * t237 + (-cos(qJ(1)) * pkin(1) - t272 * t237 + t250 * t238) * qJD(1), 0, t237 * t279 + t244 * t238, -t237 * t264 + t255, t251 * qJD(6) + t271 * t231 + t270 * t232, t231; t237 * t260 - t271 * t232 + t270 * t231 - t277 * t238 + (-sin(qJ(1)) * pkin(1) + t272 * t238 + t250 * t237) * qJD(1), 0, t244 * t237 - t238 * t279, t274 -(-t237 * t265 + t266) * qJD(6) + t270 * t230 + t271 * t229, t229; 0, 0, -qJD(3) * t276 + t245 * t241, t263 (-t270 * t261 - t271 * t263) * t242 + (t270 * t263 + (-t271 * qJD(5) - qJD(6)) * t243) * t240, -t240 * t261 - t242 * t263;];
JaD_transl  = t1;
