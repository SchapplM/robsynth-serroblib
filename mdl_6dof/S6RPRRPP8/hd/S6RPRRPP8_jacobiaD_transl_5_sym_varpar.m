% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:21
% EndTime: 2019-02-26 21:00:22
% DurationCPUTime: 0.25s
% Computational Cost: add. (168->55), mult. (518->86), div. (0->0), fcn. (433->6), ass. (0->38)
t241 = cos(qJ(3));
t237 = sin(qJ(4));
t240 = cos(qJ(4));
t266 = r_i_i_C(3) + qJ(5);
t267 = -r_i_i_C(2) + pkin(4);
t246 = t237 * t266 + t240 * t267 + pkin(3);
t238 = sin(qJ(3));
t268 = pkin(8) + r_i_i_C(1);
t254 = t268 * t238;
t272 = t241 * t246 + t254;
t253 = t268 * t241;
t271 = -pkin(3) * t238 - qJ(2) + t253;
t256 = qJD(5) * t237;
t244 = -t256 + (t237 * t267 - t240 * t266) * qJD(4);
t269 = -pkin(1) - pkin(7);
t239 = sin(qJ(1));
t265 = t239 * t237;
t264 = t239 * t240;
t242 = cos(qJ(1));
t263 = t242 * t237;
t262 = qJD(1) * t239;
t261 = qJD(1) * t242;
t260 = qJD(3) * t238;
t259 = qJD(3) * t242;
t258 = qJD(4) * t238;
t257 = qJD(4) * t241;
t255 = t240 * qJD(5);
t252 = t242 * t238 * t240;
t251 = t241 * t259;
t249 = qJD(1) * t238 + qJD(4);
t248 = t238 * t264 + t263;
t245 = t239 * qJD(3) * t241 + t242 * t249;
t243 = t238 * t256 + qJD(2) + (pkin(3) * t241 + t254) * qJD(3);
t230 = -t237 * t262 + t240 * t245 - t258 * t265;
t229 = (qJD(1) + t258) * t264 + t245 * t237;
t228 = -t240 * t251 + (t238 * t263 + t264) * qJD(4) + t248 * qJD(1);
t227 = -qJD(4) * t252 - t237 * t251 - t240 * t261 + t249 * t265;
t1 = [t239 * t255 - t267 * t228 - t266 * t227 + t243 * t242 + (t271 * t239 + t269 * t242) * qJD(1), t261, t272 * t261 + (-t244 * t241 + (-t238 * t246 + t253) * qJD(3)) * t239, qJD(5) * t248 - t229 * t267 + t230 * t266, t229, 0; -t242 * t255 + t267 * t230 + t266 * t229 + t243 * t239 + (t269 * t239 - t271 * t242) * qJD(1), t262 (t246 * t259 + t262 * t268) * t238 + (t246 * t262 + (-qJD(3) * t268 + t244) * t242) * t241 -(t252 - t265) * qJD(5) + t266 * t228 - t267 * t227, t227, 0; 0, 0, -t272 * qJD(3) + t244 * t238 (-t257 * t266 + t260 * t267) * t237 + (-t266 * t260 + (-qJD(4) * t267 + qJD(5)) * t241) * t240, -t237 * t260 + t240 * t257, 0;];
JaD_transl  = t1;
