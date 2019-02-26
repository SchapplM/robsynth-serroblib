% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:38
% EndTime: 2019-02-26 20:07:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
t221 = sin(pkin(11));
t223 = cos(pkin(11));
t226 = sin(qJ(2));
t224 = cos(pkin(6));
t228 = cos(qJ(2));
t238 = t224 * t228;
t247 = -t221 * t226 + t223 * t238;
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t243 = r_i_i_C(3) + qJ(4);
t245 = pkin(3) - r_i_i_C(2);
t246 = t243 * t225 + t227 * t245 + pkin(2);
t244 = pkin(8) + r_i_i_C(1);
t222 = sin(pkin(6));
t241 = t222 * t225;
t240 = t222 * t227;
t239 = t224 * t226;
t236 = qJD(2) * t222 * t228;
t217 = t221 * t228 + t223 * t239;
t235 = -t217 * t227 + t223 * t241;
t232 = t221 * t239 - t223 * t228;
t234 = t221 * t241 - t227 * t232;
t233 = t221 * t238 + t223 * t226;
t231 = t224 * t225 + t226 * t240;
t230 = qJD(2) * t246;
t229 = qJD(4) * t225 + (-t225 * t245 + t243 * t227) * qJD(3);
t214 = t233 * qJD(2);
t212 = t247 * qJD(2);
t210 = qJD(3) * t231 + t225 * t236;
t208 = qJD(3) * t234 - t214 * t225;
t206 = -qJD(3) * t235 + t212 * t225;
t1 = [0, -t214 * t244 - t229 * t233 + t230 * t232, t234 * qJD(4) + t243 * (-t214 * t227 + (t221 * t240 + t225 * t232) * qJD(3)) - t245 * t208, t208, 0, 0; 0, t244 * t212 - t217 * t230 + t229 * t247, -t235 * qJD(4) + t243 * (t212 * t227 + (-t217 * t225 - t223 * t240) * qJD(3)) - t245 * t206, t206, 0, 0; 0 (t229 * t228 + (-t246 * t226 + t244 * t228) * qJD(2)) * t222, t231 * qJD(4) + t243 * (t227 * t236 + (t224 * t227 - t226 * t241) * qJD(3)) - t245 * t210, t210, 0, 0;];
JaD_transl  = t1;
