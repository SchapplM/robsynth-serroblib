% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:48
% EndTime: 2019-02-26 19:49:48
% DurationCPUTime: 0.17s
% Computational Cost: add. (136->40), mult. (441->73), div. (0->0), fcn. (420->8), ass. (0->35)
t225 = sin(qJ(4));
t227 = cos(qJ(4));
t246 = r_i_i_C(3) + qJ(5);
t247 = pkin(4) - r_i_i_C(2);
t248 = t247 * t225 - t246 * t227 + qJ(3);
t221 = sin(pkin(10));
t223 = cos(pkin(10));
t226 = sin(qJ(2));
t224 = cos(pkin(6));
t228 = cos(qJ(2));
t240 = t224 * t228;
t215 = t221 * t226 - t223 * t240;
t245 = t215 * t225;
t222 = sin(pkin(6));
t244 = t222 * t225;
t243 = t222 * t227;
t242 = t222 * t228;
t241 = t224 * t226;
t239 = qJD(2) * t226;
t238 = qJD(2) * t228;
t237 = -pkin(2) - pkin(8) - r_i_i_C(1);
t236 = t221 * t239;
t235 = t222 * t239;
t234 = t223 * t238;
t217 = t221 * t240 + t223 * t226;
t233 = t217 * t225 + t221 * t243;
t232 = t221 * t228 + t223 * t241;
t231 = -t224 * t227 + t225 * t242;
t229 = -qJD(5) * t227 + qJD(3) + (t246 * t225 + t247 * t227) * qJD(4);
t214 = -t224 * t236 + t234;
t212 = t232 * qJD(2);
t209 = -t231 * qJD(4) - t227 * t235;
t207 = -qJD(4) * t245 + (qJD(4) * t222 * t223 + t212) * t227;
t205 = t233 * qJD(4) - t214 * t227;
t1 = [0, t237 * t214 - t248 * t217 * qJD(2) + t229 * (-t221 * t241 + t223 * t228) t214, t233 * qJD(5) - t247 * t205 - t246 * (-t214 * t225 + (-t217 * t227 + t221 * t244) * qJD(4)) t205, 0; 0, t237 * t212 - t248 * (-t224 * t234 + t236) + t229 * t232, t212 -(t223 * t243 - t245) * qJD(5) + t247 * t207 + t246 * (t212 * t225 + (t215 * t227 + t223 * t244) * qJD(4)) -t207, 0; 0 (t248 * t238 + (t237 * qJD(2) + t229) * t226) * t222, t235, -t231 * qJD(5) - t247 * t209 - t246 * (-t225 * t235 + (t224 * t225 + t227 * t242) * qJD(4)) t209, 0;];
JaD_transl  = t1;
