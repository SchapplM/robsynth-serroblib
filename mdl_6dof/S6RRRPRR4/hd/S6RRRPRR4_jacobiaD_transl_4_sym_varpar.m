% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:43
% EndTime: 2019-02-26 22:17:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (217->43), mult. (285->62), div. (0->0), fcn. (208->8), ass. (0->39)
t224 = qJ(2) + qJ(3);
t222 = cos(t224);
t259 = r_i_i_C(3) + qJ(4);
t244 = t259 * t222;
t221 = sin(t224);
t219 = t221 * qJD(4);
t223 = qJD(2) + qJD(3);
t225 = sin(pkin(11));
t226 = cos(pkin(11));
t260 = r_i_i_C(2) * t225;
t267 = r_i_i_C(1) * t226 + pkin(3);
t237 = t267 - t260;
t227 = sin(qJ(2));
t258 = pkin(2) * qJD(2);
t250 = t227 * t258;
t268 = (-t237 * t221 + t244) * t223 + (r_i_i_C(1) * t225 + r_i_i_C(2) * t226 + pkin(7) + pkin(8)) * qJD(1) + t219 - t250;
t257 = t222 * t223;
t266 = qJD(4) * t222 + t257 * t260;
t262 = pkin(2) * t227;
t228 = sin(qJ(1));
t256 = t222 * t228;
t230 = cos(qJ(1));
t255 = t223 * t230;
t254 = qJD(1) * t228;
t253 = qJD(1) * t230;
t251 = t221 * t260;
t248 = t221 * t253;
t246 = t221 * t254;
t245 = t259 * t221;
t243 = t259 * t228;
t241 = t266 * t230 + t267 * t246;
t240 = t267 * t223;
t239 = t267 * t230;
t238 = t266 * t228 + t253 * t244 + t248 * t260;
t234 = -t221 * t240 + t223 * t251 + t259 * t257 + t219;
t229 = cos(qJ(2));
t233 = qJD(1) * (-pkin(2) * t229 - t237 * t222 - pkin(1) - t245);
t232 = -t229 * t258 + (-t222 * t267 - t245) * t223;
t1 = [-t268 * t228 + t230 * t233 (-t244 - t251 + t262) * t254 + t232 * t230 + t241 (-t254 * t260 - t259 * t255) * t221 + (-qJD(1) * t243 - t223 * t239) * t222 + t241, t222 * t255 - t246, 0, 0; t228 * t233 + t268 * t230 (-t221 * t267 - t262) * t253 + t232 * t228 + t238, -t240 * t256 + (-qJD(1) * t239 - t223 * t243) * t221 + t238, t223 * t256 + t248, 0, 0; 0, t234 - t250, t234, t223 * t221, 0, 0;];
JaD_transl  = t1;
