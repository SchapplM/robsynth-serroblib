% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:36
% EndTime: 2019-02-26 22:10:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (217->43), mult. (285->60), div. (0->0), fcn. (208->8), ass. (0->38)
t224 = qJ(2) + qJ(3);
t222 = cos(t224);
t259 = r_i_i_C(3) + qJ(4);
t270 = t259 * t222;
t221 = sin(t224);
t219 = t221 * qJD(4);
t223 = qJD(2) + qJD(3);
t225 = sin(pkin(10));
t226 = cos(pkin(10));
t260 = r_i_i_C(2) * t225;
t268 = r_i_i_C(1) * t226 + pkin(3);
t237 = t268 - t260;
t227 = sin(qJ(2));
t258 = pkin(2) * qJD(2);
t250 = t227 * t258;
t269 = (-t237 * t221 + t270) * t223 + (t225 * r_i_i_C(1) + t226 * r_i_i_C(2) + pkin(7) + pkin(8)) * qJD(1) + t219 - t250;
t251 = t223 * t260;
t267 = (qJD(4) + t251) * t222;
t262 = pkin(2) * t227;
t228 = sin(qJ(1));
t256 = t222 * t228;
t230 = cos(qJ(1));
t255 = t223 * t230;
t254 = qJD(1) * t228;
t253 = qJD(1) * t230;
t248 = t221 * t254;
t247 = t221 * t253;
t245 = t259 * t221;
t243 = t259 * t228;
t241 = t267 * t230 + t268 * t248;
t240 = t268 * t223;
t239 = t268 * t230;
t238 = t267 * t228 + t247 * t260 + t253 * t270;
t234 = t219 + t223 * t270 + (-t240 + t251) * t221;
t229 = cos(qJ(2));
t233 = qJD(1) * (-t229 * pkin(2) - t237 * t222 - pkin(1) - t245);
t232 = -t229 * t258 + (-t222 * t268 - t245) * t223;
t1 = [-t228 * t269 + t230 * t233 (-t221 * t260 + t262 - t270) * t254 + t232 * t230 + t241 (-t254 * t260 - t259 * t255) * t221 + (-qJD(1) * t243 - t223 * t239) * t222 + t241, t222 * t255 - t248, 0, 0; t228 * t233 + t230 * t269 (-t221 * t268 - t262) * t253 + t232 * t228 + t238, -t240 * t256 + (-qJD(1) * t239 - t223 * t243) * t221 + t238, t223 * t256 + t247, 0, 0; 0, t234 - t250, t234, t223 * t221, 0, 0;];
JaD_transl  = t1;
