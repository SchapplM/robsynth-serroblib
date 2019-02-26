% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:43
% EndTime: 2019-02-26 21:33:43
% DurationCPUTime: 0.23s
% Computational Cost: add. (188->46), mult. (373->73), div. (0->0), fcn. (294->8), ass. (0->37)
t208 = sin(qJ(2));
t210 = cos(qJ(2));
t225 = sin(pkin(10)) * pkin(4) + qJ(3);
t229 = t210 * qJD(4);
t228 = pkin(2) + r_i_i_C(3) + pkin(8) + qJ(4);
t238 = t228 * t208;
t245 = (-t225 * t210 + t238) * qJD(2) - (pkin(7) + cos(pkin(10)) * pkin(4) + pkin(3)) * qJD(1) - t208 * qJD(3) - t229;
t205 = pkin(10) + qJ(5);
t203 = sin(t205);
t204 = cos(t205);
t218 = r_i_i_C(1) * t203 + r_i_i_C(2) * t204 + t225;
t243 = -t218 * t210 + t238;
t209 = sin(qJ(1));
t224 = qJD(5) * t208 + qJD(1);
t242 = t209 * t224;
t211 = cos(qJ(1));
t241 = t211 * t224;
t235 = qJD(1) * t209;
t234 = qJD(1) * t211;
t233 = qJD(2) * t208;
t232 = qJD(2) * t210;
t231 = qJD(2) * t211;
t230 = qJD(5) * t210;
t227 = t209 * t232;
t226 = t210 * t231;
t223 = -qJD(1) * t208 - qJD(5);
t221 = t228 * t210;
t217 = qJD(3) + (r_i_i_C(1) * t204 - r_i_i_C(2) * t203) * qJD(5);
t216 = t223 * t211 - t227;
t215 = t223 * t209 + t226;
t214 = qJD(1) * (-t225 * t208 - pkin(1) - t221);
t212 = -qJD(4) * t208 + t217 * t210 + (-t218 * t208 - t221) * qJD(2);
t201 = t215 * t203 + t204 * t241;
t200 = -t203 * t241 + t215 * t204;
t199 = t216 * t203 - t204 * t242;
t198 = t203 * t242 + t216 * t204;
t1 = [t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t245 * t209 + t211 * t214, t212 * t211 + t243 * t235, -t208 * t235 + t226, -t208 * t231 - t210 * t235, t200 * r_i_i_C(1) - t201 * r_i_i_C(2), 0; t201 * r_i_i_C(1) + t200 * r_i_i_C(2) + t209 * t214 - t245 * t211, t212 * t209 - t234 * t243, t208 * t234 + t227, -t209 * t233 + t210 * t234, -t198 * r_i_i_C(1) + t199 * r_i_i_C(2), 0; 0, -qJD(2) * t243 + t217 * t208 + t229, t233, t232 (-t203 * t233 + t204 * t230) * r_i_i_C(2) + (t203 * t230 + t204 * t233) * r_i_i_C(1), 0;];
JaD_transl  = t1;
