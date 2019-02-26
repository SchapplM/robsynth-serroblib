% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:08
% EndTime: 2019-02-26 20:43:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (126->47), mult. (382->75), div. (0->0), fcn. (292->6), ass. (0->36)
t203 = sin(qJ(3));
t206 = cos(qJ(3));
t220 = pkin(3) + pkin(4) + pkin(8) + r_i_i_C(3);
t227 = pkin(5) + qJ(4);
t233 = -qJD(5) + (t220 * t203 - t227 * t206 + qJ(2)) * qJD(1);
t205 = cos(qJ(6));
t217 = qJD(6) * t206 + qJD(1);
t232 = t205 * t217;
t202 = sin(qJ(6));
t231 = t217 * t202;
t230 = t220 * t206;
t210 = -(r_i_i_C(1) * t202 + r_i_i_C(2) * t205) * qJD(6) + qJD(4);
t229 = -t220 * qJD(3) + t210;
t226 = qJD(1) * t206;
t207 = cos(qJ(1));
t225 = qJD(1) * t207;
t204 = sin(qJ(1));
t224 = qJD(3) * t204;
t223 = qJD(3) * t206;
t222 = qJD(3) * t207;
t221 = qJD(6) * t203;
t219 = t203 * t224;
t218 = t203 * t222;
t216 = qJD(6) + t226;
t215 = qJD(1) * t220;
t213 = t216 * t207;
t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t202 + t227;
t211 = qJD(1) * t212;
t209 = t216 * t204 + t218;
t208 = -t206 * qJD(4) + qJD(2) + (-pkin(1) - pkin(7) + qJ(5)) * qJD(1) + (t227 * t203 + t230) * qJD(3);
t201 = qJD(1) * t204;
t200 = t205 * t213 + (-qJD(3) * t203 * t205 - t231) * t204;
t199 = t204 * t232 + (t213 - t219) * t202;
t198 = t209 * t205 + t207 * t231;
t197 = t209 * t202 - t207 * t232;
t1 = [t198 * r_i_i_C(1) - t197 * r_i_i_C(2) - t233 * t204 + t208 * t207, t225 (t207 * t215 + t212 * t224) * t206 + (t229 * t204 + t207 * t211) * t203, -t206 * t225 + t219, t201, t199 * r_i_i_C(1) + t200 * r_i_i_C(2); -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t208 * t204 + t233 * t207, t201 (t204 * t215 - t212 * t222) * t206 + (t204 * t211 - t229 * t207) * t203, -t204 * t226 - t218, -t225, t197 * r_i_i_C(1) + t198 * r_i_i_C(2); 0, 0, t210 * t206 + (-t212 * t203 - t230) * qJD(3), t223, 0 (t202 * t221 - t205 * t223) * r_i_i_C(2) + (-t202 * t223 - t205 * t221) * r_i_i_C(1);];
JaD_transl  = t1;
