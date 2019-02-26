% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiaD_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:10
% EndTime: 2019-02-26 21:26:10
% DurationCPUTime: 0.23s
% Computational Cost: add. (123->48), mult. (374->76), div. (0->0), fcn. (286->6), ass. (0->34)
t203 = sin(qJ(2));
t206 = cos(qJ(2));
t220 = pkin(2) + pkin(3) + pkin(8) + r_i_i_C(3);
t215 = t220 * t203;
t230 = pkin(4) + qJ(3);
t233 = (-t230 * t206 + t215) * qJD(2) - t203 * qJD(3);
t229 = pkin(7) - qJ(4);
t205 = cos(qJ(5));
t207 = cos(qJ(1));
t228 = t205 * t207;
t204 = sin(qJ(1));
t227 = qJD(1) * t204;
t226 = qJD(1) * t207;
t225 = qJD(2) * t203;
t224 = qJD(2) * t206;
t223 = qJD(2) * t207;
t222 = qJD(5) * t206;
t219 = t204 * t224;
t218 = t206 * t223;
t217 = qJD(5) * t203 + qJD(1);
t216 = qJD(1) * t203 + qJD(5);
t214 = t220 * t206;
t202 = sin(qJ(5));
t213 = t217 * t202;
t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t202 + t230;
t211 = qJD(3) + (-r_i_i_C(1) * t202 - r_i_i_C(2) * t205) * qJD(5);
t210 = -t230 * t203 - pkin(1) - t214;
t209 = t216 * t204 - t218;
t208 = t212 * t206 - t215;
t201 = -t216 * t228 + (-t205 * t224 + t213) * t204;
t200 = t217 * t205 * t204 + (t216 * t207 + t219) * t202;
t199 = t209 * t205 + t207 * t213;
t198 = t209 * t202 - t217 * t228;
t1 = [t201 * r_i_i_C(1) + t200 * r_i_i_C(2) - qJD(4) * t207 + t233 * t204 + (-t229 * t204 + t210 * t207) * qJD(1) (-t212 * t223 + t220 * t227) * t203 + (-t212 * t227 + (-t220 * qJD(2) + t211) * t207) * t206, -t203 * t227 + t218, -t226, r_i_i_C(1) * t198 + r_i_i_C(2) * t199, 0; -t199 * r_i_i_C(1) + t198 * r_i_i_C(2) - t204 * qJD(4) - t233 * t207 + (t210 * t204 + t229 * t207) * qJD(1), t208 * t226 + (t211 * t206 + (-t212 * t203 - t214) * qJD(2)) * t204, t203 * t226 + t219, -t227, -r_i_i_C(1) * t200 + r_i_i_C(2) * t201, 0; 0, t208 * qJD(2) + t211 * t203, t225, 0 (-t202 * t222 - t205 * t225) * r_i_i_C(2) + (-t202 * t225 + t205 * t222) * r_i_i_C(1), 0;];
JaD_transl  = t1;
