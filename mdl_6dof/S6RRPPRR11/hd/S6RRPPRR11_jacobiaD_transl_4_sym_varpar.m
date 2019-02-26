% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR11_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:20
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (137->38), mult. (400->56), div. (0->0), fcn. (372->8), ass. (0->29)
t203 = sin(pkin(11));
t204 = sin(pkin(6));
t205 = cos(pkin(11));
t223 = t204 * (r_i_i_C(1) * t205 - r_i_i_C(2) * t203 + pkin(3) + pkin(8));
t207 = sin(qJ(2));
t208 = sin(qJ(1));
t222 = t208 * t207;
t209 = cos(qJ(2));
t221 = t208 * t209;
t210 = cos(qJ(1));
t220 = t210 * t207;
t219 = t210 * t209;
t218 = qJD(2) * t207;
t217 = qJD(2) * t209;
t216 = r_i_i_C(3) + qJ(4) + pkin(2);
t206 = cos(pkin(6));
t215 = t206 * t222;
t214 = t206 * t219;
t213 = qJD(2) * t206 + qJD(1);
t212 = t203 * r_i_i_C(1) + t205 * r_i_i_C(2) + qJ(3);
t197 = t206 * t221 + t220;
t196 = t206 * t220 + t221;
t198 = t215 - t219;
t195 = -t214 + t222;
t194 = -qJD(1) * t215 - t208 * t218 + t213 * t219;
t193 = t197 * qJD(1) + t196 * qJD(2);
t192 = t196 * qJD(1) + t197 * qJD(2);
t191 = -qJD(1) * t214 - t210 * t217 + t213 * t222;
t1 = [-t195 * qJD(3) - t196 * qJD(4) - t216 * t194 - t212 * t193 + (-t210 * pkin(1) - t208 * t223) * qJD(1), -t198 * qJD(3) - t197 * qJD(4) + t216 * t191 - t212 * t192, -t191, -t192, 0, 0; t197 * qJD(3) - t198 * qJD(4) - t216 * t192 - t212 * t191 + (-t208 * pkin(1) + t210 * t223) * qJD(1), t196 * qJD(3) - t195 * qJD(4) - t216 * t193 + t212 * t194, t193, t194, 0, 0; 0 (qJD(3) * t207 + qJD(4) * t209 + (-t216 * t207 + t212 * t209) * qJD(2)) * t204, t204 * t218, t204 * t217, 0, 0;];
JaD_transl  = t1;
