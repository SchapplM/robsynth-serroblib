% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:54
% EndTime: 2019-02-26 20:40:54
% DurationCPUTime: 0.23s
% Computational Cost: add. (177->42), mult. (292->65), div. (0->0), fcn. (235->7), ass. (0->34)
t189 = pkin(9) + qJ(3);
t187 = sin(t189);
t190 = sin(pkin(10));
t191 = cos(pkin(10));
t215 = r_i_i_C(3) + qJ(5);
t218 = pkin(4) + r_i_i_C(1);
t219 = t215 * t190 + t218 * t191 + pkin(3);
t188 = cos(t189);
t216 = r_i_i_C(2) + qJ(4);
t220 = t216 * t188;
t223 = t219 * t187 - t220;
t206 = qJD(5) * t190;
t199 = t187 * qJD(4) + t188 * t206;
t222 = (-pkin(3) * t187 + t220) * qJD(3) + t199;
t193 = sin(qJ(1));
t214 = t190 * t193;
t194 = cos(qJ(1));
t213 = t190 * t194;
t212 = t193 * t191;
t211 = t194 * t191;
t210 = qJD(1) * t193;
t209 = qJD(1) * t194;
t208 = qJD(3) * t193;
t207 = qJD(3) * t194;
t205 = t187 * t208;
t204 = t187 * t207;
t203 = t216 * t187;
t200 = -t191 * qJD(5) + qJD(2);
t198 = -pkin(3) * t188 - cos(pkin(9)) * pkin(2) - pkin(1) - t203;
t195 = -t187 * t206 + qJD(4) * t188 + (-t188 * t219 - t203) * qJD(3);
t192 = -pkin(7) - qJ(2);
t184 = -t190 * t205 + (t188 * t213 - t212) * qJD(1);
t182 = t190 * t204 + (t188 * t214 + t211) * qJD(1);
t1 = [t200 * t194 + t218 * (t191 * t205 + (-t188 * t211 - t214) * qJD(1)) - t215 * t184 - t222 * t193 + (t193 * t192 + t198 * t194) * qJD(1), t209, t195 * t194 + t223 * t210, -t187 * t210 + t188 * t207, -t182, 0; t200 * t193 + t218 * (-t191 * t204 + (-t188 * t212 + t213) * qJD(1)) - t215 * t182 + t222 * t194 + (-t194 * t192 + t198 * t193) * qJD(1), t210, t195 * t193 - t209 * t223, t187 * t209 + t188 * t208, t184, 0; 0, 0, -qJD(3) * t223 + t199, qJD(3) * t187, qJD(3) * t188 * t190, 0;];
JaD_transl  = t1;
