% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_3_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_3_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_3_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:21
% EndTime: 2019-02-26 22:54:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (70->32), mult. (236->67), div. (0->0), fcn. (190->6), ass. (0->28)
t217 = pkin(2) + r_i_i_C(3);
t196 = cos(qJ(3));
t198 = cos(qJ(1));
t216 = t196 * t198;
t195 = sin(qJ(1));
t215 = qJD(1) * t195;
t214 = qJD(1) * t198;
t194 = sin(qJ(2));
t213 = qJD(2) * t194;
t197 = cos(qJ(2));
t212 = qJD(2) * t197;
t211 = qJD(2) * t198;
t210 = qJD(3) * t194;
t209 = qJD(3) * t197;
t208 = qJD(1) + t209;
t207 = qJD(1) * t197 + qJD(3);
t193 = sin(qJ(3));
t206 = r_i_i_C(1) * t196 - r_i_i_C(2) * t193;
t205 = r_i_i_C(1) * t193 + r_i_i_C(2) * t196;
t204 = t208 * t193;
t203 = qJD(2) * t206;
t200 = t194 * t211 + t207 * t195;
t199 = t217 * qJD(2) + t205 * qJD(3);
t192 = -t207 * t216 + (t196 * t213 + t204) * t195;
t191 = t208 * t196 * t195 + (-t195 * t213 + t207 * t198) * t193;
t190 = t200 * t196 + t198 * t204;
t189 = t200 * t193 - t208 * t216;
t1 = [t192 * r_i_i_C(1) + t191 * r_i_i_C(2) + t217 * (t194 * t214 + t195 * t212) (-t198 * t203 + t217 * t215) * t197 + (t199 * t198 + t206 * t215) * t194, t189 * r_i_i_C(1) + t190 * r_i_i_C(2), 0, 0, 0, 0; -t190 * r_i_i_C(1) + t189 * r_i_i_C(2) + t217 * (t194 * t215 - t197 * t211) (-t195 * t203 - t217 * t214) * t197 + (t199 * t195 - t206 * t214) * t194, -t191 * r_i_i_C(1) + t192 * r_i_i_C(2), 0, 0, 0, 0; 0, -t205 * t209 + (-t206 * t194 - t217 * t197) * qJD(2) (t193 * t210 - t196 * t212) * r_i_i_C(2) + (-t193 * t212 - t196 * t210) * r_i_i_C(1), 0, 0, 0, 0;];
JaD_transl  = t1;
