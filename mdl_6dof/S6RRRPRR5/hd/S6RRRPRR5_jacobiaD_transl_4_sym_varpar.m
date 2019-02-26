% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:22
% EndTime: 2019-02-26 22:18:22
% DurationCPUTime: 0.16s
% Computational Cost: add. (179->38), mult. (217->49), div. (0->0), fcn. (148->6), ass. (0->35)
t189 = qJ(2) + qJ(3);
t187 = cos(t189);
t218 = r_i_i_C(3) + qJ(4);
t202 = t218 * t187;
t186 = sin(t189);
t184 = t186 * qJD(4);
t188 = qJD(2) + qJD(3);
t190 = sin(qJ(2));
t217 = pkin(2) * qJD(2);
t210 = t190 * t217;
t221 = pkin(3) - r_i_i_C(2);
t226 = (-t221 * t186 + t202) * t188 + (r_i_i_C(1) + pkin(8) + pkin(7)) * qJD(1) + t184 - t210;
t191 = sin(qJ(1));
t214 = t188 * t191;
t209 = t187 * t214;
t193 = cos(qJ(1));
t212 = qJD(1) * t193;
t225 = t186 * t212 + t209;
t220 = pkin(2) * t190;
t216 = t187 * t188;
t215 = t188 * t186;
t213 = qJD(1) * t191;
t211 = qJD(4) * t187;
t208 = t193 * t216;
t205 = t186 * t213;
t207 = pkin(3) * t205 + r_i_i_C(2) * t208 + t193 * t211;
t203 = t218 * t186;
t201 = t225 * r_i_i_C(2) + t191 * t211 + t212 * t202;
t200 = -r_i_i_C(2) * t186 - t202;
t198 = -t221 * t215 + t218 * t216 + t184;
t197 = (-pkin(3) * t187 - t203) * t188;
t192 = cos(qJ(2));
t196 = qJD(1) * (-t192 * pkin(2) - t221 * t187 - pkin(1) - t203);
t195 = -t192 * t217 + t197;
t1 = [-t226 * t191 + t193 * t196, t195 * t193 + (t200 + t220) * t213 + t207, t193 * t197 + t200 * t213 + t207, -t205 + t208, 0, 0; t191 * t196 + t226 * t193 (-pkin(3) * t186 - t220) * t212 + t195 * t191 + t201, -pkin(3) * t209 + (-pkin(3) * t212 - t218 * t214) * t186 + t201, t225, 0, 0; 0, t198 - t210, t198, t215, 0, 0;];
JaD_transl  = t1;
