% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:33
% EndTime: 2019-02-26 22:04:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (179->35), mult. (217->48), div. (0->0), fcn. (148->6), ass. (0->34)
t183 = qJ(2) + qJ(3);
t181 = cos(t183);
t211 = r_i_i_C(3) + qJ(4);
t195 = t211 * t181;
t180 = sin(t183);
t178 = t180 * qJD(4);
t182 = qJD(2) + qJD(3);
t214 = pkin(3) + r_i_i_C(1);
t203 = t214 * t180;
t184 = sin(qJ(2));
t210 = pkin(2) * qJD(2);
t204 = t184 * t210;
t217 = (-t203 + t195) * t182 + (r_i_i_C(2) + pkin(8) + pkin(7)) * qJD(1) + t178 - t204;
t213 = pkin(2) * t184;
t209 = t181 * t182;
t187 = cos(qJ(1));
t208 = t182 * t187;
t185 = sin(qJ(1));
t207 = qJD(1) * t185;
t206 = qJD(1) * t187;
t205 = qJD(4) * t181;
t202 = t214 * t187;
t201 = t185 * t209;
t200 = t185 * t205 + t206 * t195;
t197 = t180 * t207;
t199 = t187 * t205 + t214 * t197;
t196 = t211 * t180;
t194 = t211 * t185;
t192 = -t214 * t181 - t196;
t191 = -t182 * t203 + t211 * t209 + t178;
t186 = cos(qJ(2));
t190 = qJD(1) * (-t186 * pkin(2) - pkin(1) + t192);
t189 = t192 * t182 - t186 * t210;
t1 = [-t217 * t185 + t187 * t190 (-t195 + t213) * t207 + t189 * t187 + t199, -t196 * t208 + (-qJD(1) * t194 - t182 * t202) * t181 + t199, t181 * t208 - t197, 0, 0; t185 * t190 + t217 * t187 (-t203 - t213) * t206 + t189 * t185 + t200, -t214 * t201 + (-qJD(1) * t202 - t182 * t194) * t180 + t200, t180 * t206 + t201, 0, 0; 0, t191 - t204, t191, t182 * t180, 0, 0;];
JaD_transl  = t1;
