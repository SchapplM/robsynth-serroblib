% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:57
% EndTime: 2019-02-26 22:03:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (288->48), mult. (249->53), div. (0->0), fcn. (170->8), ass. (0->42)
t201 = qJ(2) + qJ(3);
t196 = pkin(10) + t201;
t195 = cos(t196);
t230 = r_i_i_C(3) + qJ(5);
t215 = t230 * t195;
t194 = sin(t196);
t193 = t194 * qJD(5);
t197 = sin(t201);
t200 = qJD(2) + qJD(3);
t202 = sin(qJ(2));
t229 = pkin(2) * qJD(2);
t223 = t202 * t229;
t233 = pkin(3) * t200;
t236 = pkin(4) - r_i_i_C(2);
t241 = (-t236 * t194 + t215) * t200 + (r_i_i_C(1) + qJ(4) + pkin(8) + pkin(7)) * qJD(1) - t197 * t233 + t193 - t223;
t203 = sin(qJ(1));
t205 = cos(qJ(1));
t225 = qJD(1) * t205;
t228 = t195 * t200;
t240 = t194 * t225 + t203 * t228;
t235 = pkin(3) * t197;
t198 = cos(t201);
t234 = pkin(3) * t198;
t232 = pkin(4) * t194;
t227 = t200 * t194;
t226 = qJD(1) * t203;
t224 = qJD(5) * t195;
t221 = t205 * t228;
t218 = t194 * t226;
t220 = pkin(4) * t218 + r_i_i_C(2) * t221 + t205 * t224;
t216 = t230 * t194;
t213 = -t232 - t235;
t212 = t240 * r_i_i_C(2) + t203 * t224 + t225 * t215;
t211 = -pkin(4) * t195 - t216;
t210 = -r_i_i_C(2) * t194 - t215;
t204 = cos(qJ(2));
t209 = -t198 * t233 + t211 * t200 - t204 * t229;
t208 = t200 * (t211 - t234);
t207 = r_i_i_C(2) * t227 + t213 * t200 + t230 * t228 + t193;
t206 = qJD(4) + (-t204 * pkin(2) - t236 * t195 - pkin(1) - t216 - t234) * qJD(1);
t190 = -t202 * pkin(2) - t235;
t1 = [-t241 * t203 + t206 * t205, t209 * t205 + (-t190 + t210) * t226 + t220, t205 * t208 + (t210 + t235) * t226 + t220, t225, -t218 + t221, 0; t206 * t203 + t241 * t205 (t190 - t232) * t225 + t209 * t203 + t212, t203 * t208 + t213 * t225 + t212, t226, t240, 0; 0, t207 - t223, t207, 0, t227, 0;];
JaD_transl  = t1;
