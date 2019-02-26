% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:48
% EndTime: 2019-02-26 22:31:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (411->52), mult. (306->58), div. (0->0), fcn. (207->8), ass. (0->45)
t203 = qJ(2) + qJ(3);
t200 = qJ(4) + t203;
t196 = cos(t200);
t235 = r_i_i_C(3) + qJ(5);
t217 = t235 * t196;
t195 = sin(t200);
t194 = t195 * qJD(5);
t201 = qJD(2) + qJD(3);
t197 = qJD(4) + t201;
t204 = sin(qJ(2));
t234 = pkin(2) * qJD(2);
t225 = t204 * t234;
t198 = sin(t203);
t238 = pkin(3) * t201;
t227 = t198 * t238;
t240 = pkin(4) - r_i_i_C(2);
t245 = (-t240 * t195 + t217) * t197 + (r_i_i_C(1) + pkin(9) + pkin(8) + pkin(7)) * qJD(1) + t194 - t225 - t227;
t205 = sin(qJ(1));
t231 = t197 * t205;
t224 = t196 * t231;
t207 = cos(qJ(1));
t229 = qJD(1) * t207;
t244 = t195 * t229 + t224;
t199 = cos(t203);
t218 = t235 * t195;
t212 = (-pkin(4) * t196 - t218) * t197;
t208 = -t199 * t238 + t212;
t239 = pkin(3) * t198;
t237 = pkin(4) * t195;
t233 = t196 * t197;
t232 = t197 * t195;
t230 = qJD(1) * t205;
t228 = qJD(5) * t196;
t223 = t207 * t233;
t220 = t195 * t230;
t222 = pkin(4) * t220 + r_i_i_C(2) * t223 + t207 * t228;
t215 = t244 * r_i_i_C(2) + t205 * t228 + t229 * t217;
t214 = -r_i_i_C(2) * t195 - t217;
t213 = -t240 * t232 + t235 * t233 + t194;
t206 = cos(qJ(2));
t211 = -t206 * t234 + t208;
t210 = qJD(1) * (-t206 * pkin(2) - pkin(3) * t199 - t240 * t196 - pkin(1) - t218);
t209 = t213 - t227;
t191 = -t204 * pkin(2) - t239;
t1 = [-t245 * t205 + t207 * t210, t211 * t207 + (-t191 + t214) * t230 + t222, t208 * t207 + (t214 + t239) * t230 + t222, t207 * t212 + t214 * t230 + t222, -t220 + t223, 0; t205 * t210 + t245 * t207 (t191 - t237) * t229 + t211 * t205 + t215 (-t237 - t239) * t229 + t208 * t205 + t215, -pkin(4) * t224 + (-pkin(4) * t229 - t235 * t231) * t195 + t215, t244, 0; 0, t209 - t225, t209, t213, t232, 0;];
JaD_transl  = t1;
