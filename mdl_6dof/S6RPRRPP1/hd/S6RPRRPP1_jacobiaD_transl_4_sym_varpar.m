% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:18
% EndTime: 2019-02-26 20:56:18
% DurationCPUTime: 0.15s
% Computational Cost: add. (165->37), mult. (274->63), div. (0->0), fcn. (213->8), ass. (0->32)
t207 = sin(qJ(3));
t206 = sin(qJ(4));
t208 = cos(qJ(4));
t215 = r_i_i_C(1) * t208 - r_i_i_C(2) * t206 + pkin(3);
t209 = cos(qJ(3));
t228 = pkin(8) + r_i_i_C(3);
t230 = t228 * t209;
t211 = -t215 * t207 + t230;
t235 = qJD(1) * t211;
t234 = (-pkin(3) * t207 + t230) * qJD(3);
t217 = qJD(1) * t209 - qJD(4);
t232 = t208 * t217;
t223 = qJD(4) * t209;
t218 = -qJD(1) + t223;
t226 = qJD(3) * t207;
t229 = -t206 * t226 + t218 * t208;
t225 = qJD(3) * t209;
t224 = qJD(4) * t207;
t222 = t228 * t207;
t216 = r_i_i_C(1) * t206 + r_i_i_C(2) * t208;
t214 = t217 * t206;
t213 = -pkin(3) * t209 - pkin(2) - t222;
t212 = t218 * t206 + t208 * t226;
t210 = t216 * t224 + (-t215 * t209 - t222) * qJD(3);
t205 = qJ(1) + pkin(9);
t204 = cos(t205);
t203 = sin(t205);
t202 = t212 * t203 - t204 * t232;
t201 = t229 * t203 + t204 * t214;
t200 = t203 * t232 + t212 * t204;
t199 = t203 * t214 - t229 * t204;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, 0, t211 * qJD(3) - t216 * t223 (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
