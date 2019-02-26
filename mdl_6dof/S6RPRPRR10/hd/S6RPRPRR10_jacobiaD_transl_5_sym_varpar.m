% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:00
% EndTime: 2019-02-26 20:54:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (174->44), mult. (311->71), div. (0->0), fcn. (248->8), ass. (0->37)
t209 = sin(qJ(3));
t203 = cos(pkin(10)) * pkin(4) + pkin(3);
t206 = pkin(10) + qJ(5);
t204 = sin(t206);
t205 = cos(t206);
t216 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
t211 = cos(qJ(3));
t236 = r_i_i_C(3) + pkin(8) + qJ(4);
t239 = t236 * t211;
t245 = -(-t216 * t209 + t239) * qJD(3) - qJD(4) * t209;
t224 = t236 * t209;
t244 = t216 * t211 + t224;
t243 = -t203 * t209 - qJ(2) + t239;
t210 = sin(qJ(1));
t220 = qJD(1) * t209 + qJD(5);
t212 = cos(qJ(1));
t231 = qJD(3) * t212;
t238 = t220 * t210 - t211 * t231;
t232 = qJD(3) * t211;
t237 = t210 * t232 + t220 * t212;
t235 = qJD(1) * t210;
t234 = qJD(1) * t212;
t233 = qJD(3) * t209;
t229 = qJD(5) * t211;
t228 = t211 * qJD(4);
t221 = -qJD(5) * t209 - qJD(1);
t219 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
t218 = t221 * t210;
t217 = t221 * t212;
t215 = qJD(5) * t219;
t214 = qJD(1) * t244;
t213 = -t228 + qJD(2) + (t203 * t211 + t224) * qJD(3) + (-sin(pkin(10)) * pkin(4) - pkin(1) - pkin(7)) * qJD(1);
t202 = t204 * t218 + t237 * t205;
t201 = -t237 * t204 + t205 * t218;
t200 = t204 * t217 - t238 * t205;
t199 = t238 * t204 + t205 * t217;
t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t213 * t212 + t243 * t235, t234, t212 * t214 + (-t219 * t229 - t245) * t210, t210 * t233 - t211 * t234, t201 * r_i_i_C(1) - t202 * r_i_i_C(2), 0; t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t213 * t210 - t243 * t234, t235, t210 * t214 + (t211 * t215 + t245) * t212, -t209 * t231 - t211 * t235, -t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0; 0, 0, -t244 * qJD(3) + t209 * t215 + t228, t232 (t204 * t229 + t205 * t233) * r_i_i_C(2) + (t204 * t233 - t205 * t229) * r_i_i_C(1), 0;];
JaD_transl  = t1;
