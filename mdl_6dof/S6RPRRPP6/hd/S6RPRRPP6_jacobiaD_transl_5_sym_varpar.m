% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:26
% EndTime: 2019-02-26 20:59:26
% DurationCPUTime: 0.20s
% Computational Cost: add. (193->50), mult. (378->76), div. (0->0), fcn. (296->8), ass. (0->43)
t213 = cos(qJ(4));
t245 = t213 * pkin(4);
t205 = pkin(3) + t245;
t211 = sin(qJ(3));
t214 = cos(qJ(3));
t244 = r_i_i_C(3) + qJ(5) + pkin(8);
t249 = t244 * t214;
t256 = (-t205 * t211 - qJ(2) + t249) * qJD(1) - qJD(4) * t245;
t208 = qJ(4) + pkin(9);
t206 = sin(t208);
t207 = cos(t208);
t225 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
t222 = t205 + t225;
t255 = -(-t222 * t211 + t249) * qJD(3) - qJD(5) * t211;
t230 = t244 * t211;
t254 = t222 * t214 + t230;
t215 = cos(qJ(1));
t238 = qJD(4) * t211;
t227 = qJD(1) + t238;
t223 = t227 * t215;
t212 = sin(qJ(1));
t226 = qJD(1) * t211 + qJD(4);
t239 = qJD(3) * t215;
t248 = t226 * t212 - t214 * t239;
t240 = qJD(3) * t214;
t247 = t212 * t240 + t226 * t215;
t210 = sin(qJ(4));
t246 = pkin(4) * t210;
t243 = qJD(1) * t212;
t242 = qJD(1) * t215;
t241 = qJD(3) * t211;
t237 = qJD(4) * t214;
t235 = t214 * qJD(5);
t224 = t227 * t212;
t221 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t246;
t220 = qJD(4) * t221;
t217 = qJD(1) * t254;
t216 = -t238 * t246 - t235 + qJD(2) + (t205 * t214 + t230) * qJD(3) + (-pkin(1) - pkin(7) - t246) * qJD(1);
t204 = -t206 * t224 + t247 * t207;
t203 = -t206 * t247 - t207 * t224;
t202 = -t206 * t223 - t207 * t248;
t201 = t248 * t206 - t207 * t223;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t256 * t212 + t216 * t215, t242, t215 * t217 + (-t221 * t237 - t255) * t212, t203 * r_i_i_C(1) - t204 * r_i_i_C(2) + (-t210 * t247 - t213 * t224) * pkin(4), t212 * t241 - t214 * t242, 0; t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t216 * t212 - t256 * t215, t243, t212 * t217 + (t214 * t220 + t255) * t215, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (-t210 * t248 + t213 * t223) * pkin(4), -t211 * t239 - t214 * t243, 0; 0, 0, -t254 * qJD(3) + t211 * t220 + t235 (-t225 - t245) * t237 + t221 * t241, t240, 0;];
JaD_transl  = t1;
