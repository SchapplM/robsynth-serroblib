% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (209->47), mult. (356->70), div. (0->0), fcn. (275->8), ass. (0->35)
t211 = qJ(3) + pkin(9);
t208 = sin(t211);
t209 = cos(t211);
t213 = sin(qJ(6));
t216 = cos(qJ(6));
t222 = qJD(5) + (r_i_i_C(1) * t216 - r_i_i_C(2) * t213) * qJD(6);
t233 = pkin(4) + pkin(8) + r_i_i_C(3);
t225 = t233 * t208 + sin(qJ(3)) * pkin(3);
t226 = r_i_i_C(1) * t213 + r_i_i_C(2) * t216 + qJ(5);
t250 = -t222 * t208 + (-t226 * t209 + t225) * qJD(3);
t249 = qJD(4) + (-qJ(5) * t209 + qJ(2) + t225) * qJD(1);
t224 = t233 * t209 + cos(qJ(3)) * pkin(3);
t243 = t226 * t208 + t224;
t218 = cos(qJ(1));
t239 = t216 * t218;
t215 = sin(qJ(1));
t238 = qJD(1) * t215;
t210 = qJD(1) * t218;
t237 = qJD(3) * t208;
t236 = qJD(3) * t209;
t235 = qJD(3) * t216;
t234 = qJD(6) * t208;
t232 = t215 * t237;
t231 = t218 * t237;
t230 = qJD(6) * t209 + qJD(1);
t229 = qJD(1) * t209 + qJD(6);
t227 = t230 * t213;
t221 = t229 * t215 + t231;
t220 = qJD(1) * t243;
t219 = -qJD(5) * t209 + qJD(2) + (-pkin(1) - pkin(5) - qJ(4) - pkin(7)) * qJD(1) + (qJ(5) * t208 + t224) * qJD(3);
t207 = t221 * t213 - t230 * t239;
t206 = t221 * t216 + t218 * t227;
t205 = t230 * t216 * t215 + (t229 * t218 - t232) * t213;
t204 = -t229 * t239 + (t208 * t235 + t227) * t215;
t1 = [t207 * r_i_i_C(1) + t206 * r_i_i_C(2) - t249 * t215 + t219 * t218, t210, -t250 * t215 + t218 * t220, -t238, -t209 * t210 + t232, r_i_i_C(1) * t204 + r_i_i_C(2) * t205; -t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t219 * t215 + t249 * t218, t238, t215 * t220 + t250 * t218, t210, -t209 * t238 - t231, -r_i_i_C(1) * t206 + r_i_i_C(2) * t207; 0, 0, -t243 * qJD(3) + t222 * t209, 0, t236 (-t213 * t236 - t216 * t234) * r_i_i_C(2) + (t209 * t235 - t213 * t234) * r_i_i_C(1);];
JaD_transl  = t1;
