% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:25
% EndTime: 2019-02-26 20:40:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (243->47), mult. (378->72), div. (0->0), fcn. (288->8), ass. (0->35)
t206 = sin(qJ(3));
t208 = cos(qJ(3));
t231 = pkin(5) + qJ(4);
t222 = pkin(3) + pkin(4) + pkin(8) + r_i_i_C(3);
t234 = t222 * t206;
t239 = (-t231 * t208 + t234) * qJD(3) - t206 * qJD(4);
t205 = sin(qJ(6));
t207 = cos(qJ(6));
t214 = r_i_i_C(1) * t207 - r_i_i_C(2) * t205 + t231;
t237 = -t214 * t208 + t234;
t227 = qJD(1) * t206;
t219 = qJD(6) + t227;
t236 = t207 * t219;
t220 = qJD(6) * t206 + qJD(1);
t225 = qJD(3) * t208;
t232 = t205 * t225 + t220 * t207;
t230 = pkin(7) - qJ(5);
t204 = qJ(1) + pkin(9);
t202 = sin(t204);
t229 = qJD(1) * t202;
t203 = cos(t204);
t228 = qJD(1) * t203;
t226 = qJD(3) * t206;
t224 = qJD(6) * t208;
t217 = t222 * t208;
t215 = t219 * t205;
t213 = qJD(4) + (-r_i_i_C(1) * t205 - r_i_i_C(2) * t207) * qJD(6);
t212 = -t231 * t206 - pkin(2) - t217;
t211 = t220 * t205 - t207 * t225;
t209 = t213 * t208 + (-t214 * t206 - t217) * qJD(3);
t201 = t211 * t202 - t203 * t236;
t200 = t232 * t202 + t203 * t215;
t199 = t202 * t236 + t211 * t203;
t198 = t202 * t215 - t232 * t203;
t1 = [t201 * r_i_i_C(1) + t200 * r_i_i_C(2) - t203 * qJD(5) + t239 * t202 + (-cos(qJ(1)) * pkin(1) - t230 * t202 + t212 * t203) * qJD(1), 0, t209 * t203 + t237 * t229, -t202 * t227 + t203 * t225, -t228, t198 * r_i_i_C(1) + t199 * r_i_i_C(2); -t199 * r_i_i_C(1) + t198 * r_i_i_C(2) - t202 * qJD(5) - t239 * t203 + (-sin(qJ(1)) * pkin(1) + t230 * t203 + t212 * t202) * qJD(1), 0, t209 * t202 - t228 * t237, t202 * t225 + t203 * t227, -t229, -t200 * r_i_i_C(1) + t201 * r_i_i_C(2); 0, 0, -qJD(3) * t237 + t213 * t206, t226, 0 (-t205 * t224 - t207 * t226) * r_i_i_C(2) + (-t205 * t226 + t207 * t224) * r_i_i_C(1);];
JaD_transl  = t1;
