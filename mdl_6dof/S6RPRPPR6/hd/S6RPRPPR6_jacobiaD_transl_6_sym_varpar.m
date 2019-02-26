% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:59
% EndTime: 2019-02-26 20:42:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (269->50), mult. (341->75), div. (0->0), fcn. (269->10), ass. (0->38)
t221 = qJ(3) + pkin(9);
t216 = sin(t221);
t218 = cos(t221);
t251 = r_i_i_C(3) + pkin(8) + qJ(5);
t234 = t251 * t218 - sin(qJ(3)) * pkin(3);
t214 = cos(pkin(10)) * pkin(5) + pkin(4);
t220 = pkin(10) + qJ(6);
t215 = sin(t220);
t217 = cos(t220);
t235 = r_i_i_C(1) * t217 - r_i_i_C(2) * t215 + t214;
t263 = -(-t235 * t216 + t234) * qJD(3) - qJD(5) * t216;
t262 = -qJD(4) + (-t214 * t216 - qJ(2) + t234) * qJD(1);
t232 = t251 * t216 + cos(qJ(3)) * pkin(3);
t259 = t235 * t218 + t232;
t226 = sin(qJ(1));
t239 = qJD(1) * t216 + qJD(6);
t228 = cos(qJ(1));
t246 = qJD(3) * t228;
t255 = -t218 * t246 + t239 * t226;
t247 = qJD(3) * t226;
t254 = t218 * t247 + t239 * t228;
t249 = qJD(1) * t226;
t219 = qJD(1) * t228;
t248 = qJD(3) * t216;
t244 = qJD(6) * t218;
t243 = t218 * qJD(5);
t240 = -qJD(6) * t216 - qJD(1);
t238 = r_i_i_C(1) * t215 + r_i_i_C(2) * t217;
t237 = t240 * t226;
t236 = t240 * t228;
t231 = qJD(6) * t238;
t230 = qJD(1) * t259;
t229 = -t243 + qJD(2) + (-pkin(5) * sin(pkin(10)) - pkin(1) - qJ(4) - pkin(7)) * qJD(1) + (t214 * t218 + t232) * qJD(3);
t213 = t215 * t237 + t254 * t217;
t212 = -t254 * t215 + t217 * t237;
t211 = t215 * t236 - t255 * t217;
t210 = t255 * t215 + t217 * t236;
t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t262 * t226 + t229 * t228, t219, t228 * t230 + (-t238 * t244 - t263) * t226, -t249, t216 * t247 - t218 * t219, t212 * r_i_i_C(1) - t213 * r_i_i_C(2); t213 * r_i_i_C(1) + t212 * r_i_i_C(2) + t229 * t226 - t262 * t228, t249, t226 * t230 + (t218 * t231 + t263) * t228, t219, -t216 * t246 - t218 * t249, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2); 0, 0, -t259 * qJD(3) + t216 * t231 + t243, 0, qJD(3) * t218 (t215 * t244 + t217 * t248) * r_i_i_C(2) + (t215 * t248 - t217 * t244) * r_i_i_C(1);];
JaD_transl  = t1;
