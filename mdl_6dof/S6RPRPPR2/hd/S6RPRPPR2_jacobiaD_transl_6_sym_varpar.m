% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:47
% EndTime: 2019-02-26 20:39:47
% DurationCPUTime: 0.23s
% Computational Cost: add. (314->47), mult. (348->76), div. (0->0), fcn. (269->10), ass. (0->35)
t219 = qJ(3) + pkin(10);
t215 = sin(t219);
t217 = cos(t219);
t240 = pkin(4) + pkin(8) + r_i_i_C(3);
t231 = t240 * t215 + sin(qJ(3)) * pkin(3);
t258 = (-qJ(5) * t217 + t231) * qJD(3) - t215 * qJD(5);
t222 = sin(qJ(6));
t224 = cos(qJ(6));
t232 = r_i_i_C(1) * t222 + r_i_i_C(2) * t224 + qJ(5);
t256 = -t232 * t217 + t231;
t254 = -t240 * t217 - cos(qJ(3)) * pkin(3);
t236 = qJD(6) * t215 + qJD(1);
t243 = qJD(3) * t224;
t253 = -t217 * t243 + t236 * t222;
t244 = qJD(3) * t222;
t252 = t217 * t244 + t236 * t224;
t249 = pkin(5) + qJ(4) + pkin(7);
t220 = qJ(1) + pkin(9);
t216 = sin(t220);
t247 = qJD(1) * t216;
t218 = cos(t220);
t246 = qJD(1) * t218;
t245 = qJD(3) * t217;
t242 = qJD(6) * t217;
t235 = -qJD(1) * t215 - qJD(6);
t234 = t235 * t222;
t233 = t235 * t224;
t229 = -qJ(5) * t215 - pkin(2) + t254;
t228 = qJD(5) + (r_i_i_C(1) * t224 - r_i_i_C(2) * t222) * qJD(6);
t226 = t228 * t217 + (-t232 * t215 + t254) * qJD(3);
t213 = t216 * t234 + t252 * t218;
t212 = t216 * t233 - t253 * t218;
t211 = -t252 * t216 + t218 * t234;
t210 = t253 * t216 + t218 * t233;
t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t218 * qJD(4) + t258 * t216 + (-cos(qJ(1)) * pkin(1) - t249 * t216 + t229 * t218) * qJD(1), 0, t226 * t218 + t256 * t247, t246, -t215 * t247 + t218 * t245, t212 * r_i_i_C(1) - t213 * r_i_i_C(2); t213 * r_i_i_C(1) + t212 * r_i_i_C(2) + t216 * qJD(4) - t258 * t218 + (-sin(qJ(1)) * pkin(1) + t249 * t218 + t229 * t216) * qJD(1), 0, t226 * t216 - t246 * t256, t247, t215 * t246 + t216 * t245, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2); 0, 0, -qJD(3) * t256 + t228 * t215, 0, qJD(3) * t215 (-t215 * t244 + t224 * t242) * r_i_i_C(2) + (t215 * t243 + t222 * t242) * r_i_i_C(1);];
JaD_transl  = t1;
