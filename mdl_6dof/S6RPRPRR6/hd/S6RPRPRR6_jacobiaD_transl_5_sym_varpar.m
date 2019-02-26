% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:43
% EndTime: 2019-02-26 20:51:43
% DurationCPUTime: 0.17s
% Computational Cost: add. (257->47), mult. (309->75), div. (0->0), fcn. (248->9), ass. (0->38)
t212 = cos(pkin(11)) * pkin(4) + pkin(3);
t219 = pkin(10) + qJ(3);
t215 = sin(t219);
t239 = t215 * qJD(4);
t217 = cos(t219);
t248 = r_i_i_C(3) + pkin(8) + qJ(4);
t250 = t248 * t217;
t254 = (-t212 * t215 + t250) * qJD(3) + t239;
t218 = pkin(11) + qJ(5);
t214 = sin(t218);
t216 = cos(t218);
t229 = r_i_i_C(1) * t216 - r_i_i_C(2) * t214 + t212;
t226 = -t229 * t215 + t250;
t224 = cos(qJ(1));
t240 = qJD(5) * t217;
t233 = -qJD(1) + t240;
t252 = t224 * t233;
t232 = qJD(1) * t217 - qJD(5);
t223 = sin(qJ(1));
t243 = qJD(3) * t223;
t249 = -t215 * t243 + t232 * t224;
t246 = qJD(1) * t223;
t245 = qJD(1) * t224;
t244 = qJD(3) * t217;
t242 = qJD(3) * t224;
t241 = qJD(5) * t215;
t237 = t248 * t215;
t234 = sin(pkin(11)) * pkin(4) + pkin(7) + qJ(2);
t231 = r_i_i_C(1) * t214 + r_i_i_C(2) * t216;
t230 = t233 * t223;
t228 = -t212 * t217 - cos(pkin(10)) * pkin(2) - pkin(1) - t237;
t227 = t215 * t242 + t232 * t223;
t225 = qJD(4) * t217 + t231 * t241 + (-t229 * t217 - t237) * qJD(3);
t211 = t214 * t230 - t249 * t216;
t210 = t249 * t214 + t216 * t230;
t209 = t214 * t252 + t227 * t216;
t208 = t227 * t214 - t216 * t252;
t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t224 * qJD(2) - t254 * t223 + (-t234 * t223 + t228 * t224) * qJD(1), t245, t225 * t224 - t226 * t246, -t215 * t246 + t217 * t242, t208 * r_i_i_C(1) + t209 * r_i_i_C(2), 0; -t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t223 * qJD(2) + t254 * t224 + (t228 * t223 + t234 * t224) * qJD(1), t246, t225 * t223 + t226 * t245, t215 * t245 + t217 * t243, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2), 0; 0, 0, t226 * qJD(3) - t231 * t240 + t239, qJD(3) * t215 (t214 * t241 - t216 * t244) * r_i_i_C(2) + (-t214 * t244 - t216 * t241) * r_i_i_C(1), 0;];
JaD_transl  = t1;
