% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:48
% EndTime: 2019-02-26 21:03:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (289->50), mult. (376->75), div. (0->0), fcn. (296->9), ass. (0->41)
t229 = cos(qJ(4));
t257 = t229 * pkin(4);
t218 = pkin(3) + t257;
t223 = pkin(10) + qJ(3);
t219 = sin(t223);
t221 = cos(t223);
t244 = qJD(4) * t221 - qJD(1);
t249 = t219 * qJD(5);
t227 = sin(qJ(4));
t258 = pkin(4) * t227;
t256 = r_i_i_C(3) + qJ(5) + pkin(8);
t261 = t256 * t221;
t264 = (-t218 * t219 + t261) * qJD(3) - t244 * t258 - qJD(1) * (-pkin(7) - qJ(2)) + t249;
t224 = qJ(4) + pkin(11);
t220 = sin(t224);
t222 = cos(t224);
t242 = r_i_i_C(1) * t222 - r_i_i_C(2) * t220;
t238 = t218 + t242;
t233 = -t238 * t219 + t261;
t228 = sin(qJ(1));
t240 = t244 * t228;
t230 = cos(qJ(1));
t241 = t244 * t230;
t243 = qJD(1) * t221 - qJD(4);
t252 = qJD(3) * t228;
t260 = -t219 * t252 + t243 * t230;
t254 = qJD(1) * t228;
t253 = qJD(1) * t230;
t251 = qJD(3) * t230;
t250 = qJD(4) * t219;
t247 = t256 * t219;
t237 = r_i_i_C(1) * t220 + r_i_i_C(2) * t222 + t258;
t236 = t237 * t221;
t234 = t219 * t251 + t243 * t228;
t232 = qJD(4) * t257 + qJD(2) + (-t218 * t221 - cos(pkin(10)) * pkin(2) - pkin(1) - t247) * qJD(1);
t231 = qJD(5) * t221 + t237 * t250 + (-t238 * t221 - t247) * qJD(3);
t216 = t220 * t240 - t222 * t260;
t215 = t260 * t220 + t222 * t240;
t214 = t220 * t241 + t234 * t222;
t213 = t234 * t220 - t222 * t241;
t1 = [t216 * r_i_i_C(1) + t215 * r_i_i_C(2) - t264 * t228 + t232 * t230, t253, t231 * t230 - t233 * t254, t213 * r_i_i_C(1) + t214 * r_i_i_C(2) + (t234 * t227 - t229 * t241) * pkin(4), -t219 * t254 + t221 * t251, 0; -t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t232 * t228 + t264 * t230, t254, t231 * t228 + t233 * t253, -t215 * r_i_i_C(1) + t216 * r_i_i_C(2) + (-t227 * t260 - t229 * t240) * pkin(4), t219 * t253 + t221 * t252, 0; 0, 0, t233 * qJD(3) - qJD(4) * t236 + t249 (-t242 - t257) * t250 - qJD(3) * t236, qJD(3) * t219, 0;];
JaD_transl  = t1;
