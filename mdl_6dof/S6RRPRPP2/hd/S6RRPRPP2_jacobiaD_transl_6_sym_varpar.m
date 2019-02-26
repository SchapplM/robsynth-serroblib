% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:22
% EndTime: 2019-02-26 21:35:23
% DurationCPUTime: 0.32s
% Computational Cost: add. (396->59), mult. (672->85), div. (0->0), fcn. (559->8), ass. (0->44)
t223 = qJ(2) + pkin(9);
t221 = sin(t223);
t222 = cos(t223);
t248 = pkin(8) - r_i_i_C(3) - qJ(6);
t241 = t248 * t222 - sin(qJ(2)) * pkin(2);
t251 = t221 * qJD(6);
t225 = sin(qJ(4));
t252 = qJD(5) * t225;
t271 = (-pkin(3) * t221 + t241) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t222 * t252 - t251;
t228 = cos(qJ(4));
t249 = r_i_i_C(1) + pkin(5) + pkin(4);
t262 = r_i_i_C(2) + qJ(5);
t237 = -t262 * t225 - t249 * t228;
t235 = -pkin(3) + t237;
t232 = t235 * t221 + t241;
t266 = t249 * t225 - t262 * t228;
t269 = t266 * qJD(4) - t252;
t267 = -t248 * t221 - cos(qJ(2)) * pkin(2);
t227 = sin(qJ(1));
t261 = t227 * t225;
t230 = cos(qJ(1));
t260 = t230 * t228;
t259 = qJD(1) * t227;
t258 = qJD(1) * t230;
t257 = qJD(2) * t222;
t256 = qJD(2) * t227;
t255 = qJD(2) * t230;
t254 = qJD(4) * t228;
t253 = qJD(4) * t230;
t250 = t228 * qJD(5);
t247 = t221 * t256;
t246 = qJD(4) * t261;
t245 = t221 * t255;
t244 = t228 * t253;
t242 = t222 * t260 + t261;
t239 = t225 * t253 + t228 * t259;
t238 = t225 * t258 + t227 * t254;
t233 = -t250 + qJD(3) + (-pkin(3) * t222 - pkin(1) + t267) * qJD(1);
t231 = -qJD(6) * t222 + t269 * t221 + (t235 * t222 + t267) * qJD(2);
t209 = t242 * qJD(1) - t222 * t246 - t228 * t247 - t244;
t208 = t238 * t222 - t225 * t247 - t239;
t207 = t239 * t222 + t228 * t245 - t238;
t206 = t225 * t245 - t222 * t244 - t246 + (t222 * t261 + t260) * qJD(1);
t1 = [-t262 * t208 - t249 * t209 - t271 * t227 + t233 * t230, t231 * t230 - t232 * t259, t258, t242 * qJD(5) + t249 * t206 - t262 * t207, -t206, t221 * t259 - t222 * t255; -t262 * t206 - t249 * t207 + t233 * t227 + t271 * t230, t231 * t227 + t232 * t258, t259 -(-t227 * t222 * t228 + t230 * t225) * qJD(5) + t262 * t209 - t249 * t208, t208, -t221 * t258 - t222 * t256; 0, t232 * qJD(2) - t269 * t222 - t251, 0, -t266 * t257 + (t237 * qJD(4) + t250) * t221, t221 * t254 + t225 * t257, -qJD(2) * t221;];
JaD_transl  = t1;
