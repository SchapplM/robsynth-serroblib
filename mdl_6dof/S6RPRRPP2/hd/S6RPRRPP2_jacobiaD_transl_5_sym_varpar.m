% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:55
% EndTime: 2019-02-26 20:56:55
% DurationCPUTime: 0.23s
% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
t249 = sin(qJ(3));
t251 = cos(qJ(3));
t281 = pkin(8) + r_i_i_C(2);
t283 = t281 * t251;
t285 = (-pkin(3) * t249 + t283) * qJD(3);
t248 = sin(qJ(4));
t250 = cos(qJ(4));
t278 = r_i_i_C(3) + qJ(5);
t280 = -r_i_i_C(1) - pkin(4);
t257 = -t278 * t248 + t280 * t250;
t254 = -pkin(3) + t257;
t253 = t254 * t249 + t283;
t256 = t280 * t248 + t278 * t250;
t282 = t256 * qJD(4) + qJD(5) * t248;
t247 = qJ(1) + pkin(9);
t246 = cos(t247);
t277 = t246 * t248;
t276 = t248 * t251;
t275 = t250 * t251;
t245 = sin(t247);
t274 = qJD(1) * t245;
t273 = qJD(1) * t246;
t272 = qJD(3) * t249;
t271 = qJD(3) * t251;
t270 = qJD(4) * t248;
t269 = qJD(4) * t250;
t267 = t281 * t249;
t264 = t250 * t274;
t263 = t245 * t272;
t262 = t245 * t270;
t261 = t246 * t269;
t260 = t245 * t248 + t246 * t275;
t259 = t245 * t276 + t246 * t250;
t258 = -pkin(3) * t251 - pkin(2) - t267;
t255 = t245 * t269 + t248 * t273;
t252 = -t282 * t249 + (t254 * t251 - t267) * qJD(3);
t234 = t260 * qJD(1) - t250 * t263 - t251 * t262 - t261;
t233 = -t246 * t270 - t248 * t263 + t255 * t251 - t264;
t232 = t251 * t264 + (t250 * t272 + t251 * t270) * t246 - t255;
t231 = t259 * qJD(1) - t251 * t261 + t272 * t277 - t262;
t1 = [-t259 * qJD(5) + t280 * t234 - t278 * t233 - t245 * t285 + (-cos(qJ(1)) * pkin(1) - t245 * pkin(7) + t258 * t246) * qJD(1), 0, t252 * t246 - t253 * t274, t260 * qJD(5) - t280 * t231 - t278 * t232, -t231, 0; -(t245 * t250 - t246 * t276) * qJD(5) + t280 * t232 - t278 * t231 + t246 * t285 + (-sin(qJ(1)) * pkin(1) + t246 * pkin(7) + t258 * t245) * qJD(1), 0, t252 * t245 + t253 * t273 -(-t245 * t275 + t277) * qJD(5) + t278 * t234 + t280 * t233, t233, 0; 0, 0, t253 * qJD(3) + t282 * t251, t256 * t271 + (t257 * qJD(4) + t250 * qJD(5)) * t249, t248 * t271 + t249 * t269, 0;];
JaD_transl  = t1;
