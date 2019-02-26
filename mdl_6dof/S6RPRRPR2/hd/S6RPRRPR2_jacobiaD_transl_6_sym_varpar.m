% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR2
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
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:01:26
% EndTime: 2019-02-26 21:01:27
% DurationCPUTime: 0.27s
% Computational Cost: add. (549->63), mult. (484->90), div. (0->0), fcn. (382->12), ass. (0->52)
t250 = sin(qJ(3));
t247 = qJ(4) + pkin(11);
t235 = pkin(5) * cos(t247) + cos(qJ(4)) * pkin(4);
t233 = pkin(3) + t235;
t243 = qJ(6) + t247;
t237 = sin(t243);
t281 = r_i_i_C(2) * t237;
t238 = cos(t243);
t282 = r_i_i_C(1) * t238;
t258 = t233 - t281 + t282;
t252 = cos(qJ(3));
t278 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t285 = t278 * t252;
t254 = -t258 * t250 + t285;
t290 = qJD(1) * t254;
t234 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t247);
t230 = t234 * qJD(4);
t269 = t250 * qJD(5);
t289 = (-t233 * t250 + t285) * qJD(3) - t252 * t230 + t269;
t248 = qJ(1) + pkin(10);
t242 = cos(t248);
t246 = qJD(4) + qJD(6);
t264 = t246 * t252 - qJD(1);
t287 = t242 * t264;
t280 = r_i_i_C(2) * t238;
t262 = r_i_i_C(1) * t237 + t280;
t284 = t262 * t246 + t230;
t240 = sin(t248);
t272 = qJD(1) * t252;
t263 = -t246 + t272;
t271 = qJD(3) * t250;
t283 = -t240 * t271 + t263 * t242;
t279 = pkin(7) + t234;
t276 = t246 * t250;
t255 = t263 * t240 + t242 * t271;
t226 = t255 * t237 - t238 * t287;
t227 = t237 * t287 + t255 * t238;
t275 = t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
t260 = t264 * t240;
t228 = t283 * t237 + t238 * t260;
t229 = t237 * t260 - t283 * t238;
t274 = -t228 * r_i_i_C(1) + t229 * r_i_i_C(2);
t273 = qJD(1) * t250;
t270 = qJD(3) * t252;
t267 = t278 * t250;
t261 = t234 * t272 - t230;
t257 = -t233 * t252 - pkin(2) - t267;
t231 = t235 * qJD(4);
t256 = qJD(1) * t235 - t231 * t252 + t234 * t271;
t253 = qJD(5) * t252 + t284 * t250 + (-t258 * t252 - t267) * qJD(3);
t232 = t276 * t281;
t1 = [t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t242 * t231 - t289 * t240 + (-cos(qJ(1)) * pkin(1) - t279 * t240 + t257 * t242) * qJD(1), 0, -t240 * t290 + t253 * t242, t261 * t240 + t256 * t242 + t275, -t240 * t273 + t242 * t270, t275; -t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t240 * t231 + t289 * t242 + (-sin(qJ(1)) * pkin(1) + t279 * t242 + t257 * t240) * qJD(1), 0, t253 * t240 + t242 * t290, t256 * t240 - t261 * t242 + t274, t240 * t270 + t242 * t273, t274; 0, 0, t254 * qJD(3) - t284 * t252 + t269, t232 + (-t246 * t282 - t231) * t250 + (-t234 - t262) * t270, t271, -t270 * t280 + t232 + (-t237 * t270 - t238 * t276) * r_i_i_C(1);];
JaD_transl  = t1;
