% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:22
% EndTime: 2019-02-26 21:00:22
% DurationCPUTime: 0.23s
% Computational Cost: add. (228->59), mult. (692->91), div. (0->0), fcn. (588->6), ass. (0->41)
t241 = sin(qJ(4));
t244 = cos(qJ(4));
t263 = pkin(4) + r_i_i_C(3) + qJ(6);
t274 = r_i_i_C(2) + qJ(5);
t252 = t263 * t241 - t274 * t244;
t264 = pkin(5) + pkin(8) + r_i_i_C(1);
t280 = t264 * qJD(3) - t252 * qJD(4) + qJD(5) * t241 + qJD(6) * t244;
t253 = -t274 * t241 - t263 * t244;
t279 = -pkin(3) + t253;
t242 = sin(qJ(3));
t245 = cos(qJ(3));
t277 = -pkin(3) * t242 + t264 * t245 - qJ(2);
t275 = -pkin(1) - pkin(7);
t243 = sin(qJ(1));
t273 = t243 * t241;
t272 = t243 * t244;
t246 = cos(qJ(1));
t271 = t246 * t241;
t270 = t246 * t244;
t269 = qJD(1) * t243;
t268 = qJD(1) * t246;
t267 = qJD(3) * t242;
t266 = qJD(3) * t246;
t265 = qJD(4) * t245;
t262 = t242 * t273;
t261 = t242 * t270;
t260 = t245 * t266;
t258 = t264 * t242;
t257 = qJD(1) * t242 + qJD(4);
t255 = t242 * t272 + t271;
t254 = t242 * t271 + t272;
t250 = qJD(3) * t243 * t245 + t257 * t246;
t249 = qJD(2) + (pkin(3) * t245 + t258) * qJD(3);
t248 = qJD(3) * t279;
t234 = t261 - t273;
t231 = t262 - t270;
t230 = -qJD(4) * t262 - t241 * t269 + t250 * t244;
t229 = (qJD(4) * t242 + qJD(1)) * t272 + t250 * t241;
t228 = t255 * qJD(1) + t254 * qJD(4) - t244 * t260;
t227 = -qJD(4) * t261 - t241 * t260 - t244 * t268 + t257 * t273;
t1 = [t254 * qJD(5) + t234 * qJD(6) - t274 * t227 - t263 * t228 + t249 * t246 + (t277 * t243 + t275 * t246) * qJD(1), t268 (-t245 * t279 + t258) * t268 + (t242 * t248 + t280 * t245) * t243, qJD(5) * t255 - t231 * qJD(6) - t263 * t229 + t274 * t230, t229, t230; t231 * qJD(5) + t255 * qJD(6) + t274 * t229 + t263 * t230 + t249 * t243 + (t275 * t243 - t277 * t246) * qJD(1), t269 (t264 * t269 - t266 * t279) * t242 + (-t246 * t280 - t269 * t279) * t245, -t234 * qJD(5) + qJD(6) * t254 - t263 * t227 + t274 * t228, t227, t228; 0, 0, -t242 * t280 + t245 * t248, t252 * t267 + (t253 * qJD(4) + qJD(5) * t244 - qJD(6) * t241) * t245, -t241 * t267 + t244 * t265, -t241 * t265 - t244 * t267;];
JaD_transl  = t1;
