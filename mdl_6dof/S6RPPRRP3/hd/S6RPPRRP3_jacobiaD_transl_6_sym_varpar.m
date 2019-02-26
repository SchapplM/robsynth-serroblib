% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:19
% EndTime: 2019-02-26 20:31:19
% DurationCPUTime: 0.28s
% Computational Cost: add. (328->55), mult. (522->87), div. (0->0), fcn. (435->8), ass. (0->39)
t237 = sin(qJ(4));
t239 = cos(qJ(4));
t236 = sin(qJ(5));
t238 = cos(qJ(5));
t256 = qJD(6) * t236;
t263 = r_i_i_C(3) + qJ(6);
t264 = r_i_i_C(1) + pkin(5);
t242 = -t256 + (t264 * t236 - t263 * t238) * qJD(5);
t243 = t263 * t236 + t264 * t238 + pkin(4);
t265 = pkin(8) + r_i_i_C(2);
t268 = t265 * t239;
t273 = t242 * t239 + (t243 * t237 - t268) * qJD(4);
t253 = t265 * t237;
t271 = t243 * t239 + t253;
t269 = -pkin(4) * t237 - qJ(3) + t268;
t266 = -pkin(2) - pkin(7);
t235 = qJ(1) + pkin(9);
t234 = cos(t235);
t262 = t234 * t236;
t261 = t237 * t238;
t260 = qJD(1) * t234;
t259 = qJD(4) * t237;
t258 = qJD(4) * t239;
t257 = qJD(5) * t239;
t255 = t238 * qJD(6);
t251 = t234 * t261;
t250 = t234 * t258;
t248 = qJD(5) * t237 + qJD(1);
t247 = qJD(1) * t237 + qJD(5);
t246 = t247 * t236;
t233 = sin(t235);
t245 = t233 * t261 + t262;
t241 = t237 * t256 + qJD(3) + (pkin(4) * t239 + t253) * qJD(4);
t240 = qJD(1) * t271;
t228 = t247 * t238 * t234 + (-t248 * t236 + t238 * t258) * t233;
t227 = t234 * t246 + (t236 * t258 + t248 * t238) * t233;
t226 = -t238 * t250 + (t233 * t238 + t237 * t262) * qJD(5) + t245 * qJD(1);
t225 = -qJD(5) * t251 + t233 * t246 - t236 * t250 - t238 * t260;
t1 = [t233 * t255 - t264 * t226 - t263 * t225 + t241 * t234 + (-cos(qJ(1)) * pkin(1) + t266 * t234 + t269 * t233) * qJD(1), 0, t260, -t273 * t233 + t234 * t240, t245 * qJD(6) - t264 * t227 + t263 * t228, t227; -t234 * t255 + t264 * t228 + t263 * t227 + t241 * t233 + (-sin(qJ(1)) * pkin(1) + t266 * t233 - t269 * t234) * qJD(1), 0, qJD(1) * t233, t233 * t240 + t273 * t234 -(-t233 * t236 + t251) * qJD(6) + t263 * t226 - t264 * t225, t225; 0, 0, 0, -t271 * qJD(4) + t242 * t237 (-t263 * t257 + t264 * t259) * t236 + (-t263 * t259 + (-t264 * qJD(5) + qJD(6)) * t239) * t238, -t236 * t259 + t238 * t257;];
JaD_transl  = t1;
