% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:42
% EndTime: 2019-02-26 20:43:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (259->43), mult. (300->68), div. (0->0), fcn. (232->10), ass. (0->35)
t223 = qJ(3) + pkin(10);
t219 = sin(t223);
t221 = cos(t223);
t252 = pkin(8) + r_i_i_C(3);
t235 = t252 * t221 - sin(qJ(3)) * pkin(3);
t258 = (-pkin(4) * t219 + t235) * qJD(3);
t226 = sin(qJ(5));
t228 = cos(qJ(5));
t237 = r_i_i_C(1) * t228 - r_i_i_C(2) * t226 + pkin(4);
t231 = -t237 * t219 + t235;
t239 = qJD(1) * t221 - qJD(5);
t256 = t228 * t239;
t254 = -t252 * t219 - cos(qJ(3)) * pkin(3);
t243 = qJD(5) * t221;
t240 = -qJD(1) + t243;
t246 = qJD(3) * t226;
t253 = -t219 * t246 + t240 * t228;
t224 = qJ(1) + pkin(9);
t220 = sin(t224);
t248 = qJD(1) * t220;
t222 = cos(t224);
t247 = qJD(1) * t222;
t245 = qJD(3) * t228;
t244 = qJD(5) * t219;
t238 = r_i_i_C(1) * t226 + r_i_i_C(2) * t228;
t236 = t239 * t226;
t233 = -pkin(4) * t221 - pkin(2) + t254;
t232 = t219 * t245 + t240 * t226;
t230 = t238 * t244 + (-t237 * t221 + t254) * qJD(3);
t225 = -qJ(4) - pkin(7);
t217 = t232 * t220 - t222 * t256;
t216 = t253 * t220 + t222 * t236;
t215 = t220 * t256 + t232 * t222;
t214 = t220 * t236 - t253 * t222;
t1 = [t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t222 * qJD(4) - t220 * t258 + (-cos(qJ(1)) * pkin(1) + t220 * t225 + t233 * t222) * qJD(1), 0, t230 * t222 - t231 * t248, t247, t214 * r_i_i_C(1) + t215 * r_i_i_C(2), 0; -t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t220 * qJD(4) + t222 * t258 + (-sin(qJ(1)) * pkin(1) - t222 * t225 + t233 * t220) * qJD(1), 0, t230 * t220 + t231 * t247, t248, -t216 * r_i_i_C(1) + t217 * r_i_i_C(2), 0; 0, 0, t231 * qJD(3) - t238 * t243, 0 (-t221 * t245 + t226 * t244) * r_i_i_C(2) + (-t221 * t246 - t228 * t244) * r_i_i_C(1), 0;];
JaD_transl  = t1;
