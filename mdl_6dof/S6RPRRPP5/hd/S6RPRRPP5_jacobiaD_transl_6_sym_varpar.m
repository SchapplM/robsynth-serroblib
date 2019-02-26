% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP5
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
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:39
% EndTime: 2019-02-26 20:58:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (389->60), mult. (652->88), div. (0->0), fcn. (546->7), ass. (0->46)
t221 = pkin(9) + qJ(3);
t219 = sin(t221);
t220 = cos(t221);
t248 = t219 * qJD(6);
t223 = sin(qJ(4));
t249 = qJD(5) * t223;
t245 = pkin(8) - r_i_i_C(3) - qJ(6);
t262 = t245 * t220;
t265 = (-pkin(3) * t219 + t262) * qJD(3) + t220 * t249 - t248;
t225 = cos(qJ(4));
t246 = r_i_i_C(1) + pkin(5) + pkin(4);
t259 = r_i_i_C(2) + qJ(5);
t261 = t246 * t223 - t259 * t225;
t264 = t261 * qJD(4) - t249;
t231 = -t259 * t223 - t246 * t225;
t229 = -pkin(3) + t231;
t228 = t229 * t219 + t262;
t224 = sin(qJ(1));
t258 = t224 * t223;
t226 = cos(qJ(1));
t257 = t226 * t225;
t256 = qJD(1) * t224;
t255 = qJD(1) * t226;
t254 = qJD(3) * t220;
t253 = qJD(3) * t224;
t252 = qJD(3) * t226;
t251 = qJD(4) * t225;
t250 = qJD(4) * t226;
t247 = t225 * qJD(5);
t244 = t219 * t253;
t243 = qJD(4) * t258;
t242 = t219 * t252;
t241 = t225 * t250;
t240 = qJD(2) - t247;
t239 = t245 * t219;
t236 = t220 * t257 + t258;
t234 = t223 * t250 + t225 * t256;
t233 = t223 * t255 + t224 * t251;
t232 = -pkin(3) * t220 - cos(pkin(9)) * pkin(2) - pkin(1) - t239;
t227 = -qJD(6) * t220 + t264 * t219 + (t229 * t220 - t239) * qJD(3);
t222 = -pkin(7) - qJ(2);
t207 = t236 * qJD(1) - t220 * t243 - t225 * t244 - t241;
t206 = t233 * t220 - t223 * t244 - t234;
t205 = t234 * t220 + t225 * t242 - t233;
t204 = t223 * t242 - t220 * t241 - t243 + (t220 * t258 + t257) * qJD(1);
t1 = [t240 * t226 - t259 * t206 - t246 * t207 - t265 * t224 + (t224 * t222 + t232 * t226) * qJD(1), t255, t227 * t226 - t228 * t256, t236 * qJD(5) + t246 * t204 - t259 * t205, -t204, t219 * t256 - t220 * t252; t240 * t224 - t259 * t204 - t246 * t205 + t265 * t226 + (-t226 * t222 + t232 * t224) * qJD(1), t256, t227 * t224 + t228 * t255 -(-t224 * t220 * t225 + t226 * t223) * qJD(5) + t259 * t207 - t246 * t206, t206, -t219 * t255 - t220 * t253; 0, 0, t228 * qJD(3) - t264 * t220 - t248, -t261 * t254 + (t231 * qJD(4) + t247) * t219, t219 * t251 + t223 * t254, -qJD(3) * t219;];
JaD_transl  = t1;
