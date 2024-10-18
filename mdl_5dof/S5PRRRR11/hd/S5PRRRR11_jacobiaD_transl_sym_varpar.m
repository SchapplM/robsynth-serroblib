% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR11
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	% Symbolic code from jacobiaD_transl_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(10) + qJ(2);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(2), 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (73->19), mult. (110->37), div. (0->0), fcn. (94->6), ass. (0->19)
	t141 = sin(pkin(5));
	t152 = (pkin(7) + r_i_i_C(3)) * t141;
	t142 = cos(pkin(5));
	t143 = sin(qJ(3));
	t150 = t142 * t143;
	t144 = cos(qJ(3));
	t149 = t142 * t144;
	t140 = pkin(10) + qJ(2);
	t138 = sin(t140);
	t139 = cos(t140);
	t148 = t138 * t143 - t139 * t149;
	t147 = t138 * t144 + t139 * t150;
	t146 = t138 * t149 + t139 * t143;
	t145 = t138 * t150 - t139 * t144;
	t137 = t145 * qJD(2) + t148 * qJD(3);
	t136 = t146 * qJD(2) + t147 * qJD(3);
	t135 = t147 * qJD(2) + t146 * qJD(3);
	t134 = t148 * qJD(2) + t145 * qJD(3);
	t1 = [0, t137 * r_i_i_C(1) + t136 * r_i_i_C(2) + (-pkin(2) * t139 - t138 * t152) * qJD(2), t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0; 0, -t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(2) * t138 + t139 * t152) * qJD(2), -t136 * r_i_i_C(1) + t137 * r_i_i_C(2), 0, 0; 0, 0, (-r_i_i_C(1) * t143 - r_i_i_C(2) * t144) * t141 * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (343->48), mult. (255->76), div. (0->0), fcn. (176->12), ass. (0->41)
	t228 = qJ(3) + qJ(4);
	t222 = pkin(5) + t228;
	t251 = cos(t222) / 0.2e1;
	t227 = qJD(3) + qJD(4);
	t250 = -t227 / 0.2e1;
	t223 = pkin(5) - t228;
	t249 = cos(t223);
	t248 = sin(t222);
	t247 = pkin(3) * qJD(3);
	t224 = sin(t228);
	t246 = t224 * t227;
	t225 = cos(t228);
	t245 = t225 * t227;
	t217 = sin(t223);
	t244 = t227 * t217;
	t230 = cos(pkin(5));
	t231 = sin(qJ(3));
	t243 = t230 * t231;
	t232 = cos(qJ(3));
	t242 = t230 * t232;
	t212 = t248 * t250;
	t207 = t244 / 0.2e1 + t212;
	t209 = t249 / 0.2e1 + t251;
	t226 = pkin(10) + qJ(2);
	t220 = sin(t226);
	t221 = cos(t226);
	t202 = -t221 * t245 - t220 * t207 + (-t209 * t221 + t220 * t224) * qJD(2);
	t208 = t248 / 0.2e1 - t217 / 0.2e1;
	t206 = t209 * t227;
	t235 = -t220 * t206 - t221 * t246;
	t241 = t202 * r_i_i_C(1) + ((t208 * t221 + t220 * t225) * qJD(2) - t235) * r_i_i_C(2);
	t203 = t220 * t245 - t221 * t207 + (t209 * t220 + t221 * t224) * qJD(2);
	t237 = -t221 * t206 + t220 * t246;
	t240 = -t203 * r_i_i_C(1) + ((t208 * t220 - t221 * t225) * qJD(2) + t237) * r_i_i_C(2);
	t239 = (t227 * t251 + t249 * t250) * r_i_i_C(1) + (t212 - t244 / 0.2e1) * r_i_i_C(2);
	t238 = -t232 * pkin(3) - t225 * r_i_i_C(1) - pkin(2);
	t229 = sin(pkin(5));
	t236 = pkin(3) * t243 + t208 * r_i_i_C(1) + (-r_i_i_C(3) - pkin(7) - pkin(8)) * t229;
	t234 = t220 * t231 - t221 * t242;
	t233 = -t220 * t242 - t221 * t231;
	t1 = [0, t237 * r_i_i_C(1) + t203 * r_i_i_C(2) + t234 * t247 + (t236 * t220 + t238 * t221) * qJD(2), ((t220 * t243 - t221 * t232) * qJD(3) + t234 * qJD(2)) * pkin(3) + t241, t241, 0; 0, t235 * r_i_i_C(1) + t202 * r_i_i_C(2) + t233 * t247 + (t238 * t220 - t236 * t221) * qJD(2), ((-t220 * t232 - t221 * t243) * qJD(3) + t233 * qJD(2)) * pkin(3) + t240, t240, 0; 0, 0, -t229 * t231 * t247 + t239, t239, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:37
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (700->76), mult. (492->116), div. (0->0), fcn. (354->16), ass. (0->63)
	t251 = qJ(3) + qJ(4);
	t248 = qJ(5) + t251;
	t238 = pkin(5) + t248;
	t284 = cos(t238) / 0.2e1;
	t250 = qJD(3) + qJD(4);
	t245 = qJD(5) + t250;
	t283 = -t245 / 0.2e1;
	t239 = pkin(5) - t248;
	t282 = cos(t239);
	t281 = sin(t238);
	t246 = sin(t251);
	t280 = pkin(4) * t246;
	t249 = pkin(10) + qJ(2);
	t243 = sin(t249);
	t279 = t243 * t245;
	t244 = cos(t249);
	t278 = t244 * t245;
	t236 = sin(t239);
	t277 = t245 * t236;
	t247 = cos(t251);
	t276 = t247 * t250;
	t254 = sin(qJ(4));
	t255 = sin(qJ(3));
	t275 = t254 * t255;
	t257 = cos(qJ(3));
	t274 = t254 * t257;
	t229 = t281 * t283;
	t221 = t277 / 0.2e1 + t229;
	t225 = t282 / 0.2e1 + t284;
	t240 = sin(t248);
	t241 = cos(t248);
	t212 = -t241 * t278 - t243 * t221 + (-t225 * t244 + t240 * t243) * qJD(2);
	t224 = t281 / 0.2e1 - t236 / 0.2e1;
	t220 = t225 * t245;
	t261 = -t243 * t220 - t240 * t278;
	t273 = t212 * r_i_i_C(1) + ((t224 * t244 + t241 * t243) * qJD(2) - t261) * r_i_i_C(2);
	t213 = t241 * t279 - t244 * t221 + (t225 * t243 + t240 * t244) * qJD(2);
	t265 = -t244 * t220 + t240 * t279;
	t272 = -t213 * r_i_i_C(1) + ((t224 * t243 - t241 * t244) * qJD(2) + t265) * r_i_i_C(2);
	t271 = (t245 * t284 + t282 * t283) * r_i_i_C(1) + (t229 - t277 / 0.2e1) * r_i_i_C(2);
	t253 = cos(pkin(5));
	t256 = cos(qJ(4));
	t263 = t256 * t257 - t275;
	t270 = qJD(2) * t263 * t253 * pkin(4);
	t269 = qJD(2) * t246;
	t268 = qJD(3) * t255;
	t267 = qJD(3) * t257;
	t266 = -t257 * pkin(3) - pkin(4) * t247 - t241 * r_i_i_C(1) - pkin(2);
	t264 = -t255 * t256 - t274;
	t242 = t256 * pkin(4) + pkin(3);
	t252 = sin(pkin(5));
	t262 = t224 * r_i_i_C(1) + (pkin(4) * t274 + t255 * t242) * t253 + (-r_i_i_C(3) - pkin(8) - pkin(9) - pkin(7)) * t252;
	t260 = t264 * qJD(4);
	t259 = pkin(4) * (t264 * qJD(3) + t260);
	t258 = -t242 * t268 + (-t254 * t267 + t260) * pkin(4);
	t233 = -t255 * pkin(3) - t280;
	t228 = -pkin(3) * t267 - pkin(4) * t276;
	t227 = -pkin(3) * t268 - t250 * t280;
	t222 = (-pkin(4) * t275 + t242 * t257) * t253;
	t216 = t253 * t259;
	t215 = (t242 * t267 + (t263 * qJD(4) - t254 * t268) * pkin(4)) * t253;
	t214 = t258 * t253;
	t1 = [0, t265 * r_i_i_C(1) + t213 * r_i_i_C(2) - t243 * t227 - t244 * t215 + (t262 * t243 + t266 * t244) * qJD(2), -t243 * t214 + t244 * t228 + (-t222 * t244 - t233 * t243) * qJD(2) + t273, -t244 * t270 - t243 * t216 + (t243 * t269 - t244 * t276) * pkin(4) + t273, t273; 0, t261 * r_i_i_C(1) + t212 * r_i_i_C(2) + t244 * t227 - t243 * t215 + (t266 * t243 - t262 * t244) * qJD(2), t244 * t214 + t243 * t228 + (-t222 * t243 + t233 * t244) * qJD(2) + t272, -t243 * t270 + t244 * t216 + (-t243 * t276 - t244 * t269) * pkin(4) + t272, t272; 0, 0, t258 * t252 + t271, t252 * t259 + t271, t271;];
	JaD_transl = t1;
end