% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPP1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(9);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(7);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t30 * t40 + t32 * t37) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t30 * t37 + t32 * t40) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:41
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(9);
	t210 = sin(t212);
	t211 = cos(t212);
	t243 = pkin(8) + r_i_i_C(3);
	t225 = t243 * t211 - sin(qJ(2)) * pkin(2);
	t246 = -pkin(3) * t210 + t225;
	t214 = sin(qJ(4));
	t217 = cos(qJ(4));
	t227 = r_i_i_C(1) * t217 - r_i_i_C(2) * t214 + pkin(3);
	t221 = -t227 * t210 + t225;
	t244 = -t243 * t210 - cos(qJ(2)) * pkin(2);
	t219 = cos(qJ(1));
	t239 = t217 * t219;
	t216 = sin(qJ(1));
	t238 = qJD(1) * t216;
	t237 = qJD(1) * t219;
	t236 = qJD(2) * t216;
	t235 = qJD(2) * t217;
	t234 = qJD(2) * t219;
	t233 = qJD(4) * t210;
	t232 = qJD(4) * t211;
	t230 = -qJD(1) + t232;
	t229 = qJD(1) * t211 - qJD(4);
	t228 = r_i_i_C(1) * t214 + r_i_i_C(2) * t217;
	t226 = t230 * t214;
	t223 = -pkin(3) * t211 - pkin(1) + t244;
	t222 = t210 * t234 + t229 * t216;
	t220 = t228 * t233 + (-t227 * t211 + t244) * qJD(2);
	t213 = -qJ(3) - pkin(7);
	t208 = -t229 * t239 + (t210 * t235 + t226) * t216;
	t207 = t230 * t217 * t216 + (-t210 * t236 + t229 * t219) * t214;
	t206 = t222 * t217 + t219 * t226;
	t205 = t222 * t214 - t230 * t239;
	t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t219 * qJD(3) - t246 * t236 + (t213 * t216 + t223 * t219) * qJD(1), t220 * t219 - t221 * t238, t237, t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0, 0; -t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t216 * qJD(3) + t246 * t234 + (-t213 * t219 + t223 * t216) * qJD(1), t220 * t216 + t221 * t237, t238, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2), 0, 0; 0, t221 * qJD(2) - t228 * t232, 0, (-t211 * t235 + t214 * t233) * r_i_i_C(2) + (-qJD(2) * t211 * t214 - t217 * t233) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:41
	% EndTime: 2019-10-10 09:55:41
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (296->51), mult. (396->76), div. (0->0), fcn. (309->10), ass. (0->41)
	t230 = cos(qJ(4));
	t260 = t230 * pkin(4);
	t217 = pkin(3) + t260;
	t224 = qJ(2) + pkin(9);
	t220 = sin(t224);
	t222 = cos(t224);
	t258 = r_i_i_C(3) + qJ(5) + pkin(8);
	t241 = t258 * t222 - sin(qJ(2)) * pkin(2);
	t248 = qJD(4) * t222 - qJD(1);
	t251 = t220 * qJD(5);
	t227 = sin(qJ(4));
	t261 = pkin(4) * t227;
	t268 = (-t217 * t220 + t241) * qJD(2) - t248 * t261 - qJD(1) * (-qJ(3) - pkin(7)) + t251;
	t223 = qJ(4) + pkin(10);
	t219 = sin(t223);
	t221 = cos(t223);
	t246 = r_i_i_C(1) * t221 - r_i_i_C(2) * t219;
	t242 = t217 + t246;
	t235 = -t242 * t220 + t241;
	t229 = sin(qJ(1));
	t244 = t248 * t229;
	t232 = cos(qJ(1));
	t245 = t248 * t232;
	t265 = -t258 * t220 - cos(qJ(2)) * pkin(2);
	t247 = qJD(1) * t222 - qJD(4);
	t254 = qJD(2) * t229;
	t264 = -t220 * t254 + t247 * t232;
	t256 = qJD(1) * t229;
	t255 = qJD(1) * t232;
	t253 = qJD(2) * t232;
	t252 = qJD(4) * t220;
	t239 = r_i_i_C(1) * t219 + r_i_i_C(2) * t221 + t261;
	t238 = t239 * t222;
	t236 = t220 * t253 + t247 * t229;
	t234 = qJD(4) * t260 + qJD(3) + (-t217 * t222 - pkin(1) + t265) * qJD(1);
	t233 = qJD(5) * t222 + t239 * t252 + (-t242 * t222 + t265) * qJD(2);
	t216 = t219 * t244 - t221 * t264;
	t215 = t264 * t219 + t221 * t244;
	t214 = t219 * t245 + t236 * t221;
	t213 = t236 * t219 - t221 * t245;
	t1 = [t216 * r_i_i_C(1) + t215 * r_i_i_C(2) - t268 * t229 + t234 * t232, t233 * t232 - t235 * t256, t255, t213 * r_i_i_C(1) + t214 * r_i_i_C(2) + (t236 * t227 - t230 * t245) * pkin(4), -t220 * t256 + t222 * t253, 0; -t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t234 * t229 + t268 * t232, t233 * t229 + t235 * t255, t256, -t215 * r_i_i_C(1) + t216 * r_i_i_C(2) + (-t227 * t264 - t230 * t244) * pkin(4), t220 * t255 + t222 * t254, 0; 0, t235 * qJD(2) - qJD(4) * t238 + t251, 0, (-t246 - t260) * t252 - qJD(2) * t238, qJD(2) * t220, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:41
	% EndTime: 2019-10-10 09:55:42
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (521->73), mult. (636->102), div. (0->0), fcn. (525->10), ass. (0->53)
	t282 = cos(qJ(4));
	t320 = pkin(4) * t282;
	t269 = pkin(3) + t320;
	t276 = qJ(2) + pkin(9);
	t272 = sin(t276);
	t274 = cos(t276);
	t319 = r_i_i_C(2) + qJ(5) + pkin(8);
	t295 = t319 * t274 - sin(qJ(2)) * pkin(2);
	t298 = qJD(4) * t274 - qJD(1);
	t306 = t272 * qJD(5);
	t275 = qJ(4) + pkin(10);
	t271 = sin(t275);
	t307 = qJD(6) * t271;
	t279 = sin(qJ(4));
	t321 = pkin(4) * t279;
	t330 = (-t269 * t272 + t295) * qJD(2) - t298 * t321 - qJD(1) * (-qJ(3) - pkin(7)) + t274 * t307 + t306;
	t273 = cos(t275);
	t318 = r_i_i_C(3) + qJ(6);
	t324 = -r_i_i_C(1) - pkin(5);
	t292 = -t318 * t271 + t324 * t273;
	t290 = -t269 + t292;
	t287 = t290 * t272 + t295;
	t327 = -t319 * t272 - cos(qJ(2)) * pkin(2);
	t288 = t324 * t271 + t318 * t273 - t321;
	t325 = t288 * qJD(4) + t307;
	t284 = cos(qJ(1));
	t316 = t273 * t284;
	t281 = sin(qJ(1));
	t315 = t281 * t271;
	t314 = qJD(1) * t281;
	t313 = qJD(1) * t284;
	t312 = qJD(2) * t274;
	t311 = qJD(2) * t281;
	t310 = qJD(2) * t284;
	t309 = qJD(4) * t281;
	t308 = qJD(4) * t284;
	t305 = t273 * qJD(6);
	t304 = t272 * t311;
	t303 = t272 * t310;
	t302 = t271 * t309;
	t301 = t271 * t308;
	t300 = t273 * t308;
	t297 = qJD(1) * t274 - qJD(4);
	t296 = t298 * t282;
	t293 = t274 * t316 + t315;
	t291 = t271 * t313 + t273 * t309;
	t286 = qJD(4) * t320 - t305 + qJD(3) + (-t269 * t274 - pkin(1) + t327) * qJD(1);
	t285 = qJD(5) * t274 - t325 * t272 + (t290 * t274 + t327) * qJD(2);
	t258 = t293 * qJD(1) - t273 * t304 - t274 * t302 - t300;
	t257 = -t271 * t304 - t273 * t314 + t291 * t274 - t301;
	t256 = t274 * t301 + (t274 * t314 + t303) * t273 - t291;
	t255 = t271 * t303 - t274 * t300 - t302 + (t274 * t315 + t316) * qJD(1);
	t1 = [-t318 * t257 + t324 * t258 - t330 * t281 + t286 * t284, t285 * t284 - t287 * t314, t313, t293 * qJD(6) - t318 * t256 - t324 * t255 + (-t284 * t296 + (t297 * t281 + t303) * t279) * pkin(4), -t272 * t314 + t274 * t310, -t255; -t318 * t255 + t324 * t256 + t286 * t281 + t330 * t284, t285 * t281 + t287 * t313, t314, -(-t281 * t274 * t273 + t284 * t271) * qJD(6) + t318 * t258 + t324 * t257 + (-t281 * t296 + (-t297 * t284 + t304) * t279) * pkin(4), t272 * t313 + t274 * t311, t257; 0, t287 * qJD(2) + t325 * t274 + t306, 0, t288 * t312 + (t305 + (t292 - t320) * qJD(4)) * t272, qJD(2) * t272, qJD(4) * t272 * t273 + t271 * t312;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end