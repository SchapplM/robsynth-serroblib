% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:35
	% EndTime: 2019-12-29 18:50:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(8);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(6);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0; 0, -t33, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:32
	% EndTime: 2019-12-29 18:50:32
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(8);
	t210 = sin(t212);
	t211 = cos(t212);
	t243 = pkin(7) + r_i_i_C(3);
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
	t213 = -qJ(3) - pkin(6);
	t208 = -t229 * t239 + (t210 * t235 + t226) * t216;
	t207 = t230 * t217 * t216 + (-t210 * t236 + t229 * t219) * t214;
	t206 = t222 * t217 + t219 * t226;
	t205 = t222 * t214 - t230 * t239;
	t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t219 * qJD(3) - t246 * t236 + (t213 * t216 + t223 * t219) * qJD(1), t220 * t219 - t221 * t238, t237, t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0; -t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t216 * qJD(3) + t246 * t234 + (-t213 * t219 + t223 * t216) * qJD(1), t220 * t216 + t221 * t237, t238, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2), 0; 0, t221 * qJD(2) - t228 * t232, 0, (-t211 * t235 + t214 * t233) * r_i_i_C(2) + (-qJD(2) * t211 * t214 - t217 * t233) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:33
	% EndTime: 2019-12-29 18:50:33
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (308->56), mult. (536->83), div. (0->0), fcn. (446->8), ass. (0->45)
	t262 = qJ(2) + pkin(8);
	t260 = sin(t262);
	t261 = cos(t262);
	t303 = pkin(7) + r_i_i_C(2);
	t281 = t303 * t261 - sin(qJ(2)) * pkin(2);
	t309 = -pkin(3) * t260 + t281;
	t264 = sin(qJ(4));
	t267 = cos(qJ(4));
	t299 = r_i_i_C(3) + qJ(5);
	t304 = pkin(4) + r_i_i_C(1);
	t276 = -t299 * t264 - t304 * t267;
	t272 = -pkin(3) + t276;
	t271 = t272 * t260 + t281;
	t305 = t304 * t264 - t299 * t267;
	t308 = t305 * qJD(4) - qJD(5) * t264;
	t306 = -t303 * t260 - cos(qJ(2)) * pkin(2);
	t269 = cos(qJ(1));
	t298 = t264 * t269;
	t266 = sin(qJ(1));
	t297 = t266 * t264;
	t296 = t266 * t267;
	t295 = t267 * t269;
	t294 = qJD(1) * t266;
	t293 = qJD(1) * t269;
	t292 = qJD(2) * t261;
	t291 = qJD(2) * t266;
	t290 = qJD(2) * t269;
	t289 = qJD(4) * t267;
	t288 = qJD(4) * t269;
	t285 = t260 * t291;
	t284 = qJD(4) * t297;
	t283 = t260 * t290;
	t282 = t267 * t288;
	t279 = t261 * t295 + t297;
	t278 = t261 * t297 + t295;
	t277 = -pkin(3) * t261 - pkin(1) + t306;
	t274 = t264 * t288 + t267 * t294;
	t273 = t264 * t293 + t266 * t289;
	t270 = t308 * t260 + (t272 * t261 + t306) * qJD(2);
	t263 = -qJ(3) - pkin(6);
	t248 = t279 * qJD(1) - t261 * t284 - t267 * t285 - t282;
	t247 = t273 * t261 - t264 * t285 - t274;
	t246 = t274 * t261 + t267 * t283 - t273;
	t245 = t278 * qJD(1) - t261 * t282 + t264 * t283 - t284;
	t1 = [-t278 * qJD(5) + t269 * qJD(3) - t304 * t248 - t299 * t247 - t309 * t291 + (t266 * t263 + t277 * t269) * qJD(1), t270 * t269 - t271 * t294, t293, t279 * qJD(5) + t304 * t245 - t299 * t246, -t245; -(-t261 * t298 + t296) * qJD(5) + t266 * qJD(3) - t304 * t246 - t299 * t245 + t309 * t290 + (-t269 * t263 + t277 * t266) * qJD(1), t270 * t266 + t271 * t293, t294, -(-t261 * t296 + t298) * qJD(5) + t299 * t248 - t304 * t247, t247; 0, t271 * qJD(2) - t308 * t261, 0, -t305 * t292 + (t276 * qJD(4) + t267 * qJD(5)) * t260, t260 * t289 + t264 * t292;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end