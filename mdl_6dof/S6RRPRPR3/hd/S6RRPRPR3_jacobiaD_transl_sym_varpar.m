% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(10);
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
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:54
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
	t212 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:54
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (296->51), mult. (396->76), div. (0->0), fcn. (309->10), ass. (0->41)
	t230 = cos(qJ(4));
	t260 = t230 * pkin(4);
	t217 = pkin(3) + t260;
	t224 = qJ(2) + pkin(10);
	t220 = sin(t224);
	t222 = cos(t224);
	t258 = r_i_i_C(3) + qJ(5) + pkin(8);
	t241 = t258 * t222 - sin(qJ(2)) * pkin(2);
	t248 = qJD(4) * t222 - qJD(1);
	t251 = t220 * qJD(5);
	t227 = sin(qJ(4));
	t261 = pkin(4) * t227;
	t268 = (-t217 * t220 + t241) * qJD(2) - t248 * t261 - qJD(1) * (-qJ(3) - pkin(7)) + t251;
	t223 = qJ(4) + pkin(11);
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
	% StartTime: 2019-10-10 10:07:54
	% EndTime: 2019-10-10 10:07:54
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (544->63), mult. (506->88), div. (0->0), fcn. (399->12), ass. (0->53)
	t259 = qJ(4) + pkin(11);
	t245 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t259);
	t242 = t245 * qJD(4);
	t246 = pkin(5) * cos(t259) + cos(qJ(4)) * pkin(4);
	t244 = pkin(3) + t246;
	t260 = qJ(2) + pkin(10);
	t252 = sin(t260);
	t254 = cos(t260);
	t297 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t274 = t297 * t254 - sin(qJ(2)) * pkin(2);
	t284 = t252 * qJD(5);
	t310 = (-t244 * t252 + t274) * qJD(2) + (t245 + qJ(3) + pkin(7)) * qJD(1) - t254 * t242 + t284;
	t255 = qJ(6) + t259;
	t248 = sin(t255);
	t300 = r_i_i_C(2) * t248;
	t249 = cos(t255);
	t301 = r_i_i_C(1) * t249;
	t275 = t244 - t300 + t301;
	t269 = -t275 * t252 + t274;
	t267 = cos(qJ(1));
	t258 = qJD(4) + qJD(6);
	t281 = t254 * t258 - qJD(1);
	t308 = t267 * t281;
	t306 = -t297 * t252 - cos(qJ(2)) * pkin(2);
	t299 = r_i_i_C(2) * t249;
	t279 = r_i_i_C(1) * t248 + t299;
	t305 = t279 * t258 + t242;
	t291 = qJD(1) * t254;
	t280 = -t258 + t291;
	t264 = sin(qJ(1));
	t286 = qJD(2) * t264;
	t304 = -t252 * t286 + t280 * t267;
	t295 = t252 * t258;
	t285 = qJD(2) * t267;
	t271 = t252 * t285 + t280 * t264;
	t237 = t271 * t248 - t249 * t308;
	t238 = t248 * t308 + t271 * t249;
	t294 = t237 * r_i_i_C(1) + t238 * r_i_i_C(2);
	t277 = t281 * t264;
	t239 = t304 * t248 + t249 * t277;
	t240 = t248 * t277 - t304 * t249;
	t293 = -t239 * r_i_i_C(1) + t240 * r_i_i_C(2);
	t290 = qJD(1) * t264;
	t289 = qJD(1) * t267;
	t288 = qJD(2) * t252;
	t287 = qJD(2) * t254;
	t278 = t245 * t291 - t242;
	t243 = t246 * qJD(4);
	t272 = qJD(1) * t246 - t243 * t254 + t245 * t288;
	t270 = qJD(3) + t243 + (-t244 * t254 - pkin(1) + t306) * qJD(1);
	t268 = qJD(5) * t254 + t305 * t252 + (-t275 * t254 + t306) * qJD(2);
	t241 = t295 * t300;
	t1 = [t240 * r_i_i_C(1) + t239 * r_i_i_C(2) - t310 * t264 + t270 * t267, t268 * t267 - t269 * t290, t289, t278 * t264 + t272 * t267 + t294, -t252 * t290 + t254 * t285, t294; -t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t270 * t264 + t310 * t267, t268 * t264 + t269 * t289, t290, t272 * t264 - t278 * t267 + t293, t252 * t289 + t254 * t286, t293; 0, t269 * qJD(2) - t305 * t254 + t284, 0, t241 + (-t258 * t301 - t243) * t252 + (-t245 - t279) * t287, t288, -t287 * t299 + t241 + (-t248 * t287 - t249 * t295) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end