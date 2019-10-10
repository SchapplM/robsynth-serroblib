% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
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
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 09:59:11
	% EndTime: 2019-10-10 09:59:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(9));
	t159 = cos(pkin(9));
	t169 = r_i_i_C(1) * t159 - r_i_i_C(2) * t158 + pkin(2);
	t174 = r_i_i_C(3) + qJ(3);
	t175 = t169 * t160 - t174 * t162;
	t176 = t175 * qJD(2) - t160 * qJD(3);
	t161 = sin(qJ(1));
	t173 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t172 = qJD(1) * t163;
	t171 = qJD(2) * t163;
	t168 = t158 * r_i_i_C(1) + t159 * r_i_i_C(2) + pkin(7);
	t167 = -t174 * t160 - t169 * t162;
	t165 = -pkin(1) + t167;
	t1 = [t176 * t161 + (-t168 * t161 + t165 * t163) * qJD(1), (t169 * t173 - t174 * t171) * t160 + (-t174 * t173 + (-t169 * qJD(2) + qJD(3)) * t163) * t162, -t160 * t173 + t162 * t171, 0, 0, 0; -t176 * t163 + (t165 * t161 + t168 * t163) * qJD(1), -t175 * t172 + (t167 * qJD(2) + qJD(3) * t162) * t161, t161 * qJD(2) * t162 + t160 * t172, 0, 0, 0; 0, -t176, qJD(2) * t160, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:11
	% EndTime: 2019-10-10 09:59:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (170->42), mult. (303->71), div. (0->0), fcn. (242->8), ass. (0->37)
	t203 = cos(pkin(9)) * pkin(3) + pkin(2);
	t209 = sin(qJ(2));
	t227 = t209 * qJD(3);
	t211 = cos(qJ(2));
	t236 = r_i_i_C(3) + pkin(8) + qJ(3);
	t238 = t236 * t211;
	t242 = (-t203 * t209 + t238) * qJD(2) + t227;
	t206 = pkin(9) + qJ(4);
	t204 = sin(t206);
	t205 = cos(t206);
	t217 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
	t214 = -t217 * t209 + t238;
	t212 = cos(qJ(1));
	t228 = qJD(4) * t211;
	t221 = -qJD(1) + t228;
	t240 = t212 * t221;
	t210 = sin(qJ(1));
	t220 = qJD(1) * t211 - qJD(4);
	t232 = qJD(2) * t209;
	t237 = -t210 * t232 + t220 * t212;
	t234 = qJD(1) * t210;
	t233 = qJD(1) * t212;
	t231 = qJD(2) * t211;
	t230 = qJD(2) * t212;
	t229 = qJD(4) * t209;
	t226 = pkin(3) * sin(pkin(9)) + pkin(7);
	t224 = t236 * t209;
	t219 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
	t218 = t221 * t210;
	t216 = -t203 * t211 - pkin(1) - t224;
	t215 = t209 * t230 + t220 * t210;
	t213 = qJD(3) * t211 + t219 * t229 + (-t217 * t211 - t224) * qJD(2);
	t202 = t204 * t218 - t237 * t205;
	t201 = t237 * t204 + t205 * t218;
	t200 = t204 * t240 + t215 * t205;
	t199 = t215 * t204 - t205 * t240;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t242 * t210 + (-t226 * t210 + t216 * t212) * qJD(1), t213 * t212 - t214 * t234, -t209 * t234 + t211 * t230, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t242 * t212 + (t216 * t210 + t226 * t212) * qJD(1), t213 * t210 + t214 * t233, t209 * t233 + t210 * t231, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, t214 * qJD(2) - t219 * t228 + t227, t232, (t204 * t229 - t205 * t231) * r_i_i_C(2) + (-t204 * t231 - t205 * t229) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:12
	% EndTime: 2019-10-10 09:59:12
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (335->61), mult. (543->94), div. (0->0), fcn. (458->8), ass. (0->45)
	t256 = sin(qJ(2));
	t250 = cos(pkin(9)) * pkin(3) + pkin(2);
	t253 = pkin(9) + qJ(4);
	t251 = sin(t253);
	t252 = cos(t253);
	t290 = r_i_i_C(3) + qJ(5);
	t292 = pkin(4) - r_i_i_C(2);
	t293 = t290 * t251 + t292 * t252 + t250;
	t258 = cos(qJ(2));
	t291 = r_i_i_C(1) + pkin(8) + qJ(3);
	t295 = t291 * t258;
	t299 = t293 * t256 - t295;
	t276 = t256 * qJD(3);
	t278 = qJD(5) * t251;
	t298 = (-t250 * t256 + t295) * qJD(2) + t258 * t278 + t276;
	t296 = -t278 + (t292 * t251 - t290 * t252) * qJD(4);
	t257 = sin(qJ(1));
	t288 = t257 * t258;
	t259 = cos(qJ(1));
	t287 = t259 * t252;
	t286 = qJD(1) * t257;
	t285 = qJD(1) * t259;
	t284 = qJD(2) * t256;
	t283 = qJD(2) * t258;
	t282 = qJD(2) * t259;
	t281 = qJD(4) * t256;
	t280 = qJD(4) * t257;
	t279 = qJD(4) * t259;
	t277 = t252 * qJD(5);
	t275 = pkin(3) * sin(pkin(9)) + pkin(7);
	t274 = t257 * t284;
	t273 = t251 * t280;
	t272 = t256 * t282;
	t271 = t251 * t279;
	t270 = t252 * t279;
	t269 = t291 * t256;
	t266 = t257 * t251 + t258 * t287;
	t264 = -t250 * t258 - pkin(1) - t269;
	t263 = t251 * t285 + t252 * t280;
	t260 = qJD(3) * t258 + t296 * t256 + (-t258 * t293 - t269) * qJD(2);
	t239 = t266 * qJD(1) - t252 * t274 - t258 * t273 - t270;
	t238 = -t251 * t274 - t252 * t286 + t263 * t258 - t271;
	t237 = t258 * t271 + (t258 * t286 + t272) * t252 - t263;
	t236 = t251 * t272 - t258 * t270 - t273 + (t251 * t288 + t287) * qJD(1);
	t1 = [-t259 * t277 - t292 * t239 - t290 * t238 - t298 * t257 + (-t275 * t257 + t264 * t259) * qJD(1), t260 * t259 + t299 * t286, -t256 * t286 + t258 * t282, t266 * qJD(5) + t292 * t236 - t290 * t237, -t236, 0; -t257 * t277 - t292 * t237 - t290 * t236 + t298 * t259 + (t264 * t257 + t275 * t259) * qJD(1), t260 * t257 - t285 * t299, t256 * t285 + t257 * t283, -(t259 * t251 - t252 * t288) * qJD(5) + t290 * t239 - t292 * t238, t238, 0; 0, -qJD(2) * t299 - t296 * t258 + t276, t284, (-t290 * t281 - t292 * t283) * t251 + (t290 * t283 + (-t292 * qJD(4) + qJD(5)) * t256) * t252, t251 * t283 + t252 * t281, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:12
	% EndTime: 2019-10-10 09:59:12
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (449->66), mult. (717->99), div. (0->0), fcn. (613->8), ass. (0->48)
	t251 = cos(pkin(9)) * pkin(3) + pkin(2);
	t257 = sin(qJ(2));
	t259 = cos(qJ(2));
	t278 = pkin(5) + r_i_i_C(1) + pkin(8) + qJ(3);
	t296 = t278 * t259;
	t301 = (-t251 * t257 + t296) * qJD(2) + t257 * qJD(3);
	t254 = pkin(9) + qJ(4);
	t252 = sin(t254);
	t253 = cos(t254);
	t279 = pkin(4) + r_i_i_C(3) + qJ(6);
	t293 = r_i_i_C(2) + qJ(5);
	t294 = t279 * t252 - t293 * t253;
	t300 = -t278 * qJD(2) + t294 * qJD(4) - qJD(5) * t252 - qJD(6) * t253;
	t265 = -t293 * t252 - t279 * t253;
	t263 = -t251 + t265;
	t298 = t263 * t257 + t296;
	t258 = sin(qJ(1));
	t291 = t258 * t259;
	t260 = cos(qJ(1));
	t290 = t260 * t252;
	t289 = t260 * t253;
	t288 = qJD(1) * t258;
	t287 = qJD(1) * t260;
	t286 = qJD(2) * t257;
	t285 = qJD(2) * t259;
	t284 = qJD(2) * t260;
	t283 = qJD(4) * t257;
	t282 = qJD(4) * t258;
	t281 = qJD(4) * t260;
	t277 = pkin(3) * sin(pkin(9)) + pkin(7);
	t276 = t258 * t286;
	t275 = t257 * t284;
	t274 = t252 * t282;
	t273 = t252 * t281;
	t272 = t253 * t281;
	t268 = t258 * t252 + t259 * t289;
	t237 = t252 * t291 + t289;
	t267 = t252 * t287 + t253 * t282;
	t266 = -t251 * t259 - t278 * t257 - pkin(1);
	t262 = t263 * qJD(2) + qJD(3);
	t261 = t300 * t257 + t262 * t259;
	t239 = -t258 * t253 + t259 * t290;
	t238 = t253 * t291 - t290;
	t236 = t268 * qJD(1) - t253 * t276 - t259 * t274 - t272;
	t235 = -t252 * t276 - t253 * t288 + t267 * t259 - t273;
	t234 = t259 * t273 + (t259 * t288 + t275) * t253 - t267;
	t233 = t237 * qJD(1) + t252 * t275 - t259 * t272 - t274;
	t1 = [-t237 * qJD(5) - t238 * qJD(6) - t293 * t235 - t279 * t236 - t301 * t258 + (-t277 * t258 + t266 * t260) * qJD(1), t261 * t260 - t298 * t288, -t257 * t288 + t259 * t284, qJD(5) * t268 - t239 * qJD(6) + t279 * t233 - t293 * t234, -t233, -t234; t239 * qJD(5) + t268 * qJD(6) - t293 * t233 - t279 * t234 + t301 * t260 + (t266 * t258 + t277 * t260) * qJD(1), t261 * t258 + t298 * t287, t257 * t287 + t258 * t285, t238 * qJD(5) - t237 * qJD(6) - t279 * t235 + t293 * t236, t235, t236; 0, t262 * t257 - t300 * t259, t286, -t294 * t285 + (t265 * qJD(4) + qJD(5) * t253 - qJD(6) * t252) * t257, t252 * t285 + t253 * t283, -t252 * t283 + t253 * t285;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end