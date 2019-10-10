% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
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
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(10) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(10)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:13
	% EndTime: 2019-10-10 01:32:13
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (160->43), mult. (276->76), div. (0->0), fcn. (217->7), ass. (0->33)
	t209 = pkin(10) + qJ(3);
	t207 = sin(t209);
	t208 = cos(t209);
	t234 = pkin(8) + r_i_i_C(3);
	t224 = t234 * t208;
	t235 = -pkin(3) * t207 + t224;
	t213 = cos(qJ(4));
	t214 = cos(qJ(1));
	t232 = t213 * t214;
	t212 = sin(qJ(1));
	t231 = qJD(1) * t212;
	t230 = qJD(1) * t214;
	t229 = qJD(3) * t212;
	t228 = qJD(3) * t213;
	t227 = qJD(3) * t214;
	t226 = qJD(4) * t207;
	t225 = qJD(4) * t208;
	t223 = -qJD(1) + t225;
	t222 = qJD(1) * t208 - qJD(4);
	t211 = sin(qJ(4));
	t221 = r_i_i_C(1) * t211 + r_i_i_C(2) * t213;
	t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(3);
	t219 = t223 * t211;
	t218 = -pkin(3) * t208 - t234 * t207 - cos(pkin(10)) * pkin(2) - pkin(1);
	t217 = qJD(3) * t220;
	t216 = t207 * t227 + t222 * t212;
	t215 = -t234 * qJD(3) + t221 * qJD(4);
	t210 = -pkin(7) - qJ(2);
	t205 = -t222 * t232 + (t207 * t228 + t219) * t212;
	t204 = t223 * t213 * t212 + (-t207 * t229 + t222 * t214) * t211;
	t203 = t216 * t213 + t214 * t219;
	t202 = t216 * t211 - t223 * t232;
	t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t214 * qJD(2) - t235 * t229 + (t210 * t212 + t218 * t214) * qJD(1), t230, (-t214 * t217 - t234 * t231) * t208 + (t215 * t214 + t220 * t231) * t207, t202 * r_i_i_C(1) + t203 * r_i_i_C(2), 0, 0; -t203 * r_i_i_C(1) + t202 * r_i_i_C(2) + t212 * qJD(2) + t235 * t227 + (-t210 * t214 + t218 * t212) * qJD(1), t231, (-t212 * t217 + t234 * t230) * t208 + (t215 * t212 - t220 * t230) * t207, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2), 0, 0; 0, 0, -t221 * t225 + (-t220 * t207 + t224) * qJD(3), (-t208 * t228 + t211 * t226) * r_i_i_C(2) + (-qJD(3) * t208 * t211 - t213 * t226) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:13
	% EndTime: 2019-10-10 01:32:13
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (289->50), mult. (376->75), div. (0->0), fcn. (296->9), ass. (0->41)
	t229 = cos(qJ(4));
	t257 = t229 * pkin(4);
	t218 = pkin(3) + t257;
	t223 = pkin(10) + qJ(3);
	t219 = sin(t223);
	t221 = cos(t223);
	t244 = qJD(4) * t221 - qJD(1);
	t249 = t219 * qJD(5);
	t227 = sin(qJ(4));
	t258 = pkin(4) * t227;
	t256 = r_i_i_C(3) + qJ(5) + pkin(8);
	t261 = t256 * t221;
	t264 = (-t218 * t219 + t261) * qJD(3) - t244 * t258 - qJD(1) * (-pkin(7) - qJ(2)) + t249;
	t224 = qJ(4) + pkin(11);
	t220 = sin(t224);
	t222 = cos(t224);
	t242 = r_i_i_C(1) * t222 - r_i_i_C(2) * t220;
	t238 = t218 + t242;
	t233 = -t238 * t219 + t261;
	t228 = sin(qJ(1));
	t240 = t244 * t228;
	t230 = cos(qJ(1));
	t241 = t244 * t230;
	t243 = qJD(1) * t221 - qJD(4);
	t252 = qJD(3) * t228;
	t260 = -t219 * t252 + t243 * t230;
	t254 = qJD(1) * t228;
	t253 = qJD(1) * t230;
	t251 = qJD(3) * t230;
	t250 = qJD(4) * t219;
	t247 = t256 * t219;
	t237 = r_i_i_C(1) * t220 + r_i_i_C(2) * t222 + t258;
	t236 = t237 * t221;
	t234 = t219 * t251 + t243 * t228;
	t232 = qJD(4) * t257 + qJD(2) + (-t218 * t221 - cos(pkin(10)) * pkin(2) - pkin(1) - t247) * qJD(1);
	t231 = qJD(5) * t221 + t237 * t250 + (-t238 * t221 - t247) * qJD(3);
	t216 = t220 * t240 - t222 * t260;
	t215 = t260 * t220 + t222 * t240;
	t214 = t220 * t241 + t234 * t222;
	t213 = t234 * t220 - t222 * t241;
	t1 = [t216 * r_i_i_C(1) + t215 * r_i_i_C(2) - t264 * t228 + t232 * t230, t253, t231 * t230 - t233 * t254, t213 * r_i_i_C(1) + t214 * r_i_i_C(2) + (t234 * t227 - t229 * t241) * pkin(4), -t219 * t254 + t221 * t251, 0; -t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t232 * t228 + t264 * t230, t254, t231 * t228 + t233 * t253, -t215 * r_i_i_C(1) + t216 * r_i_i_C(2) + (-t227 * t260 - t229 * t240) * pkin(4), t219 * t253 + t221 * t252, 0; 0, 0, t233 * qJD(3) - qJD(4) * t236 + t249, (-t242 - t257) * t250 - qJD(3) * t236, qJD(3) * t219, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:13
	% EndTime: 2019-10-10 01:32:14
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (537->62), mult. (486->87), div. (0->0), fcn. (386->11), ass. (0->53)
	t256 = qJ(4) + pkin(11);
	t241 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t256);
	t238 = t241 * qJD(4);
	t242 = pkin(5) * cos(t256) + cos(qJ(4)) * pkin(4);
	t240 = pkin(3) + t242;
	t254 = pkin(10) + qJ(3);
	t247 = sin(t254);
	t249 = cos(t254);
	t278 = t247 * qJD(5);
	t291 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t298 = t291 * t249;
	t302 = (-t240 * t247 + t298) * qJD(3) + (t241 + pkin(7) + qJ(2)) * qJD(1) - t249 * t238 + t278;
	t251 = qJ(6) + t256;
	t244 = sin(t251);
	t293 = r_i_i_C(2) * t244;
	t245 = cos(t251);
	t294 = r_i_i_C(1) * t245;
	t267 = t240 - t293 + t294;
	t264 = -t267 * t247 + t298;
	t261 = cos(qJ(1));
	t255 = qJD(4) + qJD(6);
	t273 = t249 * t255 - qJD(1);
	t300 = t261 * t273;
	t292 = r_i_i_C(2) * t245;
	t271 = r_i_i_C(1) * t244 + t292;
	t297 = t271 * t255 + t238;
	t285 = qJD(1) * t249;
	t272 = -t255 + t285;
	t259 = sin(qJ(1));
	t280 = qJD(3) * t259;
	t296 = -t247 * t280 + t272 * t261;
	t289 = t247 * t255;
	t279 = qJD(3) * t261;
	t265 = t247 * t279 + t272 * t259;
	t233 = t265 * t244 - t245 * t300;
	t234 = t244 * t300 + t265 * t245;
	t288 = t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
	t269 = t273 * t259;
	t235 = t296 * t244 + t245 * t269;
	t236 = t244 * t269 - t296 * t245;
	t287 = -t235 * r_i_i_C(1) + t236 * r_i_i_C(2);
	t284 = qJD(1) * t259;
	t283 = qJD(1) * t261;
	t282 = qJD(3) * t247;
	t281 = qJD(3) * t249;
	t276 = t291 * t247;
	t270 = t241 * t285 - t238;
	t239 = t242 * qJD(4);
	t266 = qJD(1) * t242 - t239 * t249 + t241 * t282;
	t263 = qJD(2) + t239 + (-t240 * t249 - cos(pkin(10)) * pkin(2) - pkin(1) - t276) * qJD(1);
	t262 = qJD(5) * t249 + t297 * t247 + (-t267 * t249 - t276) * qJD(3);
	t237 = t289 * t293;
	t1 = [t236 * r_i_i_C(1) + t235 * r_i_i_C(2) - t302 * t259 + t263 * t261, t283, t262 * t261 - t264 * t284, t270 * t259 + t266 * t261 + t288, -t247 * t284 + t249 * t279, t288; -t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t263 * t259 + t302 * t261, t284, t262 * t259 + t264 * t283, t266 * t259 - t270 * t261 + t287, t247 * t283 + t249 * t280, t287; 0, 0, t264 * qJD(3) - t297 * t249 + t278, t237 + (-t255 * t294 - t239) * t247 + (-t241 - t271) * t281, t282, -t281 * t292 + t237 + (-t244 * t281 - t245 * t289) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end