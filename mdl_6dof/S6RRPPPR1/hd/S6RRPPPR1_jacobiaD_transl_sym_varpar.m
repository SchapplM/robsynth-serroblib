% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:09
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.11s
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
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:10
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (123->22), mult. (204->32), div. (0->0), fcn. (152->8), ass. (0->19)
	t176 = sin(pkin(10));
	t177 = cos(pkin(10));
	t175 = qJ(2) + pkin(9);
	t173 = sin(t175);
	t174 = cos(t175);
	t190 = r_i_i_C(1) * t177 - r_i_i_C(2) * t176 + pkin(3);
	t195 = r_i_i_C(3) + qJ(4);
	t186 = t190 * t173 - t195 * t174 + sin(qJ(2)) * pkin(2);
	t183 = -t186 * qJD(2) + t173 * qJD(4);
	t200 = t183 + (t176 * r_i_i_C(1) + t177 * r_i_i_C(2) + pkin(7) + qJ(3)) * qJD(1);
	t199 = -t195 * t173 - t190 * t174 - cos(qJ(2)) * pkin(2);
	t180 = sin(qJ(1));
	t194 = qJD(1) * t180;
	t182 = cos(qJ(1));
	t193 = qJD(1) * t182;
	t192 = qJD(2) * t174;
	t185 = qJD(3) + (-pkin(1) + t199) * qJD(1);
	t184 = t199 * qJD(2) + qJD(4) * t174;
	t1 = [-t200 * t180 + t185 * t182, t184 * t182 + t186 * t194, t193, -t173 * t194 + t182 * t192, 0, 0; t185 * t180 + t200 * t182, t184 * t180 - t186 * t193, t194, t173 * t193 + t180 * t192, 0, 0; 0, t183, 0, qJD(2) * t173, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (184->41), mult. (312->62), div. (0->0), fcn. (248->8), ass. (0->32)
	t189 = qJ(2) + pkin(9);
	t187 = sin(t189);
	t188 = cos(t189);
	t190 = sin(pkin(10));
	t208 = qJD(5) * t190;
	t202 = t187 * qJD(4) + t188 * t208;
	t218 = r_i_i_C(2) + qJ(4);
	t204 = t218 * t188 - sin(qJ(2)) * pkin(2);
	t227 = (-pkin(3) * t187 + t204) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t202;
	t191 = cos(pkin(10));
	t217 = r_i_i_C(3) + qJ(5);
	t222 = r_i_i_C(1) + pkin(4);
	t223 = t217 * t190 + t222 * t191 + pkin(3);
	t225 = t223 * t187 - t204;
	t224 = -t218 * t187 - cos(qJ(2)) * pkin(2);
	t194 = sin(qJ(1));
	t216 = t190 * t194;
	t196 = cos(qJ(1));
	t215 = t190 * t196;
	t214 = t194 * t191;
	t213 = t196 * t191;
	t212 = qJD(1) * t194;
	t211 = qJD(1) * t196;
	t210 = qJD(2) * t194;
	t209 = qJD(2) * t196;
	t207 = t187 * t210;
	t206 = t187 * t209;
	t199 = -t191 * qJD(5) + qJD(3) + (-pkin(3) * t188 - pkin(1) + t224) * qJD(1);
	t197 = -t187 * t208 + qJD(4) * t188 + (-t188 * t223 + t224) * qJD(2);
	t184 = -t190 * t207 + (t188 * t215 - t214) * qJD(1);
	t182 = t190 * t206 + (t188 * t216 + t213) * qJD(1);
	t1 = [t222 * (t191 * t207 + (-t188 * t213 - t216) * qJD(1)) - t217 * t184 + t199 * t196 - t227 * t194, t197 * t196 + t225 * t212, t211, -t187 * t212 + t188 * t209, -t182, 0; t222 * (-t191 * t206 + (-t188 * t214 + t215) * qJD(1)) - t217 * t182 + t199 * t194 + t227 * t196, t197 * t194 - t211 * t225, t212, t187 * t211 + t188 * t210, t184, 0; 0, -qJD(2) * t225 + t202, 0, qJD(2) * t187, qJD(2) * t188 * t190, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:10
	% EndTime: 2019-10-10 09:18:11
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (387->70), mult. (685->110), div. (0->0), fcn. (619->10), ass. (0->50)
	t261 = qJ(2) + pkin(9);
	t259 = sin(t261);
	t260 = cos(t261);
	t289 = -r_i_i_C(3) - pkin(8) + qJ(4);
	t277 = t289 * t260 - sin(qJ(2)) * pkin(2);
	t291 = t259 * qJD(4);
	t310 = (-t259 * pkin(3) + t277) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t291;
	t262 = sin(pkin(10));
	t263 = cos(pkin(10));
	t265 = sin(qJ(6));
	t268 = cos(qJ(6));
	t304 = pkin(4) + pkin(5);
	t278 = t268 * r_i_i_C(1) - t265 * r_i_i_C(2) + t304;
	t279 = t265 * r_i_i_C(1) + t268 * r_i_i_C(2) + qJ(5);
	t305 = t279 * t262 + t278 * t263 + pkin(3);
	t308 = t305 * t259 - t277;
	t307 = -t289 * t259 - cos(qJ(2)) * pkin(2);
	t281 = t262 * t265 + t263 * t268;
	t282 = t262 * t268 - t263 * t265;
	t274 = t282 * r_i_i_C(1) - t281 * r_i_i_C(2);
	t306 = t262 * qJD(5) + t274 * qJD(6);
	t267 = sin(qJ(1));
	t300 = t267 * t262;
	t299 = t267 * t263;
	t270 = cos(qJ(1));
	t298 = t270 * t262;
	t297 = t270 * t263;
	t296 = qJD(1) * t267;
	t295 = qJD(1) * t270;
	t294 = qJD(2) * t260;
	t293 = qJD(2) * t267;
	t292 = qJD(2) * t270;
	t288 = t260 * t298;
	t287 = t259 * t293;
	t286 = t259 * t292;
	t250 = t260 * t300 + t297;
	t251 = t260 * t299 - t298;
	t284 = -t250 * t268 + t251 * t265;
	t283 = t250 * t265 + t251 * t268;
	t253 = t260 * t297 + t300;
	t275 = qJD(3) + (-pkin(3) * t260 - pkin(1) + t307) * qJD(1);
	t271 = t260 * qJD(4) - t306 * t259 + (-t260 * t305 + t307) * qJD(2);
	t252 = t288 - t299;
	t249 = t253 * qJD(1) - t263 * t287;
	t248 = qJD(1) * t288 - t262 * t287 - t263 * t296;
	t247 = -t251 * qJD(1) - t263 * t286;
	t246 = t250 * qJD(1) + t262 * t286;
	t245 = -t246 * t265 + t247 * t268 + (t252 * t268 - t253 * t265) * qJD(6);
	t244 = -t246 * t268 - t247 * t265 + (-t252 * t265 - t253 * t268) * qJD(6);
	t1 = [-t250 * qJD(5) - t279 * t248 - t278 * t249 + (t284 * r_i_i_C(1) + t283 * r_i_i_C(2)) * qJD(6) + t275 * t270 - t310 * t267, t271 * t270 + t308 * t296, t295, -t259 * t296 + t260 * t292, -t246, t244 * r_i_i_C(1) - t245 * r_i_i_C(2); t245 * r_i_i_C(1) + t244 * r_i_i_C(2) - t246 * qJ(5) + t252 * qJD(5) + t304 * t247 + t275 * t267 + t310 * t270, t271 * t267 - t295 * t308, t296, t259 * t295 + t260 * t293, t248, (t248 * t268 - t249 * t265) * r_i_i_C(1) + (-t248 * t265 - t249 * t268) * r_i_i_C(2) + (-t283 * r_i_i_C(1) + t284 * r_i_i_C(2)) * qJD(6); 0, -qJD(2) * t308 + t306 * t260 + t291, 0, qJD(2) * t259, t262 * t294, (-t281 * r_i_i_C(1) - t282 * r_i_i_C(2)) * t259 * qJD(6) + t274 * t294;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end