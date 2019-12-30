% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:57
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
	% StartTime: 2019-12-29 19:53:57
	% EndTime: 2019-12-29 19:53:58
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
	% StartTime: 2019-12-29 19:53:58
	% EndTime: 2019-12-29 19:53:59
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(7) + r_i_i_C(3);
	t213 = t223 * t202;
	t224 = -pkin(2) * t199 + t213;
	t201 = cos(qJ(3));
	t203 = cos(qJ(1));
	t221 = t201 * t203;
	t200 = sin(qJ(1));
	t220 = qJD(1) * t200;
	t219 = qJD(1) * t203;
	t218 = qJD(2) * t200;
	t217 = qJD(2) * t202;
	t216 = qJD(2) * t203;
	t215 = qJD(3) * t199;
	t214 = qJD(3) * t202;
	t212 = -qJD(1) + t214;
	t211 = qJD(1) * t202 - qJD(3);
	t198 = sin(qJ(3));
	t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
	t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
	t208 = t212 * t198;
	t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
	t206 = qJD(2) * t209;
	t205 = t199 * t216 + t211 * t200;
	t204 = -t223 * qJD(2) + t210 * qJD(3);
	t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
	t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
	t195 = t205 * t201 + t203 * t208;
	t194 = t205 * t198 - t212 * t221;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(6) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(6) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:53:59
	% EndTime: 2019-12-29 19:54:00
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (164->54), mult. (510->88), div. (0->0), fcn. (427->6), ass. (0->42)
	t244 = sin(qJ(3));
	t247 = cos(qJ(3));
	t277 = r_i_i_C(3) + qJ(4);
	t279 = -r_i_i_C(2) + pkin(3);
	t281 = t244 * t279 - t247 * t277;
	t283 = -t281 * qJD(3) + qJD(4) * t244;
	t245 = sin(qJ(2));
	t248 = cos(qJ(2));
	t280 = pkin(7) + r_i_i_C(1);
	t264 = t280 * t248;
	t282 = -pkin(2) * t245 + t264;
	t255 = -t244 * t277 - t247 * t279;
	t252 = -pkin(2) + t255;
	t246 = sin(qJ(1));
	t276 = t246 * t244;
	t275 = t246 * t248;
	t249 = cos(qJ(1));
	t274 = t249 * t244;
	t273 = t249 * t247;
	t272 = qJD(1) * t246;
	t271 = qJD(1) * t249;
	t270 = qJD(2) * t246;
	t269 = qJD(2) * t248;
	t268 = qJD(2) * t249;
	t267 = qJD(3) * t247;
	t266 = qJD(3) * t249;
	t263 = t245 * t270;
	t262 = qJD(3) * t276;
	t261 = t245 * t268;
	t260 = t244 * t266;
	t259 = t247 * t266;
	t258 = t248 * t273 + t276;
	t257 = t244 * t275 + t273;
	t256 = -pkin(2) * t248 - t245 * t280 - pkin(1);
	t253 = t244 * t271 + t246 * t267;
	t251 = qJD(2) * t252;
	t250 = -t280 * qJD(2) - t283;
	t233 = qJD(1) * t258 - t247 * t263 - t248 * t262 - t259;
	t232 = -t244 * t263 - t247 * t272 + t248 * t253 - t260;
	t231 = t248 * t260 + (t248 * t272 + t261) * t247 - t253;
	t230 = qJD(1) * t257 + t244 * t261 - t248 * t259 - t262;
	t1 = [-t257 * qJD(4) - t279 * t233 - t277 * t232 - t282 * t270 + (-t246 * pkin(6) + t249 * t256) * qJD(1), (t249 * t251 - t272 * t280) * t248 + (t250 * t249 - t252 * t272) * t245, qJD(4) * t258 + t230 * t279 - t231 * t277, -t230, 0; -(t246 * t247 - t248 * t274) * qJD(4) - t279 * t231 - t277 * t230 + t282 * t268 + (t249 * pkin(6) + t246 * t256) * qJD(1), (t246 * t251 + t271 * t280) * t248 + (t246 * t250 + t252 * t271) * t245, -(-t247 * t275 + t274) * qJD(4) + t277 * t233 - t279 * t232, t232, 0; 0, t283 * t248 + (t245 * t252 + t264) * qJD(2), -t281 * t269 + (qJD(3) * t255 + t247 * qJD(4)) * t245, t244 * t269 + t245 * t267, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:54:00
	% EndTime: 2019-12-29 19:54:00
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (224->61), mult. (684->96), div. (0->0), fcn. (582->6), ass. (0->42)
	t245 = sin(qJ(3));
	t248 = cos(qJ(3));
	t266 = pkin(4) + pkin(7) + r_i_i_C(1);
	t265 = pkin(3) + r_i_i_C(3) + qJ(5);
	t278 = r_i_i_C(2) + qJ(4);
	t280 = t265 * t245 - t278 * t248;
	t251 = -t266 * qJD(2) + t280 * qJD(3) - qJD(4) * t245 - qJD(5) * t248;
	t255 = -t278 * t245 - t265 * t248;
	t253 = -pkin(2) + t255;
	t246 = sin(qJ(2));
	t249 = cos(qJ(2));
	t282 = -pkin(2) * t246 + t266 * t249;
	t247 = sin(qJ(1));
	t277 = t247 * t249;
	t250 = cos(qJ(1));
	t276 = t250 * t245;
	t275 = t250 * t248;
	t274 = qJD(1) * t247;
	t273 = qJD(1) * t250;
	t272 = qJD(2) * t247;
	t271 = qJD(2) * t249;
	t270 = qJD(2) * t250;
	t269 = qJD(3) * t245;
	t268 = qJD(3) * t248;
	t267 = qJD(3) * t250;
	t264 = t246 * t272;
	t263 = t246 * t270;
	t262 = t247 * t269;
	t261 = t245 * t267;
	t260 = t248 * t267;
	t258 = t247 * t245 + t249 * t275;
	t231 = t245 * t277 + t275;
	t257 = -pkin(2) * t249 - t266 * t246 - pkin(1);
	t256 = t245 * t273 + t247 * t268;
	t252 = qJD(2) * t253;
	t233 = -t247 * t248 + t249 * t276;
	t232 = t248 * t277 - t276;
	t230 = t258 * qJD(1) - t248 * t264 - t249 * t262 - t260;
	t229 = -t245 * t264 - t248 * t274 + t256 * t249 - t261;
	t228 = t249 * t261 + (t249 * t274 + t263) * t248 - t256;
	t227 = t231 * qJD(1) + t245 * t263 - t249 * t260 - t262;
	t1 = [-t231 * qJD(4) - t232 * qJD(5) - t278 * t229 - t265 * t230 - t282 * t272 + (-pkin(6) * t247 + t257 * t250) * qJD(1), (t250 * t252 - t266 * t274) * t249 + (t251 * t250 - t253 * t274) * t246, qJD(4) * t258 - t233 * qJD(5) + t265 * t227 - t278 * t228, -t227, -t228; t233 * qJD(4) + t258 * qJD(5) - t278 * t227 - t265 * t228 + t282 * t270 + (pkin(6) * t250 + t257 * t247) * qJD(1), (t247 * t252 + t266 * t273) * t249 + (t251 * t247 + t253 * t273) * t246, t232 * qJD(4) - t231 * qJD(5) - t265 * t229 + t278 * t230, t229, t230; 0, t246 * t252 - t249 * t251, -t280 * t271 + (t255 * qJD(3) + qJD(4) * t248 - qJD(5) * t245) * t246, t245 * t271 + t246 * t268, -t246 * t269 + t248 * t271;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end