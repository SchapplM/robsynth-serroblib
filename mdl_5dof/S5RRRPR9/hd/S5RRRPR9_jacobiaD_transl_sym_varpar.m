% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:10:51
	% EndTime: 2019-12-29 20:10:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:10:51
	% EndTime: 2019-12-29 20:10:51
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
	% StartTime: 2019-12-29 20:10:51
	% EndTime: 2019-12-29 20:10:51
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
	% StartTime: 2019-12-29 20:10:47
	% EndTime: 2019-12-29 20:10:48
	% DurationCPUTime: 0.30s
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
	% StartTime: 2019-12-29 20:10:53
	% EndTime: 2019-12-29 20:10:53
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (189->46), mult. (370->74), div. (0->0), fcn. (290->8), ass. (0->40)
	t213 = cos(qJ(3));
	t242 = t213 * pkin(3);
	t205 = pkin(2) + t242;
	t211 = sin(qJ(2));
	t214 = cos(qJ(2));
	t229 = qJD(3) * t214 - qJD(1);
	t234 = t211 * qJD(4);
	t210 = sin(qJ(3));
	t243 = pkin(3) * t210;
	t241 = r_i_i_C(3) + qJ(4) + pkin(7);
	t246 = t241 * t214;
	t249 = pkin(6) * qJD(1) + (-t205 * t211 + t246) * qJD(2) - t229 * t243 + t234;
	t208 = qJ(3) + pkin(9);
	t206 = sin(t208);
	t207 = cos(t208);
	t227 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
	t223 = t205 + t227;
	t218 = -t223 * t211 + t246;
	t212 = sin(qJ(1));
	t224 = t229 * t212;
	t215 = cos(qJ(1));
	t225 = t229 * t215;
	t228 = qJD(1) * t214 - qJD(3);
	t237 = qJD(2) * t212;
	t245 = -t211 * t237 + t228 * t215;
	t239 = qJD(1) * t212;
	t238 = qJD(1) * t215;
	t236 = qJD(2) * t215;
	t235 = qJD(3) * t211;
	t232 = t241 * t211;
	t222 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t243;
	t221 = t222 * t214;
	t219 = t211 * t236 + t228 * t212;
	t217 = qJD(3) * t242 + (-t205 * t214 - pkin(1) - t232) * qJD(1);
	t216 = qJD(4) * t214 + t222 * t235 + (-t223 * t214 - t232) * qJD(2);
	t204 = t206 * t224 - t207 * t245;
	t203 = t245 * t206 + t207 * t224;
	t202 = t206 * t225 + t219 * t207;
	t201 = t219 * t206 - t207 * t225;
	t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t249 * t212 + t217 * t215, t216 * t215 - t218 * t239, t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (t219 * t210 - t213 * t225) * pkin(3), -t211 * t239 + t214 * t236, 0; -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t217 * t212 + t249 * t215, t216 * t212 + t218 * t238, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2) + (-t210 * t245 - t213 * t224) * pkin(3), t211 * t238 + t214 * t237, 0; 0, t218 * qJD(2) - qJD(3) * t221 + t234, (-t227 - t242) * t235 - qJD(2) * t221, qJD(2) * t211, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:10:48
	% EndTime: 2019-12-29 20:10:48
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (421->58), mult. (480->85), div. (0->0), fcn. (380->10), ass. (0->51)
	t243 = qJ(3) + pkin(9);
	t232 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t243);
	t228 = t232 * qJD(3);
	t233 = pkin(4) * cos(t243) + cos(qJ(3)) * pkin(3);
	t231 = pkin(2) + t233;
	t245 = sin(qJ(2));
	t248 = cos(qJ(2));
	t266 = t245 * qJD(4);
	t277 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
	t285 = t277 * t248;
	t289 = (-t231 * t245 + t285) * qJD(2) + (pkin(6) + t232) * qJD(1) - t248 * t228 + t266;
	t239 = qJ(5) + t243;
	t235 = sin(t239);
	t280 = r_i_i_C(2) * t235;
	t236 = cos(t239);
	t281 = r_i_i_C(1) * t236;
	t255 = t231 - t280 + t281;
	t251 = -t255 * t245 + t285;
	t249 = cos(qJ(1));
	t242 = qJD(3) + qJD(5);
	t261 = t242 * t248 - qJD(1);
	t287 = t249 * t261;
	t279 = r_i_i_C(2) * t236;
	t259 = r_i_i_C(1) * t235 + t279;
	t284 = t259 * t242 + t228;
	t246 = sin(qJ(1));
	t271 = qJD(1) * t248;
	t260 = -t242 + t271;
	t269 = qJD(2) * t245;
	t283 = -t246 * t269 + t260 * t249;
	t275 = t242 * t245;
	t267 = qJD(2) * t249;
	t253 = t245 * t267 + t260 * t246;
	t224 = t253 * t235 - t236 * t287;
	t225 = t235 * t287 + t253 * t236;
	t274 = t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t257 = t261 * t246;
	t226 = t283 * t235 + t236 * t257;
	t227 = t235 * t257 - t283 * t236;
	t273 = -t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
	t272 = qJD(1) * t246;
	t270 = qJD(1) * t249;
	t268 = qJD(2) * t248;
	t264 = t277 * t245;
	t258 = t232 * t271 - t228;
	t229 = t233 * qJD(3);
	t254 = qJD(1) * t233 - t229 * t248 + t232 * t269;
	t252 = t229 + (-t231 * t248 - pkin(1) - t264) * qJD(1);
	t250 = qJD(4) * t248 + t284 * t245 + (-t255 * t248 - t264) * qJD(2);
	t230 = t275 * t280;
	t1 = [t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t289 * t246 + t252 * t249, t250 * t249 - t251 * t272, t258 * t246 + t254 * t249 + t274, -t245 * t272 + t248 * t267, t274; -t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + t252 * t246 + t289 * t249, t250 * t246 + t251 * t270, t254 * t246 - t258 * t249 + t273, t245 * t270 + t246 * t268, t273; 0, t251 * qJD(2) - t284 * t248 + t266, t230 + (-t242 * t281 - t229) * t245 + (-t232 - t259) * t268, t269, -t268 * t279 + t230 + (-t235 * t268 - t236 * t275) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end