% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP9
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
% Datum: 2019-12-29 18:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:55:59
	% EndTime: 2019-12-29 18:55:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:55:50
	% EndTime: 2019-12-29 18:55:50
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
	% StartTime: 2019-12-29 18:55:42
	% EndTime: 2019-12-29 18:55:43
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
	% StartTime: 2019-12-29 18:56:01
	% EndTime: 2019-12-29 18:56:01
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (55->22), mult. (178->38), div. (0->0), fcn. (133->6), ass. (0->17)
	t160 = sin(qJ(2));
	t162 = cos(qJ(2));
	t158 = sin(pkin(8));
	t159 = cos(pkin(8));
	t169 = r_i_i_C(1) * t159 - r_i_i_C(2) * t158 + pkin(2);
	t174 = r_i_i_C(3) + qJ(3);
	t175 = t169 * t160 - t174 * t162;
	t176 = t175 * qJD(2) - t160 * qJD(3);
	t161 = sin(qJ(1));
	t173 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t172 = qJD(1) * t163;
	t171 = qJD(2) * t163;
	t168 = t158 * r_i_i_C(1) + t159 * r_i_i_C(2) + pkin(6);
	t167 = -t174 * t160 - t169 * t162;
	t165 = -pkin(1) + t167;
	t1 = [t176 * t161 + (-t168 * t161 + t165 * t163) * qJD(1), (t169 * t173 - t174 * t171) * t160 + (-t174 * t173 + (-t169 * qJD(2) + qJD(3)) * t163) * t162, -t160 * t173 + t162 * t171, 0, 0; -t176 * t163 + (t165 * t161 + t168 * t163) * qJD(1), -t175 * t172 + (t167 * qJD(2) + qJD(3) * t162) * t161, t161 * qJD(2) * t162 + t160 * t172, 0, 0; 0, -t176, qJD(2) * t160, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:56:06
	% EndTime: 2019-12-29 18:56:07
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (170->42), mult. (303->71), div. (0->0), fcn. (242->8), ass. (0->37)
	t203 = cos(pkin(8)) * pkin(3) + pkin(2);
	t209 = sin(qJ(2));
	t227 = t209 * qJD(3);
	t211 = cos(qJ(2));
	t236 = r_i_i_C(3) + pkin(7) + qJ(3);
	t238 = t236 * t211;
	t242 = (-t203 * t209 + t238) * qJD(2) + t227;
	t206 = pkin(8) + qJ(4);
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
	t226 = pkin(3) * sin(pkin(8)) + pkin(6);
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
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t242 * t210 + (-t226 * t210 + t216 * t212) * qJD(1), t213 * t212 - t214 * t234, -t209 * t234 + t211 * t230, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t242 * t212 + (t216 * t210 + t226 * t212) * qJD(1), t213 * t210 + t214 * t233, t209 * t233 + t210 * t231, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0; 0, t214 * qJD(2) - t219 * t228 + t227, t232, (t204 * t229 - t205 * t231) * r_i_i_C(2) + (-t204 * t231 - t205 * t229) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:56:02
	% EndTime: 2019-12-29 18:56:03
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (335->61), mult. (543->94), div. (0->0), fcn. (458->8), ass. (0->45)
	t259 = sin(qJ(2));
	t253 = cos(pkin(8)) * pkin(3) + pkin(2);
	t256 = pkin(8) + qJ(4);
	t254 = sin(t256);
	t255 = cos(t256);
	t293 = r_i_i_C(3) + qJ(5);
	t295 = r_i_i_C(1) + pkin(4);
	t296 = t293 * t254 + t295 * t255 + t253;
	t261 = cos(qJ(2));
	t294 = r_i_i_C(2) + pkin(7) + qJ(3);
	t298 = t294 * t261;
	t302 = t296 * t259 - t298;
	t279 = t259 * qJD(3);
	t281 = qJD(5) * t254;
	t301 = (-t253 * t259 + t298) * qJD(2) + t261 * t281 + t279;
	t299 = -t281 + (t295 * t254 - t293 * t255) * qJD(4);
	t260 = sin(qJ(1));
	t291 = t260 * t261;
	t262 = cos(qJ(1));
	t290 = t262 * t255;
	t289 = qJD(1) * t260;
	t288 = qJD(1) * t262;
	t287 = qJD(2) * t259;
	t286 = qJD(2) * t261;
	t285 = qJD(2) * t262;
	t284 = qJD(4) * t259;
	t283 = qJD(4) * t260;
	t282 = qJD(4) * t262;
	t280 = t255 * qJD(5);
	t278 = pkin(3) * sin(pkin(8)) + pkin(6);
	t277 = t260 * t287;
	t276 = t254 * t283;
	t275 = t259 * t285;
	t274 = t254 * t282;
	t273 = t255 * t282;
	t272 = t294 * t259;
	t269 = t260 * t254 + t261 * t290;
	t267 = -t253 * t261 - pkin(1) - t272;
	t266 = t254 * t288 + t255 * t283;
	t263 = qJD(3) * t261 + t299 * t259 + (-t261 * t296 - t272) * qJD(2);
	t242 = t269 * qJD(1) - t255 * t277 - t261 * t276 - t273;
	t241 = -t254 * t277 - t255 * t289 + t266 * t261 - t274;
	t240 = t261 * t274 + (t261 * t289 + t275) * t255 - t266;
	t239 = t254 * t275 - t261 * t273 - t276 + (t254 * t291 + t290) * qJD(1);
	t1 = [-t262 * t280 - t295 * t242 - t293 * t241 - t301 * t260 + (-t278 * t260 + t267 * t262) * qJD(1), t263 * t262 + t302 * t289, -t259 * t289 + t261 * t285, t269 * qJD(5) + t295 * t239 - t293 * t240, -t239; -t260 * t280 - t295 * t240 - t293 * t239 + t301 * t262 + (t267 * t260 + t278 * t262) * qJD(1), t263 * t260 - t288 * t302, t259 * t288 + t260 * t286, -(t262 * t254 - t255 * t291) * qJD(5) + t293 * t242 - t295 * t241, t241; 0, -qJD(2) * t302 - t299 * t261 + t279, t287, (-t293 * t284 - t295 * t286) * t254 + (t293 * t286 + (-t295 * qJD(4) + qJD(5)) * t259) * t255, t254 * t286 + t255 * t284;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end