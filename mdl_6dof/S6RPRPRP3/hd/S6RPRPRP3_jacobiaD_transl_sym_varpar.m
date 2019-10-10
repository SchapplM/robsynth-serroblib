% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
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
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:30
	% EndTime: 2019-10-10 00:32:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:31
	% EndTime: 2019-10-10 00:32:32
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (113->22), mult. (182->34), div. (0->0), fcn. (135->8), ass. (0->19)
	t167 = sin(qJ(3));
	t168 = cos(qJ(3));
	t165 = sin(pkin(10));
	t166 = cos(pkin(10));
	t176 = r_i_i_C(1) * t166 - r_i_i_C(2) * t165 + pkin(3);
	t180 = r_i_i_C(3) + qJ(4);
	t173 = t176 * t167 - t180 * t168;
	t182 = qJD(1) * t173;
	t181 = t173 * qJD(3) - t167 * qJD(4);
	t179 = qJD(1) * t167;
	t178 = qJD(3) * t168;
	t175 = t165 * r_i_i_C(1) + t166 * r_i_i_C(2) + pkin(7);
	t174 = -t180 * t167 - t176 * t168;
	t171 = -pkin(2) + t174;
	t170 = t174 * qJD(3) + qJD(4) * t168;
	t164 = qJ(1) + pkin(9);
	t163 = cos(t164);
	t162 = sin(t164);
	t1 = [t181 * t162 + (-cos(qJ(1)) * pkin(1) - t175 * t162 + t171 * t163) * qJD(1), 0, t162 * t182 + t170 * t163, -t162 * t179 + t163 * t178, 0, 0; -t181 * t163 + (-sin(qJ(1)) * pkin(1) + t175 * t163 + t171 * t162) * qJD(1), 0, t170 * t162 - t163 * t182, t162 * t178 + t163 * t179, 0, 0; 0, 0, -t181, qJD(3) * t167, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:32
	% EndTime: 2019-10-10 00:32:32
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (262->46), mult. (307->74), div. (0->0), fcn. (244->10), ass. (0->36)
	t216 = sin(qJ(3));
	t207 = cos(pkin(10)) * pkin(4) + pkin(3);
	t212 = pkin(10) + qJ(5);
	t208 = sin(t212);
	t210 = cos(t212);
	t222 = r_i_i_C(1) * t210 - r_i_i_C(2) * t208 + t207;
	t217 = cos(qJ(3));
	t239 = r_i_i_C(3) + pkin(8) + qJ(4);
	t240 = t239 * t217;
	t219 = -t222 * t216 + t240;
	t244 = qJD(1) * t219;
	t231 = t216 * qJD(4);
	t243 = (-t207 * t216 + t240) * qJD(3) + t231;
	t213 = qJ(1) + pkin(9);
	t211 = cos(t213);
	t237 = t210 * t211;
	t236 = qJD(1) * t216;
	t235 = qJD(3) * t216;
	t234 = qJD(3) * t217;
	t233 = qJD(5) * t216;
	t232 = qJD(5) * t217;
	t230 = pkin(4) * sin(pkin(10)) + pkin(7);
	t229 = t239 * t216;
	t226 = -qJD(1) + t232;
	t225 = qJD(1) * t217 - qJD(5);
	t224 = r_i_i_C(1) * t208 + r_i_i_C(2) * t210;
	t223 = t226 * t208;
	t221 = -t207 * t217 - pkin(2) - t229;
	t209 = sin(t213);
	t220 = t225 * t209 + t211 * t235;
	t218 = qJD(4) * t217 + t224 * t233 + (-t222 * t217 - t229) * qJD(3);
	t206 = -t225 * t237 + (t210 * t235 + t223) * t209;
	t205 = t226 * t210 * t209 + (-t209 * t235 + t225 * t211) * t208;
	t204 = t220 * t210 + t211 * t223;
	t203 = t220 * t208 - t226 * t237;
	t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) - t243 * t209 + (-cos(qJ(1)) * pkin(1) - t230 * t209 + t221 * t211) * qJD(1), 0, -t209 * t244 + t218 * t211, -t209 * t236 + t211 * t234, t203 * r_i_i_C(1) + t204 * r_i_i_C(2), 0; -t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t243 * t211 + (-sin(qJ(1)) * pkin(1) + t230 * t211 + t221 * t209) * qJD(1), 0, t218 * t209 + t211 * t244, t209 * t234 + t211 * t236, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0; 0, 0, t219 * qJD(3) - t224 * t232 + t231, t235, (t208 * t233 - t210 * t234) * r_i_i_C(2) + (-t208 * t234 - t210 * t233) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:32:32
	% EndTime: 2019-10-10 00:32:33
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (499->64), mult. (547->96), div. (0->0), fcn. (460->10), ass. (0->46)
	t265 = sin(qJ(3));
	t256 = cos(pkin(10)) * pkin(4) + pkin(3);
	t261 = pkin(10) + qJ(5);
	t257 = sin(t261);
	t259 = cos(t261);
	t297 = r_i_i_C(3) + qJ(6);
	t299 = r_i_i_C(1) + pkin(5);
	t300 = t297 * t257 + t299 * t259 + t256;
	t266 = cos(qJ(3));
	t298 = r_i_i_C(2) + pkin(8) + qJ(4);
	t302 = t298 * t266;
	t306 = t300 * t265 - t302;
	t283 = t265 * qJD(4);
	t285 = qJD(6) * t257;
	t305 = (-t256 * t265 + t302) * qJD(3) + t266 * t285 + t283;
	t303 = -t285 + (t299 * t257 - t297 * t259) * qJD(5);
	t262 = qJ(1) + pkin(9);
	t258 = sin(t262);
	t295 = t258 * t266;
	t260 = cos(t262);
	t294 = t260 * t259;
	t293 = qJD(1) * t258;
	t292 = qJD(1) * t260;
	t291 = qJD(1) * t265;
	t290 = qJD(3) * t265;
	t289 = qJD(3) * t266;
	t288 = qJD(5) * t258;
	t287 = qJD(5) * t260;
	t286 = qJD(5) * t265;
	t284 = t259 * qJD(6);
	t282 = pkin(4) * sin(pkin(10)) + pkin(7);
	t281 = t258 * t290;
	t280 = t257 * t288;
	t279 = t260 * t290;
	t278 = t257 * t287;
	t277 = t259 * t287;
	t276 = t298 * t265;
	t273 = t258 * t257 + t266 * t294;
	t271 = -t256 * t266 - pkin(2) - t276;
	t270 = t257 * t292 + t259 * t288;
	t267 = qJD(4) * t266 + t303 * t265 + (-t266 * t300 - t276) * qJD(3);
	t245 = t273 * qJD(1) - t259 * t281 - t266 * t280 - t277;
	t244 = -t257 * t281 - t259 * t293 + t270 * t266 - t278;
	t243 = t266 * t278 + (t266 * t293 + t279) * t259 - t270;
	t242 = t257 * t279 - t266 * t277 - t280 + (t257 * t295 + t294) * qJD(1);
	t1 = [-t260 * t284 - t299 * t245 - t297 * t244 - t305 * t258 + (-cos(qJ(1)) * pkin(1) - t282 * t258 + t271 * t260) * qJD(1), 0, t267 * t260 + t306 * t293, -t258 * t291 + t260 * t289, t273 * qJD(6) + t299 * t242 - t297 * t243, -t242; -t258 * t284 - t299 * t243 - t297 * t242 + t305 * t260 + (-sin(qJ(1)) * pkin(1) + t282 * t260 + t271 * t258) * qJD(1), 0, t267 * t258 - t292 * t306, t258 * t289 + t260 * t291, -(t260 * t257 - t259 * t295) * qJD(6) + t297 * t245 - t299 * t244, t244; 0, 0, -qJD(3) * t306 - t303 * t266 + t283, t290, (-t297 * t286 - t299 * t289) * t257 + (t297 * t289 + (-t299 * qJD(5) + qJD(6)) * t265) * t259, t257 * t289 + t259 * t286;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end