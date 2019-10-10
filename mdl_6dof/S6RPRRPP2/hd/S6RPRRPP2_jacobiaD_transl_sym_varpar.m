% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
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
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
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
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 01:11:27
	% EndTime: 2019-10-10 01:11:28
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (165->37), mult. (274->63), div. (0->0), fcn. (213->8), ass. (0->32)
	t207 = sin(qJ(3));
	t206 = sin(qJ(4));
	t208 = cos(qJ(4));
	t215 = r_i_i_C(1) * t208 - r_i_i_C(2) * t206 + pkin(3);
	t209 = cos(qJ(3));
	t228 = pkin(8) + r_i_i_C(3);
	t230 = t228 * t209;
	t211 = -t215 * t207 + t230;
	t235 = qJD(1) * t211;
	t234 = (-pkin(3) * t207 + t230) * qJD(3);
	t217 = qJD(1) * t209 - qJD(4);
	t232 = t208 * t217;
	t223 = qJD(4) * t209;
	t218 = -qJD(1) + t223;
	t226 = qJD(3) * t207;
	t229 = -t206 * t226 + t218 * t208;
	t225 = qJD(3) * t209;
	t224 = qJD(4) * t207;
	t222 = t228 * t207;
	t216 = r_i_i_C(1) * t206 + r_i_i_C(2) * t208;
	t214 = t217 * t206;
	t213 = -pkin(3) * t209 - pkin(2) - t222;
	t212 = t218 * t206 + t208 * t226;
	t210 = t216 * t224 + (-t215 * t209 - t222) * qJD(3);
	t205 = qJ(1) + pkin(9);
	t204 = cos(t205);
	t203 = sin(t205);
	t202 = t212 * t203 - t204 * t232;
	t201 = t229 * t203 + t204 * t214;
	t200 = t203 * t232 + t212 * t204;
	t199 = t203 * t214 - t229 * t204;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, 0, t211 * qJD(3) - t216 * t223, (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:28
	% EndTime: 2019-10-10 01:11:28
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
	t250 = sin(qJ(3));
	t252 = cos(qJ(3));
	t281 = pkin(8) + r_i_i_C(2);
	t284 = t281 * t252;
	t287 = (-pkin(3) * t250 + t284) * qJD(3);
	t249 = sin(qJ(4));
	t251 = cos(qJ(4));
	t279 = r_i_i_C(3) + qJ(5);
	t282 = pkin(4) + r_i_i_C(1);
	t283 = t282 * t249 - t279 * t251;
	t286 = t283 * qJD(4) - qJD(5) * t249;
	t258 = -t279 * t249 - t282 * t251;
	t255 = -pkin(3) + t258;
	t254 = t255 * t250 + t284;
	t278 = t249 * t252;
	t277 = t251 * t252;
	t248 = qJ(1) + pkin(9);
	t246 = sin(t248);
	t276 = qJD(1) * t246;
	t247 = cos(t248);
	t275 = qJD(1) * t247;
	t274 = qJD(3) * t250;
	t273 = qJD(3) * t252;
	t272 = qJD(4) * t249;
	t271 = qJD(4) * t251;
	t269 = t281 * t250;
	t266 = t251 * t276;
	t265 = t249 * t274;
	t264 = t251 * t274;
	t263 = t246 * t272;
	t262 = t247 * t271;
	t261 = t246 * t249 + t247 * t277;
	t260 = t246 * t278 + t247 * t251;
	t259 = -pkin(3) * t252 - pkin(2) - t269;
	t256 = t246 * t271 + t249 * t275;
	t253 = t286 * t250 + (t255 * t252 - t269) * qJD(3);
	t235 = t261 * qJD(1) - t246 * t264 - t252 * t263 - t262;
	t234 = -t246 * t265 - t247 * t272 + t256 * t252 - t266;
	t233 = t252 * t266 + (t252 * t272 + t264) * t247 - t256;
	t232 = t260 * qJD(1) + t247 * t265 - t252 * t262 - t263;
	t1 = [-t260 * qJD(5) - t282 * t235 - t279 * t234 - t246 * t287 + (-cos(qJ(1)) * pkin(1) - t246 * pkin(7) + t259 * t247) * qJD(1), 0, t253 * t247 - t254 * t276, t261 * qJD(5) + t282 * t232 - t279 * t233, -t232, 0; -(t246 * t251 - t247 * t278) * qJD(5) - t282 * t233 - t279 * t232 + t247 * t287 + (-sin(qJ(1)) * pkin(1) + t247 * pkin(7) + t259 * t246) * qJD(1), 0, t253 * t246 + t254 * t275, -(-t246 * t277 + t247 * t249) * qJD(5) + t279 * t235 - t282 * t234, t234, 0; 0, 0, t254 * qJD(3) - t286 * t252, -t283 * t273 + (t258 * qJD(4) + t251 * qJD(5)) * t250, t249 * t273 + t250 * t271, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:27
	% EndTime: 2019-10-10 01:11:28
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (407->60), mult. (650->89), div. (0->0), fcn. (542->8), ass. (0->44)
	t221 = sin(qJ(3));
	t223 = cos(qJ(3));
	t244 = t221 * qJD(6);
	t220 = sin(qJ(4));
	t245 = qJD(5) * t220;
	t241 = pkin(8) - r_i_i_C(3) - qJ(6);
	t258 = t241 * t223;
	t261 = (-pkin(3) * t221 + t258) * qJD(3) + t223 * t245 - t244;
	t222 = cos(qJ(4));
	t242 = pkin(4) + pkin(5) + r_i_i_C(1);
	t255 = r_i_i_C(2) + qJ(5);
	t257 = t242 * t220 - t255 * t222;
	t260 = t257 * qJD(4) - t245;
	t228 = -t255 * t220 - t242 * t222;
	t226 = -pkin(3) + t228;
	t225 = t226 * t221 + t258;
	t219 = qJ(1) + pkin(9);
	t217 = sin(t219);
	t254 = t217 * t220;
	t253 = t222 * t223;
	t252 = qJD(1) * t217;
	t218 = cos(t219);
	t251 = qJD(1) * t218;
	t250 = qJD(1) * t221;
	t249 = qJD(3) * t221;
	t248 = qJD(3) * t223;
	t247 = qJD(4) * t220;
	t246 = qJD(4) * t222;
	t243 = t222 * qJD(5);
	t240 = t222 * t252;
	t239 = t220 * t249;
	t238 = t222 * t249;
	t237 = t217 * t247;
	t236 = t218 * t246;
	t235 = t241 * t221;
	t232 = t218 * t253 + t254;
	t230 = t217 * t246 + t220 * t251;
	t229 = -pkin(3) * t223 - pkin(2) - t235;
	t224 = -qJD(6) * t223 + t260 * t221 + (t226 * t223 - t235) * qJD(3);
	t206 = t232 * qJD(1) - t217 * t238 - t223 * t237 - t236;
	t205 = -t217 * t239 - t218 * t247 + t230 * t223 - t240;
	t204 = t223 * t240 + (t223 * t247 + t238) * t218 - t230;
	t203 = t218 * t239 - t223 * t236 - t237 + (t218 * t222 + t223 * t254) * qJD(1);
	t1 = [-t218 * t243 - t255 * t205 - t242 * t206 - t261 * t217 + (-cos(qJ(1)) * pkin(1) - t217 * pkin(7) + t229 * t218) * qJD(1), 0, t224 * t218 - t225 * t252, t232 * qJD(5) + t242 * t203 - t255 * t204, -t203, t217 * t250 - t218 * t248; -t217 * t243 - t255 * t203 - t242 * t204 + t261 * t218 + (-sin(qJ(1)) * pkin(1) + t218 * pkin(7) + t229 * t217) * qJD(1), 0, t224 * t217 + t225 * t251, -(-t217 * t253 + t218 * t220) * qJD(5) + t255 * t206 - t242 * t205, t205, -t217 * t248 - t218 * t250; 0, 0, t225 * qJD(3) - t260 * t223 - t244, -t257 * t248 + (t228 * qJD(4) + t243) * t221, t220 * t248 + t221 * t246, -t249;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end