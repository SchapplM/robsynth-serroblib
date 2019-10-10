% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR2
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
% Datum: 2019-10-10 01:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
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
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
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
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 0.20s
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
	t205 = qJ(1) + pkin(10);
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
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (297->54), mult. (374->82), div. (0->0), fcn. (292->10), ass. (0->46)
	t218 = sin(qJ(3));
	t219 = cos(qJ(4));
	t249 = t219 * pkin(4);
	t209 = pkin(3) + t249;
	t214 = qJ(4) + pkin(11);
	t210 = sin(t214);
	t212 = cos(t214);
	t233 = r_i_i_C(1) * t212 - r_i_i_C(2) * t210;
	t229 = t209 + t233;
	t220 = cos(qJ(3));
	t248 = r_i_i_C(3) + qJ(5) + pkin(8);
	t251 = t248 * t220;
	t222 = -t229 * t218 + t251;
	t256 = qJD(1) * t222;
	t241 = t218 * qJD(5);
	t242 = qJD(4) * t220;
	t217 = sin(qJ(4));
	t250 = pkin(4) * t217;
	t255 = (-t209 * t218 + t251) * qJD(3) - t242 * t250 + t241;
	t215 = qJ(1) + pkin(10);
	t213 = cos(t215);
	t234 = qJD(1) * t220 - qJD(4);
	t231 = t234 * t213;
	t235 = -qJD(1) + t242;
	t252 = t235 * t212;
	t246 = qJD(1) * t218;
	t245 = qJD(3) * t218;
	t244 = qJD(3) * t220;
	t243 = qJD(4) * t218;
	t240 = qJD(4) * t249;
	t239 = pkin(7) + t250;
	t238 = t248 * t218;
	t232 = t235 * t210;
	t211 = sin(t215);
	t230 = t234 * t211;
	t228 = r_i_i_C(1) * t210 + r_i_i_C(2) * t212 + t250;
	t227 = -t209 * t220 - pkin(2) - t238;
	t225 = t228 * t220;
	t224 = t213 * t245 + t230;
	t223 = t217 * t245 - t235 * t219;
	t221 = qJD(5) * t220 + t228 * t243 + (-t229 * t220 - t238) * qJD(3);
	t208 = -t212 * t231 + (t212 * t245 + t232) * t211;
	t207 = t211 * t252 + (-t211 * t245 + t231) * t210;
	t206 = t224 * t212 + t213 * t232;
	t205 = t224 * t210 - t213 * t252;
	t1 = [t213 * t240 + t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t255 * t211 + (-cos(qJ(1)) * pkin(1) - t239 * t211 + t227 * t213) * qJD(1), 0, -t211 * t256 + t221 * t213, t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t223 * t213 + t217 * t230) * pkin(4), -t211 * t246 + t213 * t244, 0; t211 * t240 - t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t255 * t213 + (-sin(qJ(1)) * pkin(1) + t239 * t213 + t227 * t211) * qJD(1), 0, t221 * t211 + t213 * t256, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2) + (t223 * t211 - t217 * t231) * pkin(4), t211 * t244 + t213 * t246, 0; 0, 0, t222 * qJD(3) - qJD(4) * t225 + t241, (-t233 - t249) * t243 - qJD(3) * t225, t245, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:15
	% EndTime: 2019-10-10 01:25:15
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (549->63), mult. (484->90), div. (0->0), fcn. (382->12), ass. (0->52)
	t250 = sin(qJ(3));
	t247 = qJ(4) + pkin(11);
	t235 = pkin(5) * cos(t247) + cos(qJ(4)) * pkin(4);
	t233 = pkin(3) + t235;
	t243 = qJ(6) + t247;
	t237 = sin(t243);
	t281 = r_i_i_C(2) * t237;
	t238 = cos(t243);
	t282 = r_i_i_C(1) * t238;
	t258 = t233 - t281 + t282;
	t252 = cos(qJ(3));
	t278 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t285 = t278 * t252;
	t254 = -t258 * t250 + t285;
	t290 = qJD(1) * t254;
	t234 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t247);
	t230 = t234 * qJD(4);
	t269 = t250 * qJD(5);
	t289 = (-t233 * t250 + t285) * qJD(3) - t252 * t230 + t269;
	t248 = qJ(1) + pkin(10);
	t242 = cos(t248);
	t246 = qJD(4) + qJD(6);
	t264 = t246 * t252 - qJD(1);
	t287 = t242 * t264;
	t280 = r_i_i_C(2) * t238;
	t262 = r_i_i_C(1) * t237 + t280;
	t284 = t262 * t246 + t230;
	t240 = sin(t248);
	t272 = qJD(1) * t252;
	t263 = -t246 + t272;
	t271 = qJD(3) * t250;
	t283 = -t240 * t271 + t263 * t242;
	t279 = pkin(7) + t234;
	t276 = t246 * t250;
	t255 = t263 * t240 + t242 * t271;
	t226 = t255 * t237 - t238 * t287;
	t227 = t237 * t287 + t255 * t238;
	t275 = t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
	t260 = t264 * t240;
	t228 = t283 * t237 + t238 * t260;
	t229 = t237 * t260 - t283 * t238;
	t274 = -t228 * r_i_i_C(1) + t229 * r_i_i_C(2);
	t273 = qJD(1) * t250;
	t270 = qJD(3) * t252;
	t267 = t278 * t250;
	t261 = t234 * t272 - t230;
	t257 = -t233 * t252 - pkin(2) - t267;
	t231 = t235 * qJD(4);
	t256 = qJD(1) * t235 - t231 * t252 + t234 * t271;
	t253 = qJD(5) * t252 + t284 * t250 + (-t258 * t252 - t267) * qJD(3);
	t232 = t276 * t281;
	t1 = [t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t242 * t231 - t289 * t240 + (-cos(qJ(1)) * pkin(1) - t279 * t240 + t257 * t242) * qJD(1), 0, -t240 * t290 + t253 * t242, t261 * t240 + t256 * t242 + t275, -t240 * t273 + t242 * t270, t275; -t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t240 * t231 + t289 * t242 + (-sin(qJ(1)) * pkin(1) + t279 * t242 + t257 * t240) * qJD(1), 0, t253 * t240 + t242 * t290, t256 * t240 - t261 * t242 + t274, t240 * t270 + t242 * t273, t274; 0, 0, t254 * qJD(3) - t284 * t252 + t269, t232 + (-t246 * t282 - t231) * t250 + (-t234 - t262) * t270, t271, -t270 * t280 + t232 + (-t237 * t270 - t238 * t276) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end