% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP3
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
% Datum: 2019-10-10 01:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
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
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
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
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
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
	% StartTime: 2019-10-10 01:13:12
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 0.27s
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
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, r_i_i_C(1) * t199 + r_i_i_C(2) * t200, 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -r_i_i_C(1) * t201 + r_i_i_C(2) * t202, 0, 0; 0, 0, t211 * qJD(3) - t216 * t223, (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:12
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
	t248 = sin(qJ(3));
	t250 = cos(qJ(3));
	t279 = pkin(8) + r_i_i_C(1);
	t282 = t279 * t250;
	t285 = (-pkin(3) * t248 + t282) * qJD(3);
	t247 = sin(qJ(4));
	t249 = cos(qJ(4));
	t277 = r_i_i_C(3) + qJ(5);
	t280 = pkin(4) - r_i_i_C(2);
	t281 = t247 * t280 - t277 * t249;
	t284 = t281 * qJD(4) - qJD(5) * t247;
	t256 = -t277 * t247 - t249 * t280;
	t253 = -pkin(3) + t256;
	t252 = t248 * t253 + t282;
	t246 = qJ(1) + pkin(9);
	t245 = cos(t246);
	t276 = t245 * t247;
	t275 = t247 * t250;
	t274 = t249 * t250;
	t244 = sin(t246);
	t273 = qJD(1) * t244;
	t272 = qJD(1) * t245;
	t271 = qJD(3) * t248;
	t270 = qJD(3) * t250;
	t269 = qJD(4) * t247;
	t268 = qJD(4) * t249;
	t266 = t279 * t248;
	t263 = t249 * t273;
	t262 = t244 * t271;
	t261 = t244 * t269;
	t260 = t245 * t268;
	t259 = t244 * t247 + t245 * t274;
	t258 = t244 * t275 + t245 * t249;
	t257 = -pkin(3) * t250 - pkin(2) - t266;
	t254 = t244 * t268 + t247 * t272;
	t251 = t284 * t248 + (t250 * t253 - t266) * qJD(3);
	t233 = qJD(1) * t259 - t249 * t262 - t250 * t261 - t260;
	t232 = -t245 * t269 - t247 * t262 + t250 * t254 - t263;
	t231 = t250 * t263 + (t249 * t271 + t250 * t269) * t245 - t254;
	t230 = qJD(1) * t258 - t250 * t260 + t271 * t276 - t261;
	t1 = [-t258 * qJD(5) - t280 * t233 - t277 * t232 - t244 * t285 + (-cos(qJ(1)) * pkin(1) - t244 * pkin(7) + t257 * t245) * qJD(1), 0, t251 * t245 - t252 * t273, t259 * qJD(5) + t230 * t280 - t277 * t231, -t230, 0; -(t244 * t249 - t245 * t275) * qJD(5) - t280 * t231 - t277 * t230 + t245 * t285 + (-sin(qJ(1)) * pkin(1) + t245 * pkin(7) + t257 * t244) * qJD(1), 0, t244 * t251 + t252 * t272, -(-t244 * t274 + t276) * qJD(5) + t277 * t233 - t280 * t232, t232, 0; 0, 0, t252 * qJD(3) - t284 * t250, -t281 * t270 + (qJD(4) * t256 + t249 * qJD(5)) * t248, t247 * t270 + t248 * t268, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:12
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (432->62), mult. (688->93), div. (0->0), fcn. (584->8), ass. (0->43)
	t250 = sin(qJ(3));
	t252 = cos(qJ(3));
	t270 = pkin(5) + pkin(8) + r_i_i_C(1);
	t283 = t270 * t252;
	t287 = (-pkin(3) * t250 + t283) * qJD(3);
	t249 = sin(qJ(4));
	t251 = cos(qJ(4));
	t269 = pkin(4) + r_i_i_C(3) + qJ(6);
	t279 = r_i_i_C(2) + qJ(5);
	t281 = t269 * t249 - t279 * t251;
	t286 = -t270 * qJD(3) + t281 * qJD(4) - qJD(5) * t249 - qJD(6) * t251;
	t257 = -t279 * t249 - t269 * t251;
	t255 = -pkin(3) + t257;
	t284 = t255 * t250 + t283;
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
	t268 = t251 * t276;
	t267 = t249 * t274;
	t266 = t251 * t274;
	t265 = t246 * t272;
	t264 = t247 * t271;
	t260 = t246 * t249 + t247 * t277;
	t232 = t246 * t278 + t247 * t251;
	t259 = -pkin(3) * t252 - t270 * t250 - pkin(2);
	t258 = t246 * t271 + t249 * t275;
	t254 = qJD(3) * t255;
	t253 = t286 * t250 + t252 * t254;
	t234 = -t246 * t251 + t247 * t278;
	t233 = t246 * t277 - t247 * t249;
	t231 = t260 * qJD(1) - t246 * t266 - t252 * t265 - t264;
	t230 = -t246 * t267 - t247 * t272 + t258 * t252 - t268;
	t229 = t252 * t268 + (t252 * t272 + t266) * t247 - t258;
	t228 = t232 * qJD(1) + t247 * t267 - t252 * t264 - t265;
	t1 = [-t232 * qJD(5) - t233 * qJD(6) - t279 * t230 - t269 * t231 - t246 * t287 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t246 + t259 * t247) * qJD(1), 0, t253 * t247 - t284 * t276, qJD(5) * t260 - t234 * qJD(6) + t269 * t228 - t279 * t229, -t228, -t229; t234 * qJD(5) + t260 * qJD(6) - t279 * t228 - t269 * t229 + t247 * t287 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t247 + t259 * t246) * qJD(1), 0, t253 * t246 + t284 * t275, t233 * qJD(5) - t232 * qJD(6) - t269 * t230 + t279 * t231, t230, t231; 0, 0, t250 * t254 - t286 * t252, -t281 * t273 + (t257 * qJD(4) + qJD(5) * t251 - qJD(6) * t249) * t250, t249 * t273 + t250 * t271, -t250 * t272 + t251 * t273;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end