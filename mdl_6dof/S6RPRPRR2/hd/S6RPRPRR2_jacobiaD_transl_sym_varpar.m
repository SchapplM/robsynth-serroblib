% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
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
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
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
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
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
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->19), mult. (94->27), div. (0->0), fcn. (61->8), ass. (0->15)
	t33 = qJ(3) + pkin(11);
	t29 = sin(t33);
	t31 = cos(t33);
	t47 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t29 - cos(qJ(3)) * pkin(3);
	t45 = r_i_i_C(3) + qJ(4) + pkin(7);
	t34 = qJ(1) + pkin(10);
	t32 = cos(t34);
	t44 = qJD(1) * t32;
	t42 = -pkin(2) + t47;
	t41 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t29 + r_i_i_C(2) * t31;
	t30 = sin(t34);
	t40 = t41 * t30;
	t39 = qJD(3) * t47;
	t38 = t41 * qJD(3);
	t1 = [t32 * qJD(4) + qJD(3) * t40 + (-cos(qJ(1)) * pkin(1) - t45 * t30 + t42 * t32) * qJD(1), 0, qJD(1) * t40 + t32 * t39, t44, 0, 0; t30 * qJD(4) - t32 * t38 + (-sin(qJ(1)) * pkin(1) + t45 * t32 + t42 * t30) * qJD(1), 0, t30 * t39 - t41 * t44, qJD(1) * t30, 0, 0; 0, 0, -t38, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:59
	% EndTime: 2019-10-10 00:48:00
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (259->43), mult. (300->68), div. (0->0), fcn. (232->10), ass. (0->35)
	t223 = qJ(3) + pkin(11);
	t219 = sin(t223);
	t221 = cos(t223);
	t252 = pkin(8) + r_i_i_C(3);
	t235 = t252 * t221 - sin(qJ(3)) * pkin(3);
	t258 = (-pkin(4) * t219 + t235) * qJD(3);
	t226 = sin(qJ(5));
	t228 = cos(qJ(5));
	t237 = r_i_i_C(1) * t228 - r_i_i_C(2) * t226 + pkin(4);
	t231 = -t237 * t219 + t235;
	t239 = qJD(1) * t221 - qJD(5);
	t256 = t228 * t239;
	t254 = -t252 * t219 - cos(qJ(3)) * pkin(3);
	t243 = qJD(5) * t221;
	t240 = -qJD(1) + t243;
	t246 = qJD(3) * t226;
	t253 = -t219 * t246 + t240 * t228;
	t224 = qJ(1) + pkin(10);
	t220 = sin(t224);
	t248 = qJD(1) * t220;
	t222 = cos(t224);
	t247 = qJD(1) * t222;
	t245 = qJD(3) * t228;
	t244 = qJD(5) * t219;
	t238 = r_i_i_C(1) * t226 + r_i_i_C(2) * t228;
	t236 = t239 * t226;
	t233 = -pkin(4) * t221 - pkin(2) + t254;
	t232 = t219 * t245 + t240 * t226;
	t230 = t238 * t244 + (-t237 * t221 + t254) * qJD(3);
	t225 = -qJ(4) - pkin(7);
	t217 = t232 * t220 - t222 * t256;
	t216 = t253 * t220 + t222 * t236;
	t215 = t220 * t256 + t232 * t222;
	t214 = t220 * t236 - t253 * t222;
	t1 = [t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t222 * qJD(4) - t220 * t258 + (-cos(qJ(1)) * pkin(1) + t220 * t225 + t233 * t222) * qJD(1), 0, t230 * t222 - t231 * t248, t247, t214 * r_i_i_C(1) + t215 * r_i_i_C(2), 0; -t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t220 * qJD(4) + t222 * t258 + (-sin(qJ(1)) * pkin(1) - t222 * t225 + t233 * t220) * qJD(1), 0, t230 * t220 + t231 * t247, t248, -t216 * r_i_i_C(1) + t217 * r_i_i_C(2), 0; 0, 0, t231 * qJD(3) - t238 * t243, 0, (-t221 * t245 + t226 * t244) * r_i_i_C(2) + (-t221 * t246 - t228 * t244) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:59
	% EndTime: 2019-10-10 00:48:00
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (504->60), mult. (450->90), div. (0->0), fcn. (353->12), ass. (0->52)
	t267 = cos(qJ(5));
	t252 = pkin(5) * t267 + pkin(4);
	t261 = qJ(3) + pkin(11);
	t254 = sin(t261);
	t256 = cos(t261);
	t300 = r_i_i_C(3) + pkin(9) + pkin(8);
	t276 = t300 * t256 - sin(qJ(3)) * pkin(3);
	t265 = sin(qJ(5));
	t299 = pkin(5) * qJD(5);
	t289 = t265 * t299;
	t312 = (-t252 * t254 + t276) * qJD(3) - t256 * t289;
	t263 = qJ(5) + qJ(6);
	t259 = cos(t263);
	t258 = sin(t263);
	t302 = r_i_i_C(2) * t258;
	t277 = r_i_i_C(1) * t259 + t252 - t302;
	t271 = -t277 * t254 + t276;
	t260 = qJD(5) + qJD(6);
	t293 = qJD(1) * t256;
	t283 = -t260 + t293;
	t310 = t259 * t283;
	t309 = t265 * (-qJD(5) + t293);
	t307 = -t300 * t254 - cos(qJ(3)) * pkin(3);
	t284 = t256 * t260 - qJD(1);
	t291 = qJD(3) * t258;
	t306 = -t254 * t291 + t259 * t284;
	t279 = r_i_i_C(1) * t258 + r_i_i_C(2) * t259;
	t305 = t260 * t279 + t289;
	t303 = pkin(5) * t265;
	t297 = t259 * t260;
	t262 = qJ(1) + pkin(10);
	t255 = sin(t262);
	t257 = cos(t262);
	t278 = t283 * t258;
	t247 = t255 * t278 - t257 * t306;
	t290 = qJD(3) * t259;
	t273 = t254 * t290 + t258 * t284;
	t248 = t255 * t310 + t257 * t273;
	t296 = r_i_i_C(1) * t247 + r_i_i_C(2) * t248;
	t249 = t255 * t306 + t257 * t278;
	t250 = t255 * t273 - t257 * t310;
	t295 = -r_i_i_C(1) * t249 + r_i_i_C(2) * t250;
	t294 = qJD(1) * t255;
	t292 = qJD(1) * t257;
	t288 = t267 * t299;
	t285 = qJ(4) + pkin(7) + t303;
	t280 = qJD(4) + t288;
	t274 = -t252 * t256 - pkin(2) + t307;
	t272 = qJD(3) * t254 * t265 + (-qJD(5) * t256 + qJD(1)) * t267;
	t270 = t305 * t254 + (-t256 * t277 + t307) * qJD(3);
	t251 = t254 * t260 * t302;
	t1 = [t250 * r_i_i_C(1) + t249 * r_i_i_C(2) + t280 * t257 - t312 * t255 + (-cos(qJ(1)) * pkin(1) - t285 * t255 + t274 * t257) * qJD(1), 0, t270 * t257 - t271 * t294, t292, (t255 * t309 + t257 * t272) * pkin(5) + t296, t296; -t248 * r_i_i_C(1) + t247 * r_i_i_C(2) + t280 * t255 + t312 * t257 + (-sin(qJ(1)) * pkin(1) + t285 * t257 + t274 * t255) * qJD(1), 0, t255 * t270 + t271 * t292, t294, (t255 * t272 - t257 * t309) * pkin(5) + t295, t295; 0, 0, t271 * qJD(3) - t256 * t305, 0, t251 + (-r_i_i_C(1) * t297 - t288) * t254 + (-t279 - t303) * t256 * qJD(3), -t256 * r_i_i_C(2) * t290 + t251 + (-t254 * t297 - t256 * t291) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end