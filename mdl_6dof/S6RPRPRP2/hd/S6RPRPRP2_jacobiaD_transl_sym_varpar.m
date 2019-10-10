% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP2
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
% Datum: 2019-10-10 00:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
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
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
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
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
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
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->19), mult. (94->27), div. (0->0), fcn. (61->8), ass. (0->15)
	t33 = qJ(3) + pkin(10);
	t29 = sin(t33);
	t31 = cos(t33);
	t47 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t29 - cos(qJ(3)) * pkin(3);
	t45 = r_i_i_C(3) + qJ(4) + pkin(7);
	t34 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 00:30:48
	% EndTime: 2019-10-10 00:30:48
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (259->43), mult. (300->68), div. (0->0), fcn. (232->10), ass. (0->35)
	t223 = qJ(3) + pkin(10);
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
	t224 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 00:30:49
	% EndTime: 2019-10-10 00:30:49
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (472->59), mult. (540->85), div. (0->0), fcn. (448->10), ass. (0->44)
	t272 = qJ(3) + pkin(10);
	t268 = sin(t272);
	t270 = cos(t272);
	t310 = pkin(8) + r_i_i_C(2);
	t290 = t310 * t270 - sin(qJ(3)) * pkin(3);
	t317 = (-pkin(4) * t268 + t290) * qJD(3);
	t275 = sin(qJ(5));
	t277 = cos(qJ(5));
	t306 = r_i_i_C(3) + qJ(6);
	t311 = pkin(5) + r_i_i_C(1);
	t285 = -t306 * t275 - t311 * t277;
	t281 = -pkin(4) + t285;
	t280 = t281 * t268 + t290;
	t312 = t311 * t275 - t306 * t277;
	t315 = t312 * qJD(5) - qJD(6) * t275;
	t313 = -t310 * t268 - cos(qJ(3)) * pkin(3);
	t273 = qJ(1) + pkin(9);
	t269 = sin(t273);
	t305 = t269 * t275;
	t304 = t269 * t277;
	t271 = cos(t273);
	t303 = t271 * t275;
	t302 = t271 * t277;
	t301 = qJD(1) * t269;
	t300 = qJD(1) * t271;
	t299 = qJD(3) * t275;
	t298 = qJD(5) * t275;
	t297 = qJD(5) * t277;
	t294 = t268 * t299;
	t293 = qJD(3) * t268 * t277;
	t292 = t269 * t298;
	t291 = t271 * t297;
	t288 = t270 * t302 + t305;
	t287 = t270 * t305 + t302;
	t286 = -pkin(4) * t270 - pkin(2) + t313;
	t283 = t269 * t297 + t275 * t300;
	t282 = t271 * t298 + t277 * t301;
	t279 = t315 * t268 + (t281 * t270 + t313) * qJD(3);
	t274 = -qJ(4) - pkin(7);
	t256 = t288 * qJD(1) - t269 * t293 - t270 * t292 - t291;
	t255 = -t269 * t294 + t283 * t270 - t282;
	t254 = t282 * t270 + t271 * t293 - t283;
	t253 = t287 * qJD(1) - t270 * t291 + t271 * t294 - t292;
	t1 = [-t287 * qJD(6) + t271 * qJD(4) - t311 * t256 - t306 * t255 - t269 * t317 + (-cos(qJ(1)) * pkin(1) + t269 * t274 + t286 * t271) * qJD(1), 0, t279 * t271 - t280 * t301, t300, t288 * qJD(6) + t311 * t253 - t306 * t254, -t253; -(-t270 * t303 + t304) * qJD(6) + t269 * qJD(4) - t311 * t254 - t306 * t253 + t271 * t317 + (-sin(qJ(1)) * pkin(1) - t271 * t274 + t286 * t269) * qJD(1), 0, t279 * t269 + t280 * t300, t301, -(-t270 * t304 + t303) * qJD(6) + t306 * t256 - t311 * t255, t255; 0, 0, t280 * qJD(3) - t315 * t270, 0, -t312 * t270 * qJD(3) + (t285 * qJD(5) + t277 * qJD(6)) * t268, t268 * t297 + t270 * t299;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end