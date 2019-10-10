% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (149->33), mult. (146->42), div. (0->0), fcn. (95->8), ass. (0->32)
	t52 = sin(qJ(2));
	t50 = qJD(2) + qJD(3);
	t51 = qJ(2) + qJ(3);
	t46 = pkin(10) + t51;
	t44 = sin(t46);
	t45 = cos(t46);
	t59 = r_i_i_C(1) * t44 + r_i_i_C(2) * t45;
	t47 = sin(t51);
	t75 = pkin(3) * t47;
	t78 = t59 + t75;
	t56 = t78 * t50;
	t67 = pkin(2) * qJD(2);
	t76 = t52 * t67 + t56;
	t48 = cos(t51);
	t72 = r_i_i_C(1) * t45;
	t77 = -pkin(3) * t48 - t72;
	t71 = r_i_i_C(2) * t44;
	t69 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t68 = t48 * t50;
	t53 = sin(qJ(1));
	t66 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t65 = qJD(1) * t55;
	t64 = t50 * t72;
	t63 = t50 * t71;
	t62 = t55 * t63 + t59 * t66;
	t54 = cos(qJ(2));
	t60 = -pkin(3) * t68 - t54 * t67 - t64;
	t58 = -t54 * pkin(2) - pkin(1) + t71 + t77;
	t43 = -t52 * pkin(2) - t75;
	t38 = t53 * t63;
	t1 = [t55 * qJD(4) + t76 * t53 + (-t69 * t53 + t58 * t55) * qJD(1), -t43 * t66 + t60 * t55 + t62, -t55 * t64 + (t47 * t66 - t55 * t68) * pkin(3) + t62, t65, 0, 0; t53 * qJD(4) - t76 * t55 + (t58 * t53 + t69 * t55) * qJD(1), t38 + t60 * t53 + (t43 - t59) * t65, t50 * t53 * t77 - t65 * t78 + t38, t66, 0, 0; 0, -t76, -t56, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:05
	% EndTime: 2019-10-10 11:17:06
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (340->50), mult. (317->57), div. (0->0), fcn. (230->10), ass. (0->45)
	t236 = qJ(2) + qJ(3);
	t231 = pkin(10) + t236;
	t230 = cos(t231);
	t270 = r_i_i_C(3) + qJ(5);
	t257 = t270 * t230;
	t229 = sin(t231);
	t228 = t229 * qJD(5);
	t232 = sin(t236);
	t235 = qJD(2) + qJD(3);
	t237 = sin(pkin(11));
	t238 = cos(pkin(11));
	t271 = r_i_i_C(2) * t237;
	t280 = r_i_i_C(1) * t238 + pkin(4);
	t251 = t280 - t271;
	t239 = sin(qJ(2));
	t269 = pkin(2) * qJD(2);
	t263 = t239 * t269;
	t273 = pkin(3) * t235;
	t281 = (-t251 * t229 + t257) * t235 + (t237 * r_i_i_C(1) + t238 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(4)) * qJD(1) - t232 * t273 + t228 - t263;
	t264 = t235 * t271;
	t279 = (qJD(5) + t264) * t230;
	t275 = pkin(3) * t232;
	t233 = cos(t236);
	t274 = pkin(3) * t233;
	t268 = t230 * t235;
	t240 = sin(qJ(1));
	t267 = qJD(1) * t240;
	t242 = cos(qJ(1));
	t266 = qJD(1) * t242;
	t261 = t229 * t267;
	t260 = t229 * t266;
	t258 = t270 * t229;
	t254 = t279 * t242 + t280 * t261;
	t253 = t280 * t229;
	t252 = t279 * t240 + t266 * t257 + t260 * t271;
	t249 = -t229 * t271 - t257;
	t248 = -t253 - t275;
	t247 = -t230 * t280 - t258;
	t241 = cos(qJ(2));
	t246 = -t233 * t273 + t247 * t235 - t241 * t269;
	t245 = t235 * (t247 - t274);
	t244 = t229 * t264 + t248 * t235 + t270 * t268 + t228;
	t243 = qJD(4) + (-t241 * pkin(2) - t251 * t230 - pkin(1) - t258 - t274) * qJD(1);
	t225 = -t239 * pkin(2) - t275;
	t1 = [-t281 * t240 + t243 * t242, (-t225 + t249) * t267 + t246 * t242 + t254, (t249 + t275) * t267 + t242 * t245 + t254, t266, t242 * t268 - t261, 0; t243 * t240 + t281 * t242, (t225 - t253) * t266 + t246 * t240 + t252, t240 * t245 + t248 * t266 + t252, t267, t240 * t268 + t260, 0; 0, t244 - t263, t244, 0, t235 * t229, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:05
	% EndTime: 2019-10-10 11:17:06
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (557->75), mult. (469->99), div. (0->0), fcn. (360->12), ass. (0->65)
	t279 = pkin(11) + qJ(6);
	t273 = sin(t279);
	t274 = cos(t279);
	t281 = qJ(2) + qJ(3);
	t275 = pkin(10) + t281;
	t270 = sin(t275);
	t318 = qJD(6) * t270;
	t271 = cos(t275);
	t280 = qJD(2) + qJD(3);
	t325 = t271 * t280;
	t342 = t273 * t325 + t274 * t318;
	t272 = cos(pkin(11)) * pkin(5) + pkin(4);
	t341 = r_i_i_C(1) * t274 + t272;
	t269 = t270 * qJD(5);
	t276 = sin(t281);
	t284 = sin(qJ(2));
	t327 = pkin(2) * qJD(2);
	t315 = t284 * t327;
	t283 = -pkin(9) - qJ(5);
	t328 = r_i_i_C(3) - t283;
	t332 = pkin(3) * t280;
	t340 = (-t270 * t272 + t328 * t271) * t280 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(4)) * qJD(1) - t276 * t332 + t269 - t315;
	t309 = t273 * t318;
	t339 = r_i_i_C(1) * t309 + t342 * r_i_i_C(2) + qJD(5) * t271;
	t287 = cos(qJ(1));
	t317 = qJD(6) * t271;
	t301 = -qJD(1) + t317;
	t338 = t287 * t301;
	t300 = qJD(1) * t271 - qJD(6);
	t285 = sin(qJ(1));
	t323 = t280 * t285;
	t314 = t270 * t323;
	t336 = t300 * t287 - t314;
	t334 = pkin(3) * t276;
	t277 = cos(t281);
	t333 = pkin(3) * t277;
	t330 = r_i_i_C(2) * t273;
	t329 = r_i_i_C(3) * t271;
	t324 = t271 * t283;
	t322 = t280 * t287;
	t321 = qJD(1) * t285;
	t320 = qJD(1) * t287;
	t316 = t270 * t330;
	t313 = t270 * t322;
	t311 = t270 * t321;
	t310 = t270 * t320;
	t299 = t301 * t285;
	t298 = -t316 - t329;
	t297 = t283 * t314 + t339 * t285 + t310 * t330 + t320 * t329;
	t296 = -r_i_i_C(3) * t270 - t271 * t341;
	t295 = -t270 * t341 - t324;
	t294 = t283 * t313 + t339 * t287 + t341 * t311 + t321 * t324;
	t293 = t300 * t285 + t313;
	t292 = t295 - t334;
	t286 = cos(qJ(2));
	t291 = qJD(4) + (-t286 * pkin(2) - t328 * t270 - t271 * t272 - pkin(1) - t333) * qJD(1);
	t290 = -t277 * t332 + t296 * t280 - t286 * t327;
	t289 = t280 * (t296 - t333);
	t288 = r_i_i_C(3) * t325 + t269 + (-r_i_i_C(1) * t273 - r_i_i_C(2) * t274) * t317 + (t316 + t292) * t280;
	t266 = -t284 * pkin(2) - t334;
	t247 = t273 * t299 - t336 * t274;
	t246 = t336 * t273 + t274 * t299;
	t245 = t273 * t338 + t293 * t274;
	t244 = t293 * t273 - t274 * t338;
	t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) - t340 * t285 + t291 * t287, (-t266 + t298) * t321 + t290 * t287 + t294, (t298 + t334) * t321 + t287 * t289 + t294, t320, t271 * t322 - t311, t244 * r_i_i_C(1) + t245 * r_i_i_C(2); -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t291 * t285 + t340 * t287, t290 * t285 + (t266 + t295) * t320 + t297, t285 * t289 + t292 * t320 + t297, t321, t271 * t323 + t310, -t246 * r_i_i_C(1) + t247 * r_i_i_C(2); 0, t288 - t315, t288, 0, t280 * t270, (-t274 * t325 + t309) * r_i_i_C(2) - t342 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end