% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
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
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) - r_i_i_C(1);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->11), mult. (36->20), div. (0->0), fcn. (26->4), ass. (0->8)
	t58 = -pkin(1) - pkin(2);
	t57 = cos(qJ(1));
	t56 = sin(qJ(1));
	t55 = cos(pkin(9));
	t54 = sin(pkin(9));
	t53 = (t54 * t57 - t55 * t56) * qJD(1);
	t52 = (t54 * t56 + t55 * t57) * qJD(1);
	t1 = [-t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t57 * qJD(2) + (-qJ(2) * t56 + t58 * t57) * qJD(1), qJD(1) * t57, 0, 0, 0, 0; t53 * r_i_i_C(1) + t52 * r_i_i_C(2) + t56 * qJD(2) + (qJ(2) * t57 + t58 * t56) * qJD(1), qJD(1) * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->23), mult. (138->41), div. (0->0), fcn. (116->6), ass. (0->17)
	t47 = sin(qJ(4));
	t49 = cos(qJ(4));
	t51 = (r_i_i_C(1) * t47 + r_i_i_C(2) * t49) * qJD(4);
	t57 = -pkin(1) - pkin(2);
	t56 = -pkin(7) - r_i_i_C(3);
	t55 = qJD(4) * t47;
	t54 = qJD(4) * t49;
	t45 = sin(pkin(9));
	t46 = cos(pkin(9));
	t48 = sin(qJ(1));
	t50 = cos(qJ(1));
	t42 = t45 * t50 - t48 * t46;
	t43 = t48 * t45 + t46 * t50;
	t52 = r_i_i_C(1) * t49 - r_i_i_C(2) * t47 + pkin(3);
	t41 = t42 * qJD(1);
	t40 = t43 * qJD(1);
	t1 = [t50 * qJD(2) + t56 * t41 - t42 * t51 - t52 * t40 + (-qJ(2) * t48 + t57 * t50) * qJD(1), qJD(1) * t50, 0, (-t41 * t49 + t43 * t55) * r_i_i_C(2) + (-t41 * t47 - t43 * t54) * r_i_i_C(1), 0, 0; t48 * qJD(2) + t56 * t40 - t43 * t51 + t52 * t41 + (qJ(2) * t50 + t57 * t48) * qJD(1), qJD(1) * t48, 0, (-t40 * t49 - t42 * t55) * r_i_i_C(2) + (-t40 * t47 + t42 * t54) * r_i_i_C(1), 0, 0; 0, 0, 0, t51, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:57
	% EndTime: 2019-10-09 23:52:58
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (171->48), mult. (482->79), div. (0->0), fcn. (461->8), ass. (0->37)
	t235 = sin(qJ(4));
	t234 = sin(qJ(5));
	t236 = cos(qJ(5));
	t243 = t236 * r_i_i_C(1) - t234 * r_i_i_C(2) + pkin(4);
	t237 = cos(qJ(4));
	t260 = pkin(8) + r_i_i_C(3);
	t250 = t260 * t237;
	t239 = -t243 * t235 + t250;
	t267 = t239 * qJD(4);
	t249 = t260 * t235;
	t266 = t243 * t237 + t249;
	t256 = sin(pkin(9));
	t257 = cos(pkin(9));
	t258 = sin(qJ(1));
	t259 = cos(qJ(1));
	t229 = t259 * t256 - t258 * t257;
	t263 = qJD(1) * t258;
	t262 = qJD(1) * t259;
	t261 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t255 = qJD(4) * t235;
	t254 = qJD(4) * t237;
	t253 = qJD(5) * t234;
	t252 = qJD(5) * t236;
	t251 = qJD(5) * t237;
	t246 = t234 * r_i_i_C(1) + t236 * r_i_i_C(2);
	t228 = -t258 * t256 - t259 * t257;
	t226 = t228 * qJD(1);
	t245 = t228 * t251 + t226;
	t227 = t229 * qJD(1);
	t244 = t229 * t251 + t227;
	t242 = qJD(5) * t246;
	t241 = qJD(5) * t228 + t226 * t237 - t229 * t255;
	t240 = qJD(5) * t229 + t227 * t237 + t228 * t255;
	t238 = t266 * qJD(4) - t235 * t242;
	t225 = t245 * t234 + t240 * t236;
	t224 = -t240 * t234 + t245 * t236;
	t1 = [(-t227 * t234 + t228 * t252) * r_i_i_C(1) + (-t227 * t236 - t228 * t253) * r_i_i_C(2) - t227 * pkin(7) - qJ(2) * t263 - (-pkin(3) - t266) * t226 + (-t237 * t242 + t267) * t229 + t261 * t259, t262, 0, t239 * t227 + t238 * t228, t224 * r_i_i_C(1) - t225 * r_i_i_C(2), 0; t226 * pkin(7) + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) + (pkin(4) * t235 - t250) * t228 * qJD(4) + qJ(2) * t262 + (pkin(4) * t237 + pkin(3) + t249) * t227 + t261 * t258, t263, 0, -t239 * t226 + t238 * t229, (t244 * r_i_i_C(1) + t241 * r_i_i_C(2)) * t236 + (t241 * r_i_i_C(1) - t244 * r_i_i_C(2)) * t234, 0; 0, 0, 0, t246 * t251 - t267, (-t235 * t253 + t236 * t254) * r_i_i_C(2) + (t234 * t254 + t235 * t252) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:58
	% EndTime: 2019-10-09 23:52:59
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (324->60), mult. (890->94), div. (0->0), fcn. (893->8), ass. (0->44)
	t292 = sin(qJ(4));
	t294 = cos(qJ(4));
	t325 = pkin(8) + r_i_i_C(2);
	t331 = t325 * t294;
	t335 = (-pkin(4) * t292 + t331) * qJD(4);
	t310 = t325 * t292;
	t334 = t294 * pkin(4) + pkin(3) + t310;
	t291 = sin(qJ(5));
	t293 = cos(qJ(5));
	t320 = r_i_i_C(3) + qJ(6);
	t324 = r_i_i_C(1) + pkin(5);
	t300 = t320 * t291 + t324 * t293;
	t298 = pkin(4) + t300;
	t332 = t298 * t292 - t331;
	t318 = sin(pkin(9));
	t319 = cos(pkin(9));
	t322 = sin(qJ(1));
	t323 = cos(qJ(1));
	t286 = t323 * t318 - t322 * t319;
	t330 = qJD(1) * t322;
	t329 = qJD(1) * t323;
	t299 = t324 * t291 - t320 * t293;
	t328 = t299 * qJD(5) - qJD(6) * t291;
	t327 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
	t285 = -t322 * t318 - t323 * t319;
	t283 = t285 * qJD(1);
	t315 = qJD(4) * t292;
	t306 = t286 * t315;
	t297 = -qJD(5) * t285 - t283 * t294 + t306;
	t284 = t286 * qJD(1);
	t313 = qJD(5) * t294;
	t303 = t286 * t313 + t284;
	t326 = t303 * t291 + t297 * t293;
	t317 = t291 * t294;
	t316 = t293 * t294;
	t314 = qJD(4) * t294;
	t307 = t285 * t315;
	t302 = t285 * t291 + t286 * t316;
	t301 = t285 * t316 - t286 * t291;
	t295 = -t328 * t292 + (t298 * t294 + t310) * qJD(4);
	t278 = (t285 * t313 + t283) * t291 + (qJD(5) * t286 + t284 * t294 + t307) * t293;
	t277 = -t301 * qJD(5) - t283 * t293 + t284 * t317 + t291 * t307;
	t273 = -t302 * qJD(5) - t283 * t317 - t284 * t293 + t291 * t306;
	t1 = [-(t285 * t293 - t286 * t317) * qJD(6) - t284 * pkin(7) - t324 * t326 + t320 * (-t297 * t291 + t303 * t293) + t286 * t335 - qJ(2) * t330 + t334 * t283 + t327 * t323, t329, 0, -t284 * t332 + t295 * t285, -t301 * qJD(6) - t324 * t277 + t320 * t278, t277; -(t285 * t317 + t286 * t293) * qJD(6) + t283 * pkin(7) + t324 * t278 + t320 * t277 - t285 * t335 + qJ(2) * t329 + t334 * t284 + t327 * t322, t330, 0, t283 * t332 + t295 * t286, -t302 * qJD(6) - t324 * t273 + t320 * t326, t273; 0, 0, 0, t332 * qJD(4) + t328 * t294, t299 * t314 + (t300 * qJD(5) - qJD(6) * t293) * t292, -t292 * qJD(5) * t293 - t291 * t314;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end