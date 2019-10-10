% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR1
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
% Datum: 2019-10-10 00:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
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
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (173->33), mult. (136->42), div. (0->0), fcn. (88->10), ass. (0->28)
	t55 = qJ(3) + pkin(11);
	t45 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t55);
	t54 = qJD(3) + qJD(5);
	t52 = qJ(5) + t55;
	t47 = cos(t52);
	t72 = r_i_i_C(2) * t47;
	t46 = sin(t52);
	t74 = r_i_i_C(1) * t46;
	t61 = t72 + t74;
	t59 = t61 * t54;
	t75 = t45 * qJD(3) - t59;
	t73 = r_i_i_C(2) * t46;
	t71 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
	t70 = t47 * t54;
	t56 = qJ(1) + pkin(10);
	t49 = sin(t56);
	t69 = qJD(1) * t49;
	t51 = cos(t56);
	t68 = qJD(1) * t51;
	t67 = r_i_i_C(1) * t70;
	t66 = t54 * t73;
	t64 = qJD(1) * t72;
	t65 = t49 * t64 + t51 * t66 + t69 * t74;
	t62 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t55);
	t63 = t62 * qJD(3) - t67;
	t60 = -r_i_i_C(1) * t47 - pkin(2) + t62 + t73;
	t38 = t49 * t66;
	t1 = [t51 * qJD(4) - t75 * t49 + (-cos(qJ(1)) * pkin(1) - t71 * t49 + t60 * t51) * qJD(1), 0, -t45 * t69 + t63 * t51 + t65, t68, -t51 * t67 + t65, 0; t49 * qJD(4) + t75 * t51 + (-sin(qJ(1)) * pkin(1) + t71 * t51 + t60 * t49) * qJD(1), 0, t38 + t63 * t49 + (t45 - t61) * t68, t69, -t51 * t64 + t38 + (-t46 * t68 - t49 * t70) * r_i_i_C(1), 0; 0, 0, t75, 0, -t59, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:15
	% EndTime: 2019-10-10 00:46:16
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (510->64), mult. (416->95), div. (0->0), fcn. (312->12), ass. (0->52)
	t277 = qJ(3) + pkin(11);
	t274 = qJ(5) + t277;
	t268 = sin(t274);
	t281 = cos(qJ(6));
	t325 = r_i_i_C(1) * t281 + pkin(5);
	t330 = t268 * t325;
	t269 = cos(t274);
	t309 = qJD(6) * t281;
	t276 = qJD(3) + qJD(5);
	t279 = sin(qJ(6));
	t314 = t276 * t279;
	t329 = t268 * t309 + t269 * t314;
	t319 = pkin(9) + r_i_i_C(3);
	t327 = t319 * t269;
	t328 = (-pkin(5) * t268 + t327) * t276;
	t265 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t277);
	t258 = t265 * qJD(3);
	t326 = t258 + t328;
	t310 = qJD(6) * t279;
	t298 = t268 * t310;
	t322 = r_i_i_C(1) * t298 + r_i_i_C(2) * t329;
	t292 = qJD(1) * t269 - qJD(6);
	t321 = t281 * t292;
	t293 = qJD(6) * t269 - qJD(1);
	t304 = t268 * t314;
	t320 = t293 * t281 - t304;
	t316 = r_i_i_C(2) * t279;
	t313 = t276 * t281;
	t278 = qJ(1) + pkin(10);
	t271 = sin(t278);
	t312 = qJD(1) * t271;
	t273 = cos(t278);
	t311 = qJD(1) * t273;
	t308 = qJD(1) * t316;
	t307 = t319 * t268;
	t305 = t319 * t276;
	t303 = t268 * t313;
	t291 = t325 * t276;
	t290 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t277);
	t289 = t322 * t273 + t312 * t330;
	t288 = t292 * t279;
	t287 = t273 * t268 * t308 + t322 * t271 + t311 * t327;
	t286 = -pkin(5) * t269 - pkin(2) + t290 - t307;
	t285 = t293 * t279 + t303;
	t284 = t290 * qJD(3) + (-t269 * t325 - t307) * t276;
	t283 = (-t269 * t310 - t303) * r_i_i_C(1) + (-t269 * t309 + t304) * r_i_i_C(2) + t328;
	t275 = -pkin(8) - qJ(4) - pkin(7);
	t249 = t285 * t271 - t273 * t321;
	t248 = t320 * t271 + t273 * t288;
	t247 = t271 * t321 + t285 * t273;
	t246 = t271 * t288 - t320 * t273;
	t1 = [t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t273 * qJD(4) - t326 * t271 + (-cos(qJ(1)) * pkin(1) + t271 * t275 + t286 * t273) * qJD(1), 0, (-t268 * t316 - t265 - t327) * t312 + t284 * t273 + t289, t311, (-t271 * t308 - t273 * t305) * t268 + (-t273 * t291 - t319 * t312) * t269 + t289, t246 * r_i_i_C(1) + t247 * r_i_i_C(2); -t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t271 * qJD(4) + t326 * t273 + (-sin(qJ(1)) * pkin(1) - t273 * t275 + t286 * t271) * qJD(1), 0, (t265 - t330) * t311 + t284 * t271 + t287, t312, -t271 * t269 * t291 + (-t271 * t305 - t311 * t325) * t268 + t287, -t248 * r_i_i_C(1) + t249 * r_i_i_C(2); 0, 0, t258 + t283, 0, t283, (-t269 * t313 + t298) * r_i_i_C(2) - t329 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end