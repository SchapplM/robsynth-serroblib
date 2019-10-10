% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
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
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(11)) + r_i_i_C(2) * sin(pkin(11)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:21
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(11) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(11)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (124->30), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->28)
	t44 = qJD(3) + qJD(4);
	t43 = pkin(11) + qJ(3);
	t41 = qJ(4) + t43;
	t38 = cos(t41);
	t61 = r_i_i_C(2) * t38;
	t37 = sin(t41);
	t63 = r_i_i_C(1) * t37;
	t51 = t61 + t63;
	t49 = t51 * t44;
	t39 = sin(t43);
	t64 = pkin(3) * t39;
	t65 = qJD(3) * t64 + t49;
	t62 = r_i_i_C(2) * t37;
	t60 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
	t59 = t38 * t44;
	t45 = sin(qJ(1));
	t58 = qJD(1) * t45;
	t46 = cos(qJ(1));
	t57 = qJD(1) * t46;
	t40 = cos(t43);
	t56 = qJD(3) * t40;
	t55 = r_i_i_C(1) * t59;
	t54 = t44 * t62;
	t52 = qJD(1) * t61;
	t50 = -r_i_i_C(1) * t38 - pkin(3) * t40 - cos(pkin(11)) * pkin(2) - pkin(1) + t62;
	t48 = t45 * t52 + t58 * t63 + (t54 - t55) * t46;
	t32 = t45 * t54;
	t1 = [t46 * qJD(2) + t65 * t45 + (-t60 * t45 + t50 * t46) * qJD(1), t57, (t39 * t58 - t46 * t56) * pkin(3) + t48, t48, 0, 0; t45 * qJD(2) - t65 * t46 + (t50 * t45 + t60 * t46) * qJD(1), t58, t32 + (-pkin(3) * t56 - t55) * t45 + (-t51 - t64) * t57, -t46 * t52 + t32 + (-t37 * t57 - t45 * t59) * r_i_i_C(1), 0, 0; 0, 0, -t65, -t49, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:21
	% EndTime: 2019-10-10 09:04:22
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (269->40), mult. (176->50), div. (0->0), fcn. (115->9), ass. (0->38)
	t56 = pkin(11) + qJ(3);
	t51 = sin(t56);
	t57 = qJD(3) + qJD(4);
	t53 = qJD(5) + t57;
	t54 = qJ(4) + t56;
	t50 = qJ(5) + t54;
	t47 = cos(t50);
	t79 = r_i_i_C(2) * t47;
	t46 = sin(t50);
	t81 = r_i_i_C(1) * t46;
	t64 = t79 + t81;
	t62 = t64 * t53;
	t48 = sin(t54);
	t82 = pkin(4) * t48;
	t60 = -t57 * t82 - t62;
	t75 = pkin(3) * qJD(3);
	t84 = -t51 * t75 + t60;
	t77 = t47 * t53;
	t70 = r_i_i_C(1) * t77;
	t49 = cos(t54);
	t76 = t49 * t57;
	t83 = -pkin(4) * t76 - t70;
	t80 = r_i_i_C(2) * t46;
	t78 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7) + qJ(2);
	t58 = sin(qJ(1));
	t74 = qJD(1) * t58;
	t59 = cos(qJ(1));
	t73 = qJD(1) * t59;
	t69 = t53 * t80;
	t66 = qJD(1) * t79;
	t68 = t58 * t66 + t59 * t69 + t74 * t81;
	t52 = cos(t56);
	t65 = -t52 * t75 + t83;
	t63 = -r_i_i_C(1) * t47 - pkin(4) * t49 - pkin(3) * t52 - cos(pkin(11)) * pkin(2) - pkin(1) + t80;
	t61 = -t59 * t70 + t68;
	t43 = -pkin(3) * t51 - t82;
	t41 = t58 * t69;
	t1 = [qJD(2) * t59 - t84 * t58 + (-t78 * t58 + t63 * t59) * qJD(1), t73, -t43 * t74 + t65 * t59 + t68, (t48 * t74 - t59 * t76) * pkin(4) + t61, t61, 0; t58 * qJD(2) + t84 * t59 + (t63 * t58 + t78 * t59) * qJD(1), t74, t41 + t65 * t58 + (t43 - t64) * t73, t41 + t83 * t58 + (-t64 - t82) * t73, -t59 * t66 + t41 + (-t46 * t73 - t58 * t77) * r_i_i_C(1), 0; 0, 0, t84, t60, -t62, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:04:23
	% EndTime: 2019-10-10 09:04:23
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (738->75), mult. (530->107), div. (0->0), fcn. (392->11), ass. (0->65)
	t275 = pkin(11) + qJ(3);
	t273 = qJ(4) + t275;
	t269 = qJ(5) + t273;
	t265 = sin(t269);
	t279 = cos(qJ(6));
	t333 = r_i_i_C(1) * t279 + pkin(5);
	t294 = t333 * t265;
	t277 = sin(qJ(6));
	t315 = qJD(6) * t279;
	t266 = cos(t269);
	t276 = qJD(3) + qJD(4);
	t272 = qJD(5) + t276;
	t323 = t266 * t272;
	t335 = t265 * t315 + t277 * t323;
	t330 = pkin(10) + r_i_i_C(3);
	t308 = t330 * t266;
	t270 = sin(t275);
	t324 = pkin(3) * qJD(3);
	t311 = t270 * t324;
	t267 = sin(t273);
	t328 = pkin(4) * t276;
	t314 = t267 * t328;
	t327 = pkin(5) * t265;
	t334 = (t308 - t327) * t272 - t311 - t314;
	t316 = qJD(6) * t277;
	t301 = t265 * t316;
	t331 = r_i_i_C(1) * t301 + t335 * r_i_i_C(2);
	t268 = cos(t273);
	t293 = t333 * t266;
	t309 = t330 * t265;
	t283 = (-t293 - t309) * t272 - t268 * t328;
	t329 = pkin(4) * t267;
	t325 = r_i_i_C(2) * t277;
	t278 = sin(qJ(1));
	t322 = t272 * t278;
	t321 = t272 * t279;
	t280 = cos(qJ(1));
	t320 = t272 * t280;
	t319 = t279 * t280;
	t318 = qJD(1) * t278;
	t317 = qJD(1) * t280;
	t312 = t265 * t325;
	t310 = qJD(1) * t325;
	t307 = t330 * t278;
	t306 = t265 * t321;
	t296 = qJD(6) * t266 - qJD(1);
	t295 = qJD(1) * t266 - qJD(6);
	t292 = t333 * t280;
	t291 = t331 * t280 + t318 * t294;
	t290 = t296 * t277;
	t289 = t280 * t265 * t310 + t331 * t278 + t317 * t308;
	t288 = -t308 - t312;
	t271 = cos(t275);
	t287 = -pkin(5) * t266 - pkin(4) * t268 - pkin(3) * t271 - cos(pkin(11)) * pkin(2) - pkin(1) - t309;
	t286 = t265 * t320 + t295 * t278;
	t284 = -t271 * t324 + t283;
	t282 = -t266 * r_i_i_C(2) * t315 + (-t266 * t316 - t306) * r_i_i_C(1) + t330 * t323 + (-t327 + t312) * t272;
	t281 = t282 - t314;
	t274 = -pkin(9) - pkin(8) - pkin(7) - qJ(2);
	t259 = -pkin(3) * t270 - t329;
	t246 = -t295 * t319 + (t290 + t306) * t278;
	t245 = t296 * t279 * t278 + (-t265 * t322 + t295 * t280) * t277;
	t244 = t286 * t279 + t280 * t290;
	t243 = t286 * t277 - t296 * t319;
	t1 = [t246 * r_i_i_C(1) + t245 * r_i_i_C(2) + t280 * qJD(2) - t334 * t278 + (t274 * t278 + t287 * t280) * qJD(1), t317, (-t259 + t288) * t318 + t284 * t280 + t291, (t288 + t329) * t318 + t283 * t280 + t291, (-t278 * t310 - t330 * t320) * t265 + (-qJD(1) * t307 - t272 * t292) * t266 + t291, t243 * r_i_i_C(1) + t244 * r_i_i_C(2); -t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t278 * qJD(2) + t334 * t280 + (-t274 * t280 + t287 * t278) * qJD(1), t318, (t259 - t294) * t317 + t284 * t278 + t289, (-t294 - t329) * t317 + t283 * t278 + t289, -t293 * t322 + (-qJD(1) * t292 - t272 * t307) * t265 + t289, -t245 * r_i_i_C(1) + t246 * r_i_i_C(2); 0, 0, t281 - t311, t281, t282, (-t266 * t321 + t301) * r_i_i_C(2) - t335 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end