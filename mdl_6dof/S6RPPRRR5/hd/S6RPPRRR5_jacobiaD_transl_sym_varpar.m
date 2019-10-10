% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
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
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (13->9), mult. (28->12), div. (0->0), fcn. (18->2), ass. (0->7)
	t12 = r_i_i_C(2) + qJ(2);
	t8 = sin(qJ(1));
	t11 = qJD(1) * t8;
	t10 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t9 = cos(qJ(1));
	t7 = qJD(1) * t9;
	t1 = [t9 * qJD(2) - t8 * qJD(3) + (t10 * t9 - t12 * t8) * qJD(1), t7, -t11, 0, 0, 0; t8 * qJD(2) + t9 * qJD(3) + (t10 * t8 + t12 * t9) * qJD(1), t11, t7, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (28->20), mult. (80->31), div. (0->0), fcn. (52->4), ass. (0->13)
	t18 = sin(qJ(4));
	t20 = cos(qJ(4));
	t23 = (r_i_i_C(1) * t20 - r_i_i_C(2) * t18) * qJD(4);
	t29 = qJD(3) + t23;
	t19 = sin(qJ(1));
	t28 = qJD(1) * t19;
	t21 = cos(qJ(1));
	t17 = qJD(1) * t21;
	t27 = qJD(4) * t19;
	t26 = qJD(4) * t21;
	t25 = pkin(7) + r_i_i_C(3) - qJ(2);
	t22 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20 - pkin(1) - qJ(3);
	t1 = [t21 * qJD(2) - t29 * t19 + (t19 * t25 + t21 * t22) * qJD(1), t17, -t28, (t18 * t28 - t20 * t26) * r_i_i_C(2) + (-t18 * t26 - t20 * t28) * r_i_i_C(1), 0, 0; t19 * qJD(2) + t29 * t21 + (t19 * t22 - t21 * t25) * qJD(1), t28, t17, (-t17 * t18 - t20 * t27) * r_i_i_C(2) + (t17 * t20 - t18 * t27) * r_i_i_C(1), 0, 0; 0, 0, 0, -t23, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:19
	% EndTime: 2019-10-10 00:08:19
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (90->30), mult. (134->39), div. (0->0), fcn. (87->6), ass. (0->28)
	t43 = cos(qJ(4));
	t39 = qJD(4) + qJD(5);
	t40 = qJ(4) + qJ(5);
	t37 = sin(t40);
	t59 = r_i_i_C(2) * t37;
	t38 = cos(t40);
	t60 = r_i_i_C(1) * t38;
	t51 = (t59 - t60) * t39;
	t55 = pkin(4) * qJD(4);
	t65 = -t43 * t55 + t51;
	t63 = qJD(3) - t65;
	t62 = pkin(4) * t43;
	t61 = r_i_i_C(1) * t37;
	t58 = r_i_i_C(2) * t38;
	t42 = sin(qJ(1));
	t57 = t39 * t42;
	t44 = cos(qJ(1));
	t56 = t39 * t44;
	t54 = qJD(1) * t42;
	t36 = qJD(1) * t44;
	t53 = r_i_i_C(3) - qJ(2) + pkin(8) + pkin(7);
	t49 = -t58 - t61;
	t41 = sin(qJ(4));
	t47 = -pkin(4) * t41 - pkin(1) - qJ(3) + t49;
	t46 = t49 * t39 - t41 * t55;
	t34 = t36 * t60;
	t33 = t54 * t59;
	t1 = [t44 * qJD(2) - t63 * t42 + (t53 * t42 + t47 * t44) * qJD(1), t36, -t54, t33 + (-t60 - t62) * t54 + t46 * t44, -t56 * t58 + t33 + (-t37 * t56 - t38 * t54) * r_i_i_C(1), 0; t42 * qJD(2) + t63 * t44 + (t47 * t42 - t53 * t44) * qJD(1), t54, t36, t34 + (-t59 + t62) * t36 + t46 * t42, -t57 * t61 + t34 + (-t37 * t36 - t38 * t57) * r_i_i_C(2), 0; 0, 0, 0, t65, t51, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:08:20
	% EndTime: 2019-10-10 00:08:20
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (279->60), mult. (414->86), div. (0->0), fcn. (311->8), ass. (0->50)
	t307 = pkin(9) + r_i_i_C(3);
	t263 = qJ(4) + qJ(5);
	t261 = cos(t263);
	t267 = cos(qJ(6));
	t312 = (r_i_i_C(1) * t267 + pkin(5)) * t261;
	t262 = qJD(4) + qJD(5);
	t260 = sin(t263);
	t292 = t307 * t260;
	t268 = cos(qJ(4));
	t303 = pkin(4) * qJD(4);
	t293 = t268 * t303;
	t311 = (pkin(5) * t261 + t292) * t262 + qJD(3) + t293;
	t306 = pkin(4) * t268;
	t302 = t260 * t262;
	t301 = t261 * t262;
	t300 = t262 * t267;
	t269 = cos(qJ(1));
	t299 = t267 * t269;
	t298 = qJ(2) - pkin(8) - pkin(7);
	t266 = sin(qJ(1));
	t297 = qJD(1) * t266;
	t259 = qJD(1) * t269;
	t264 = sin(qJ(6));
	t296 = qJD(6) * t264;
	t295 = qJD(6) * t267;
	t294 = r_i_i_C(2) * t261 * t264;
	t291 = t264 * t302;
	t290 = t266 * t301;
	t289 = t269 * t301;
	t286 = t261 * t295;
	t285 = r_i_i_C(2) * t291;
	t284 = qJD(1) * t294;
	t283 = qJD(6) * t260 + qJD(1);
	t282 = qJD(1) * t260 + qJD(6);
	t281 = t266 * t284 + t269 * t285 + t307 * t289;
	t279 = t283 * t264;
	t278 = t259 * t312 + t266 * t285 + t307 * (t260 * t259 + t290);
	t277 = t260 * t300 + t261 * t296;
	t276 = -t292 - t312;
	t275 = t282 * t266 - t289;
	t265 = sin(qJ(4));
	t274 = -pkin(4) * t265 - pkin(5) * t260 + t307 * t261 - pkin(1) - qJ(3);
	t273 = (t276 + t294) * t262 + (r_i_i_C(1) * t296 + r_i_i_C(2) * t295) * t260;
	t272 = -pkin(5) * t302 - t277 * r_i_i_C(1) - r_i_i_C(2) * t286;
	t271 = -t265 * t303 + t272;
	t244 = -t282 * t299 + (-t261 * t300 + t279) * t266;
	t243 = t283 * t267 * t266 + (t282 * t269 + t290) * t264;
	t242 = t275 * t267 + t269 * t279;
	t241 = t275 * t264 - t283 * t299;
	t1 = [t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t269 * qJD(2) - t311 * t266 + (-t298 * t266 + t274 * t269) * qJD(1), t259, -t297, t271 * t269 + (t276 - t306) * t297 + t281, t272 * t269 + t276 * t297 + t281, t241 * r_i_i_C(1) + t242 * r_i_i_C(2); -t242 * r_i_i_C(1) + t241 * r_i_i_C(2) + t266 * qJD(2) + t311 * t269 + (t274 * t266 + t298 * t269) * qJD(1), t297, t259, (-t294 + t306) * t259 + t271 * t266 + t278, t272 * t266 - t269 * t284 + t278, -t243 * r_i_i_C(1) + t244 * r_i_i_C(2); 0, 0, 0, t273 - t293, t273, t277 * r_i_i_C(2) + (-t286 + t291) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end