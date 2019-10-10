% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
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
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (23->16), mult. (72->29), div. (0->0), fcn. (46->4), ass. (0->13)
	t17 = sin(qJ(3));
	t19 = cos(qJ(3));
	t29 = (r_i_i_C(1) * t19 - r_i_i_C(2) * t17) * qJD(3);
	t18 = sin(qJ(1));
	t28 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t27 = qJD(1) * t20;
	t26 = qJD(3) * t18;
	t25 = qJD(3) * t20;
	t24 = -pkin(1) - pkin(7) - r_i_i_C(3);
	t22 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19 + qJ(2);
	t21 = qJD(2) + t29;
	t1 = [t21 * t20 + (-t18 * t22 + t20 * t24) * qJD(1), t27, (-t17 * t27 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 + t19 * t27) * r_i_i_C(1), 0, 0, 0; t21 * t18 + (t18 * t24 + t20 * t22) * qJD(1), t28, (-t17 * t28 + t19 * t25) * r_i_i_C(2) + (t17 * t25 + t19 * t28) * r_i_i_C(1), 0, 0, 0; 0, 0, -t29, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (55->18), mult. (102->25), div. (0->0), fcn. (67->6), ass. (0->16)
	t25 = qJ(3) + pkin(10);
	t22 = sin(t25);
	t23 = cos(t25);
	t35 = sin(qJ(3)) * pkin(3) + r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t39 = qJD(3) * t35;
	t34 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22;
	t38 = t34 * qJD(3);
	t28 = sin(qJ(1));
	t37 = qJD(1) * t28;
	t36 = -pkin(1) - r_i_i_C(3) - qJ(4) - pkin(7);
	t33 = qJ(2) + t35;
	t32 = qJD(1) * t34;
	t31 = qJD(2) + t38;
	t30 = cos(qJ(1));
	t24 = qJD(1) * t30;
	t1 = [-t28 * qJD(4) + t31 * t30 + (-t33 * t28 + t36 * t30) * qJD(1), t24, -t28 * t39 + t30 * t32, -t37, 0, 0; t30 * qJD(4) + t31 * t28 + (t36 * t28 + t33 * t30) * qJD(1), t37, t28 * t32 + t30 * t39, t24, 0, 0; 0, 0, -t38, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:39
	% EndTime: 2019-10-10 00:56:39
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (136->32), mult. (144->44), div. (0->0), fcn. (94->8), ass. (0->28)
	t51 = qJ(3) + pkin(10);
	t41 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t51);
	t70 = t41 * qJD(3);
	t47 = qJ(5) + t51;
	t44 = cos(t47);
	t69 = r_i_i_C(1) * t44;
	t43 = sin(t47);
	t68 = r_i_i_C(2) * t43;
	t50 = qJD(3) + qJD(5);
	t67 = t43 * t50;
	t66 = t44 * t50;
	t53 = sin(qJ(1));
	t65 = qJD(1) * t53;
	t55 = cos(qJ(1));
	t48 = qJD(1) * t55;
	t64 = -pkin(1) - r_i_i_C(3) - pkin(8) - qJ(4) - pkin(7);
	t63 = r_i_i_C(1) * t67;
	t61 = qJD(1) * t69;
	t62 = t53 * t61 + (t66 * r_i_i_C(2) + t63) * t55;
	t60 = -r_i_i_C(1) * t66 + r_i_i_C(2) * t67;
	t40 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t51);
	t59 = -r_i_i_C(1) * t43 - r_i_i_C(2) * t44;
	t58 = qJD(1) * (t41 - t68);
	t57 = qJ(2) + t40 - t59;
	t56 = qJD(2) + t70 + (-t68 + t69) * t50;
	t39 = t55 * t61;
	t34 = t40 * qJD(3);
	t1 = [-t53 * qJD(4) + t56 * t55 + (-t57 * t53 + t64 * t55) * qJD(1), t48, t39 + t55 * t58 + (t59 * t50 - t34) * t53, -t65, -t53 * t63 + t39 + (-t43 * t48 - t53 * t66) * r_i_i_C(2), 0; t55 * qJD(4) + t56 * t53 + (t64 * t53 + t57 * t55) * qJD(1), t65, t55 * t34 + t53 * t58 + t62, t48, -t65 * t68 + t62, 0; 0, 0, t60 - t70, 0, t60, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:56:40
	% EndTime: 2019-10-10 00:56:40
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (399->66), mult. (424->99), div. (0->0), fcn. (318->10), ass. (0->56)
	t307 = pkin(9) + r_i_i_C(3);
	t268 = cos(qJ(6));
	t304 = r_i_i_C(1) * t268;
	t314 = pkin(5) + t304;
	t264 = qJ(3) + pkin(10);
	t254 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t264);
	t260 = qJ(5) + t264;
	t257 = cos(t260);
	t294 = t307 * t257;
	t256 = sin(t260);
	t306 = pkin(5) * t256;
	t312 = t294 - qJ(2) - t254 - t306;
	t270 = cos(qJ(1));
	t280 = qJD(1) * t256 + qJD(6);
	t311 = t280 * t270;
	t310 = t307 * t256;
	t255 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t264);
	t309 = t255 * qJD(3);
	t267 = sin(qJ(1));
	t263 = qJD(3) + qJD(5);
	t299 = t263 * t270;
	t308 = -t257 * t299 + t280 * t267;
	t305 = pkin(5) * t257;
	t265 = sin(qJ(6));
	t303 = r_i_i_C(2) * t265;
	t302 = -pkin(1) - pkin(8) - qJ(4) - pkin(7);
	t301 = t256 * t265;
	t300 = t263 * t268;
	t298 = qJD(1) * t267;
	t261 = qJD(1) * t270;
	t297 = qJD(6) * t265;
	t296 = qJD(6) * t268;
	t295 = t257 * t303;
	t293 = t263 * t301;
	t291 = t257 * t263 * t267;
	t287 = t257 * t261;
	t286 = t257 * t297;
	t285 = t257 * t296;
	t284 = r_i_i_C(2) * t293;
	t283 = qJD(1) * t257 * t304;
	t282 = r_i_i_C(2) * t285;
	t281 = -qJD(6) * t256 - qJD(1);
	t278 = t281 * t270;
	t277 = qJD(1) * (t255 - t295);
	t276 = pkin(5) * t287 + t267 * t284 + t270 * t283 + t307 * (t256 * t261 + t291);
	t275 = t256 * t300 + t286;
	t274 = t267 * t283 + t298 * t305 + (r_i_i_C(1) * t286 + t282) * t270 + (t307 * t298 + t314 * t299) * t256;
	t273 = qJD(2) + t309 + (t305 + t310) * t263;
	t272 = (-t314 * t257 + t295 - t310) * t263 + (r_i_i_C(1) * t297 + r_i_i_C(2) * t296) * t256;
	t271 = -t275 * r_i_i_C(1) - t263 * t306 - t282;
	t240 = t254 * qJD(3);
	t233 = t268 * t311 + (t257 * t300 + t281 * t265) * t267;
	t232 = t281 * t268 * t267 + (-t291 - t311) * t265;
	t231 = t265 * t278 - t308 * t268;
	t230 = t308 * t265 + t268 * t278;
	t1 = [t231 * r_i_i_C(1) + t230 * r_i_i_C(2) - t267 * qJD(4) + t273 * t270 + (t312 * t267 + t302 * t270) * qJD(1), t261, t270 * t277 + (-t240 + t271) * t267 + t276, -t298, t271 * t267 - t287 * t303 + t276, t232 * r_i_i_C(1) - t233 * r_i_i_C(2); t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t270 * qJD(4) + t273 * t267 + (t302 * t267 - t312 * t270) * qJD(1), t298, t267 * t277 + (t240 + (-r_i_i_C(2) * t301 - t294) * t263) * t270 + t274, t261, -t270 * t284 + (-t298 * t303 - t307 * t299) * t257 + t274, -t230 * r_i_i_C(1) + t231 * r_i_i_C(2); 0, 0, t272 - t309, 0, t272, t275 * r_i_i_C(2) + (-t285 + t293) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end