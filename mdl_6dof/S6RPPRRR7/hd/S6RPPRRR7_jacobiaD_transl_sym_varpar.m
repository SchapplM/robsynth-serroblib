% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
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
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(10) + qJ(4);
	t21 = sin(t24);
	t22 = cos(t24);
	t36 = (r_i_i_C(1) * t22 - r_i_i_C(2) * t21) * qJD(4);
	t27 = sin(qJ(1));
	t35 = qJD(1) * t27;
	t28 = cos(qJ(1));
	t23 = qJD(1) * t28;
	t34 = qJD(4) * t27;
	t33 = qJD(4) * t28;
	t32 = -pkin(1) - r_i_i_C(3) - pkin(7) - qJ(3);
	t30 = sin(pkin(10)) * pkin(3) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t30 * t27 + t32 * t28) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0, 0; t28 * qJD(3) + t29 * t27 + (t32 * t27 + t30 * t28) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0, 0; 0, 0, 0, -t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:46
	% EndTime: 2019-10-10 00:11:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (129->32), mult. (132->44), div. (0->0), fcn. (87->7), ass. (0->28)
	t48 = pkin(10) + qJ(4);
	t44 = cos(t48);
	t67 = pkin(4) * t44;
	t45 = qJ(5) + t48;
	t42 = cos(t45);
	t66 = r_i_i_C(1) * t42;
	t41 = sin(t45);
	t65 = r_i_i_C(2) * t41;
	t49 = qJD(4) + qJD(5);
	t64 = t41 * t49;
	t63 = t42 * t49;
	t50 = sin(qJ(1));
	t62 = qJD(1) * t50;
	t51 = cos(qJ(1));
	t46 = qJD(1) * t51;
	t43 = sin(t48);
	t61 = qJD(4) * t43;
	t60 = -pkin(1) - r_i_i_C(3) - pkin(8) - pkin(7) - qJ(3);
	t59 = r_i_i_C(1) * t64;
	t58 = qJD(4) * t67;
	t57 = qJD(1) * t66;
	t56 = -r_i_i_C(1) * t63 + r_i_i_C(2) * t64;
	t55 = -r_i_i_C(1) * t41 - r_i_i_C(2) * t42;
	t54 = qJ(2) + pkin(4) * t43 + sin(pkin(10)) * pkin(3) - t55;
	t53 = t50 * t57 - t62 * t65 + (t63 * r_i_i_C(2) + t59) * t51;
	t52 = t58 + qJD(2) + (-t65 + t66) * t49;
	t39 = t51 * t57;
	t1 = [-t50 * qJD(3) + t52 * t51 + (-t54 * t50 + t60 * t51) * qJD(1), t46, -t62, t39 + (-t65 + t67) * t46 + (-pkin(4) * t61 + t55 * t49) * t50, -t50 * t59 + t39 + (-t41 * t46 - t50 * t63) * r_i_i_C(2), 0; t51 * qJD(3) + t52 * t50 + (t60 * t50 + t54 * t51) * qJD(1), t62, t46, (t44 * t62 + t51 * t61) * pkin(4) + t53, t53, 0; 0, 0, 0, t56 - t58, t56, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:11:48
	% EndTime: 2019-10-10 00:11:48
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (392->67), mult. (412->99), div. (0->0), fcn. (311->9), ass. (0->57)
	t306 = pkin(9) + r_i_i_C(3);
	t265 = cos(qJ(6));
	t303 = r_i_i_C(1) * t265;
	t312 = pkin(5) + t303;
	t261 = pkin(10) + qJ(4);
	t256 = sin(t261);
	t258 = qJ(5) + t261;
	t255 = cos(t258);
	t290 = t306 * t255;
	t254 = sin(t258);
	t305 = pkin(5) * t254;
	t310 = t290 - qJ(2) - pkin(4) * t256 - sin(pkin(10)) * pkin(3) - t305;
	t266 = cos(qJ(1));
	t276 = qJD(1) * t254 + qJD(6);
	t309 = t276 * t266;
	t308 = t306 * t254;
	t264 = sin(qJ(1));
	t262 = qJD(4) + qJD(5);
	t297 = t262 * t266;
	t307 = -t255 * t297 + t276 * t264;
	t304 = pkin(5) * t255;
	t263 = sin(qJ(6));
	t302 = r_i_i_C(2) * t263;
	t301 = -pkin(1) - pkin(8) - pkin(7) - qJ(3);
	t300 = pkin(4) * qJD(4);
	t299 = t254 * t263;
	t298 = t262 * t265;
	t296 = qJD(1) * t264;
	t259 = qJD(1) * t266;
	t295 = qJD(6) * t263;
	t294 = qJD(6) * t265;
	t293 = t255 * t302;
	t292 = t256 * t300;
	t257 = cos(t261);
	t291 = t257 * t300;
	t289 = t262 * t299;
	t287 = t255 * t262 * t264;
	t283 = t255 * t259;
	t282 = t255 * t295;
	t281 = t255 * t294;
	t280 = r_i_i_C(2) * t289;
	t279 = qJD(1) * t255 * t303;
	t278 = r_i_i_C(2) * t281;
	t277 = -qJD(6) * t254 - qJD(1);
	t274 = t277 * t266;
	t273 = pkin(5) * t283 + t264 * t280 + t266 * t279 + t306 * (t254 * t259 + t287);
	t272 = qJD(1) * (pkin(4) * t257 - t293);
	t271 = t254 * t298 + t282;
	t270 = t264 * t279 + t296 * t304 + (r_i_i_C(1) * t282 + t278) * t266 + (t306 * t296 + t312 * t297) * t254;
	t269 = t291 + qJD(2) + (t304 + t308) * t262;
	t268 = (-t312 * t255 + t293 - t308) * t262 + (r_i_i_C(1) * t295 + r_i_i_C(2) * t294) * t254;
	t267 = -t271 * r_i_i_C(1) - t262 * t305 - t278;
	t234 = t265 * t309 + (t255 * t298 + t277 * t263) * t264;
	t233 = t277 * t265 * t264 + (-t287 - t309) * t263;
	t232 = t263 * t274 - t307 * t265;
	t231 = t307 * t263 + t265 * t274;
	t1 = [t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t264 * qJD(3) + t269 * t266 + (t310 * t264 + t301 * t266) * qJD(1), t259, -t296, t266 * t272 + (t267 - t292) * t264 + t273, t267 * t264 - t283 * t302 + t273, t233 * r_i_i_C(1) - t234 * r_i_i_C(2); t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t266 * qJD(3) + t269 * t264 + (t301 * t264 - t310 * t266) * qJD(1), t296, t259, t264 * t272 + (t292 + (-r_i_i_C(2) * t299 - t290) * t262) * t266 + t270, -t266 * t280 + (-t296 * t302 - t306 * t297) * t255 + t270, -t231 * r_i_i_C(1) + t232 * r_i_i_C(2); 0, 0, 0, t268 - t291, t268, t271 * r_i_i_C(2) + (-t281 + t289) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end