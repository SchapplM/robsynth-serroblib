% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
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
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(7) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(10) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(10)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t28, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (124->30), mult. (120->40), div. (0->0), fcn. (79->7), ass. (0->28)
	t44 = qJD(3) + qJD(4);
	t43 = pkin(10) + qJ(3);
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
	t50 = -r_i_i_C(1) * t38 - pkin(3) * t40 - cos(pkin(10)) * pkin(2) - pkin(1) + t62;
	t48 = t45 * t52 + t58 * t63 + (t54 - t55) * t46;
	t32 = t45 * t54;
	t1 = [t46 * qJD(2) + t65 * t45 + (-t60 * t45 + t50 * t46) * qJD(1), t57, (t39 * t58 - t46 * t56) * pkin(3) + t48, t48, 0, 0; t45 * qJD(2) - t65 * t46 + (t50 * t45 + t60 * t46) * qJD(1), t58, t32 + (-pkin(3) * t56 - t55) * t45 + (-t51 - t64) * t57, -t46 * t52 + t32 + (-t37 * t57 - t45 * t59) * r_i_i_C(1), 0, 0; 0, 0, -t65, -t49, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:47
	% EndTime: 2019-10-10 01:28:47
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (315->47), mult. (291->61), div. (0->0), fcn. (214->9), ass. (0->39)
	t228 = pkin(10) + qJ(3);
	t226 = qJ(4) + t228;
	t223 = cos(t226);
	t261 = r_i_i_C(3) + qJ(5);
	t272 = t261 * t223;
	t222 = sin(t226);
	t221 = t222 * qJD(5);
	t229 = qJD(3) + qJD(4);
	t230 = sin(pkin(11));
	t231 = cos(pkin(11));
	t262 = r_i_i_C(2) * t230;
	t270 = r_i_i_C(1) * t231 + pkin(4);
	t239 = t270 - t262;
	t224 = sin(t228);
	t260 = pkin(3) * qJD(3);
	t252 = t224 * t260;
	t271 = (-t239 * t222 + t272) * t229 + (t230 * r_i_i_C(1) + t231 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(2)) * qJD(1) + t221 - t252;
	t253 = t229 * t262;
	t269 = (qJD(5) + t253) * t223;
	t264 = pkin(3) * t224;
	t232 = sin(qJ(1));
	t258 = t223 * t232;
	t233 = cos(qJ(1));
	t257 = t229 * t233;
	t256 = qJD(1) * t232;
	t255 = qJD(1) * t233;
	t250 = t222 * t256;
	t249 = t222 * t255;
	t247 = t261 * t222;
	t245 = t261 * t232;
	t243 = t269 * t233 + t270 * t250;
	t242 = t270 * t229;
	t241 = t270 * t233;
	t240 = t269 * t232 + t249 * t262 + t255 * t272;
	t236 = t221 + t229 * t272 + (-t242 + t253) * t222;
	t225 = cos(t228);
	t235 = -t225 * t260 + (-t223 * t270 - t247) * t229;
	t234 = qJD(2) + (-t239 * t223 - pkin(3) * t225 - cos(pkin(10)) * pkin(2) - pkin(1) - t247) * qJD(1);
	t1 = [-t271 * t232 + t234 * t233, t255, (-t222 * t262 + t264 - t272) * t256 + t235 * t233 + t243, (-t256 * t262 - t261 * t257) * t222 + (-qJD(1) * t245 - t229 * t241) * t223 + t243, t223 * t257 - t250, 0; t234 * t232 + t271 * t233, t256, (-t222 * t270 - t264) * t255 + t235 * t232 + t240, -t242 * t258 + (-qJD(1) * t241 - t229 * t245) * t222 + t240, t229 * t258 + t249, 0; 0, 0, t236 - t252, t236, t229 * t222, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:47
	% EndTime: 2019-10-10 01:28:48
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (532->72), mult. (443->102), div. (0->0), fcn. (344->11), ass. (0->57)
	t271 = pkin(11) + qJ(6);
	t265 = sin(t271);
	t267 = cos(t271);
	t272 = pkin(10) + qJ(3);
	t269 = qJ(4) + t272;
	t262 = sin(t269);
	t305 = qJD(6) * t262;
	t263 = cos(t269);
	t273 = qJD(3) + qJD(4);
	t312 = t263 * t273;
	t327 = t265 * t312 + t267 * t305;
	t264 = cos(pkin(11)) * pkin(5) + pkin(4);
	t326 = r_i_i_C(1) * t267 + t264;
	t261 = t262 * qJD(5);
	t266 = sin(t272);
	t314 = pkin(3) * qJD(3);
	t303 = t266 * t314;
	t275 = -pkin(9) - qJ(5);
	t315 = r_i_i_C(3) - t275;
	t325 = (-t262 * t264 + t315 * t263) * t273 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(2)) * qJD(1) + t261 - t303;
	t297 = t265 * t305;
	t324 = r_i_i_C(1) * t297 + t327 * r_i_i_C(2) + qJD(5) * t263;
	t277 = cos(qJ(1));
	t290 = qJD(6) * t263 - qJD(1);
	t323 = t277 * t290;
	t289 = qJD(1) * t263 - qJD(6);
	t276 = sin(qJ(1));
	t310 = t273 * t276;
	t302 = t262 * t310;
	t321 = t289 * t277 - t302;
	t319 = pkin(3) * t266;
	t317 = r_i_i_C(2) * t265;
	t316 = r_i_i_C(3) * t263;
	t311 = t263 * t275;
	t309 = t273 * t277;
	t308 = qJD(1) * t276;
	t307 = qJD(1) * t277;
	t304 = t262 * t317;
	t301 = t262 * t309;
	t299 = t262 * t308;
	t298 = t262 * t307;
	t288 = t326 * t273;
	t287 = t290 * t276;
	t285 = t275 * t302 + t324 * t276 + t298 * t317 + t307 * t316;
	t284 = -t262 * t326 - t311;
	t283 = t275 * t301 + t324 * t277 + t326 * t299 + t308 * t311;
	t282 = (-r_i_i_C(3) * t262 - t263 * t326) * t273;
	t281 = t289 * t276 + t301;
	t268 = cos(t272);
	t280 = qJD(2) + (-t315 * t262 - t263 * t264 - pkin(3) * t268 - cos(pkin(10)) * pkin(2) - pkin(1)) * qJD(1);
	t279 = -t268 * t314 + t282;
	t278 = t273 * t304 + r_i_i_C(3) * t312 + t261 - t262 * t288 + (-t273 * t275 + (-r_i_i_C(1) * t265 - r_i_i_C(2) * t267) * qJD(6)) * t263;
	t242 = t265 * t287 - t321 * t267;
	t241 = t321 * t265 + t267 * t287;
	t240 = t265 * t323 + t281 * t267;
	t239 = t281 * t265 - t267 * t323;
	t1 = [t242 * r_i_i_C(1) + t241 * r_i_i_C(2) - t325 * t276 + t280 * t277, t307, (-t304 - t316 + t319) * t308 + t279 * t277 + t283, (-r_i_i_C(3) * t309 - t308 * t317) * t262 + (-r_i_i_C(3) * t308 - t277 * t288) * t263 + t283, t263 * t309 - t299, t239 * r_i_i_C(1) + t240 * r_i_i_C(2); -t240 * r_i_i_C(1) + t239 * r_i_i_C(2) + t280 * t276 + t325 * t277, t308, t279 * t276 + (t284 - t319) * t307 + t285, t276 * t282 + t284 * t307 + t285, t263 * t310 + t298, -t241 * r_i_i_C(1) + t242 * r_i_i_C(2); 0, 0, t278 - t303, t278, t273 * t262, (-t267 * t312 + t297) * r_i_i_C(2) - t327 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end