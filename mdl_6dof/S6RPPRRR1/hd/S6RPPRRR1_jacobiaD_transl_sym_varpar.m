% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
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
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
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
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(11)) + r_i_i_C(2) * sin(pkin(11)) - pkin(2);
	t15 = qJ(1) + pkin(10);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->22), mult. (74->34), div. (0->0), fcn. (48->7), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(7) + qJ(3);
	t31 = qJ(1) + pkin(10);
	t27 = sin(t31);
	t39 = qJD(1) * t27;
	t29 = cos(t31);
	t38 = qJD(1) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(11) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(11)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [t29 * qJD(3) + t35 * t37 + (-cos(qJ(1)) * pkin(1) - t40 * t27 + t34 * t29) * qJD(1), 0, t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0, 0; t27 * qJD(3) - t29 * t33 + (-sin(qJ(1)) * pkin(1) + t40 * t29 + t34 * t27) * qJD(1), 0, t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0, 0; 0, 0, 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (166->33), mult. (124->42), div. (0->0), fcn. (81->9), ass. (0->29)
	t50 = qJD(4) + qJD(5);
	t49 = pkin(11) + qJ(4);
	t47 = qJ(5) + t49;
	t42 = cos(t47);
	t66 = r_i_i_C(2) * t42;
	t41 = sin(t47);
	t68 = r_i_i_C(1) * t41;
	t56 = t66 + t68;
	t54 = t56 * t50;
	t43 = sin(t49);
	t69 = pkin(4) * t43;
	t70 = qJD(4) * t69 + t54;
	t67 = r_i_i_C(2) * t41;
	t65 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
	t64 = t42 * t50;
	t51 = qJ(1) + pkin(10);
	t44 = sin(t51);
	t63 = qJD(1) * t44;
	t46 = cos(t51);
	t62 = qJD(1) * t46;
	t45 = cos(t49);
	t61 = qJD(4) * t45;
	t60 = r_i_i_C(1) * t64;
	t59 = t50 * t67;
	t57 = qJD(1) * t66;
	t55 = -r_i_i_C(1) * t42 - pkin(4) * t45 - cos(pkin(11)) * pkin(3) - pkin(2) + t67;
	t53 = t44 * t57 + t63 * t68 + (t59 - t60) * t46;
	t36 = t44 * t59;
	t1 = [t46 * qJD(3) + t70 * t44 + (-cos(qJ(1)) * pkin(1) - t65 * t44 + t55 * t46) * qJD(1), 0, t62, (t43 * t63 - t46 * t61) * pkin(4) + t53, t53, 0; t44 * qJD(3) - t70 * t46 + (-sin(qJ(1)) * pkin(1) + t65 * t46 + t55 * t44) * qJD(1), 0, t63, t36 + (-pkin(4) * t61 - t60) * t44 + (-t56 - t69) * t62, -t46 * t57 + t36 + (-t41 * t62 - t44 * t64) * r_i_i_C(1), 0; 0, 0, 0, -t70, -t54, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:25
	% EndTime: 2019-10-10 00:01:25
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (503->63), mult. (404->95), div. (0->0), fcn. (305->11), ass. (0->54)
	t271 = pkin(11) + qJ(4);
	t269 = qJ(5) + t271;
	t263 = sin(t269);
	t275 = cos(qJ(6));
	t320 = r_i_i_C(1) * t275 + pkin(5);
	t325 = t263 * t320;
	t264 = cos(t269);
	t302 = qJD(6) * t275;
	t272 = qJD(4) + qJD(5);
	t274 = sin(qJ(6));
	t307 = t272 * t274;
	t324 = t263 * t302 + t264 * t307;
	t314 = pkin(9) + r_i_i_C(3);
	t322 = t314 * t264;
	t323 = (-pkin(5) * t263 + t322) * t272;
	t265 = sin(t271);
	t309 = pkin(4) * qJD(4);
	t301 = t265 * t309;
	t321 = -t301 + t323;
	t303 = qJD(6) * t274;
	t290 = t263 * t303;
	t317 = r_i_i_C(1) * t290 + t324 * r_i_i_C(2);
	t284 = qJD(1) * t264 - qJD(6);
	t316 = t275 * t284;
	t285 = qJD(6) * t264 - qJD(1);
	t296 = t263 * t307;
	t315 = t285 * t275 - t296;
	t313 = pkin(4) * t265;
	t310 = r_i_i_C(2) * t274;
	t306 = t272 * t275;
	t273 = qJ(1) + pkin(10);
	t266 = sin(t273);
	t305 = qJD(1) * t266;
	t268 = cos(t273);
	t304 = qJD(1) * t268;
	t300 = qJD(1) * t310;
	t299 = t314 * t263;
	t297 = t314 * t272;
	t295 = t263 * t306;
	t283 = t320 * t272;
	t282 = t317 * t268 + t305 * t325;
	t281 = t284 * t274;
	t280 = t268 * t263 * t300 + t317 * t266 + t304 * t322;
	t267 = cos(t271);
	t279 = -pkin(5) * t264 - pkin(4) * t267 - cos(pkin(11)) * pkin(3) - pkin(2) - t299;
	t278 = t285 * t274 + t295;
	t277 = -t267 * t309 + (-t264 * t320 - t299) * t272;
	t276 = (-t264 * t303 - t295) * r_i_i_C(1) + (-t264 * t302 + t296) * r_i_i_C(2) + t323;
	t270 = -pkin(8) - pkin(7) - qJ(3);
	t247 = t278 * t266 - t268 * t316;
	t246 = t315 * t266 + t268 * t281;
	t245 = t266 * t316 + t278 * t268;
	t244 = t266 * t281 - t315 * t268;
	t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t268 * qJD(3) - t321 * t266 + (-cos(qJ(1)) * pkin(1) + t266 * t270 + t279 * t268) * qJD(1), 0, t304, (-t263 * t310 + t313 - t322) * t305 + t277 * t268 + t282, (-t266 * t300 - t268 * t297) * t263 + (-t268 * t283 - t314 * t305) * t264 + t282, t244 * r_i_i_C(1) + t245 * r_i_i_C(2); -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t266 * qJD(3) + t321 * t268 + (-sin(qJ(1)) * pkin(1) - t268 * t270 + t279 * t266) * qJD(1), 0, t305, (-t313 - t325) * t304 + t277 * t266 + t280, -t266 * t264 * t283 + (-t266 * t297 - t304 * t320) * t263 + t280, -t246 * r_i_i_C(1) + t247 * r_i_i_C(2); 0, 0, 0, t276 - t301, t276, (-t264 * t306 + t290) * r_i_i_C(2) - t324 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end