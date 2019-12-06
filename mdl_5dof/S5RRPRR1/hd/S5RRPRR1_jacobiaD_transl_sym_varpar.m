% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiaD_transl_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (15->13), mult. (56->29), div. (0->0), fcn. (36->4), ass. (0->12)
	t16 = sin(qJ(1));
	t25 = qJD(1) * t16;
	t18 = cos(qJ(1));
	t24 = qJD(1) * t18;
	t23 = qJD(2) * t16;
	t22 = qJD(2) * t18;
	t15 = sin(qJ(2));
	t17 = cos(qJ(2));
	t21 = -r_i_i_C(1) * t17 + r_i_i_C(2) * t15;
	t20 = r_i_i_C(1) * t15 + r_i_i_C(2) * t17;
	t19 = t20 * qJD(2);
	t1 = [t20 * t23 + (-r_i_i_C(3) * t16 + t21 * t18) * qJD(1), (t15 * t22 + t17 * t25) * r_i_i_C(2) + (t15 * t25 - t17 * t22) * r_i_i_C(1), 0, 0, 0; -t18 * t19 + (r_i_i_C(3) * t18 + t21 * t16) * qJD(1), (t15 * t23 - t17 * t24) * r_i_i_C(2) + (-t15 * t24 - t17 * t23) * r_i_i_C(1), 0, 0, 0; 0, -t19, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (28->13), mult. (90->23), div. (0->0), fcn. (59->4), ass. (0->13)
	t18 = sin(qJ(2));
	t20 = cos(qJ(2));
	t30 = pkin(1) + r_i_i_C(1);
	t31 = r_i_i_C(2) * t20 + t30 * t18;
	t28 = r_i_i_C(3) + qJ(3);
	t21 = cos(qJ(1));
	t27 = qJD(1) * t21;
	t26 = r_i_i_C(2) * t18 - t30 * t20;
	t19 = sin(qJ(1));
	t24 = t31 * t19;
	t23 = qJD(2) * t26;
	t22 = t31 * qJD(2);
	t1 = [t21 * qJD(3) + qJD(2) * t24 + (-t28 * t19 + t26 * t21) * qJD(1), qJD(1) * t24 + t21 * t23, t27, 0, 0; t19 * qJD(3) - t21 * t22 + (t26 * t19 + t28 * t21) * qJD(1), t19 * t23 - t27 * t31, qJD(1) * t19, 0, 0; 0, -t22, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (90->28), mult. (118->39), div. (0->0), fcn. (77->6), ass. (0->29)
	t40 = sin(qJ(2));
	t37 = qJD(2) + qJD(4);
	t38 = qJ(2) + qJ(4);
	t36 = cos(t38);
	t61 = r_i_i_C(2) * t36;
	t35 = sin(t38);
	t63 = r_i_i_C(1) * t35;
	t49 = t61 + t63;
	t48 = t49 * t37;
	t44 = pkin(2) + pkin(1);
	t55 = qJD(2) * t44;
	t64 = t40 * t55 + t48;
	t62 = r_i_i_C(2) * t35;
	t60 = r_i_i_C(3) + pkin(3) + qJ(3);
	t59 = t36 * t37;
	t58 = t40 * t44;
	t41 = sin(qJ(1));
	t57 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t56 = qJD(1) * t43;
	t54 = r_i_i_C(1) * t59;
	t53 = t37 * t62;
	t51 = qJD(1) * t61;
	t52 = t41 * t51 + t43 * t53 + t57 * t63;
	t42 = cos(qJ(2));
	t47 = -t42 * t55 - t54;
	t46 = -r_i_i_C(1) * t36 - t42 * t44 + t62;
	t31 = t41 * t53;
	t1 = [t43 * qJD(3) + t64 * t41 + (-t41 * t60 + t43 * t46) * qJD(1), t43 * t47 + t57 * t58 + t52, t56, -t43 * t54 + t52, 0; t41 * qJD(3) - t64 * t43 + (t41 * t46 + t43 * t60) * qJD(1), t31 + t47 * t41 + (-t49 - t58) * t56, t57, -t43 * t51 + t31 + (-t35 * t56 - t41 * t59) * r_i_i_C(1), 0; 0, -t64, 0, -t48, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:05
	% EndTime: 2019-12-05 18:26:06
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (243->54), mult. (358->88), div. (0->0), fcn. (275->8), ass. (0->51)
	t255 = qJ(2) + qJ(4);
	t252 = sin(t255);
	t257 = sin(qJ(5));
	t260 = cos(qJ(5));
	t289 = qJD(5) * t260;
	t253 = cos(t255);
	t254 = qJD(2) + qJD(4);
	t297 = t253 * t254;
	t304 = t252 * t289 + t257 * t297;
	t300 = pkin(4) + r_i_i_C(3);
	t303 = t300 * t297;
	t288 = r_i_i_C(2) * t252 * t257;
	t302 = t300 * t253 + t288;
	t290 = qJD(5) * t257;
	t279 = t252 * t290;
	t301 = r_i_i_C(1) * t279 + t304 * r_i_i_C(2);
	t299 = r_i_i_C(1) * t260;
	t259 = sin(qJ(1));
	t298 = t252 * t259;
	t296 = t254 * t260;
	t258 = sin(qJ(2));
	t263 = pkin(2) + pkin(1);
	t295 = t258 * t263;
	t262 = cos(qJ(1));
	t294 = t260 * t262;
	t293 = qJD(1) * t259;
	t292 = qJD(1) * t262;
	t291 = qJD(2) * t263;
	t287 = qJD(1) * t299;
	t286 = t300 * t252;
	t285 = t300 * t254;
	t284 = t252 * t296;
	t282 = t253 * t296;
	t280 = t258 * t291;
	t274 = qJD(5) * t253 - qJD(1);
	t273 = qJD(1) * t253 - qJD(5);
	t272 = t301 * t262 + t287 * t298;
	t271 = t274 * t257;
	t270 = t301 * t259 + t302 * t292;
	t261 = cos(qJ(2));
	t269 = -t261 * t263 - t286;
	t267 = (-t253 * t299 - t286) * t254;
	t266 = t252 * t254 * t262 + t273 * t259;
	t265 = -t261 * t291 + t267;
	t264 = -t253 * r_i_i_C(2) * t289 + t254 * t288 + (-t253 * t290 - t284) * r_i_i_C(1) + t303;
	t256 = -pkin(3) - qJ(3);
	t238 = -t273 * t294 + (t271 + t284) * t259;
	t237 = t274 * t260 * t259 + (-t254 * t298 + t273 * t262) * t257;
	t236 = t266 * t260 + t262 * t271;
	t235 = t266 * t257 - t274 * t294;
	t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) + t262 * qJD(3) + (-t253 * t285 + t280) * t259 + (t256 * t259 + t269 * t262) * qJD(1), t265 * t262 + (-t302 + t295) * t293 + t272, t292, t262 * t267 - t293 * t302 + t272, t235 * r_i_i_C(1) + t236 * r_i_i_C(2); -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t259 * qJD(3) + (-t280 + t303) * t262 + (-t256 * t262 + t269 * t259) * qJD(1), (-t252 * t299 - t295) * t292 + t265 * t259 + t270, t293, -t259 * r_i_i_C(1) * t282 + (-t259 * t285 - t262 * t287) * t252 + t270, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2); 0, t264 - t280, 0, t264, (t279 - t282) * r_i_i_C(2) - t304 * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end