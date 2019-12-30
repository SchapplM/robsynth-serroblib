% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t18 * t28 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:58
	% DurationCPUTime: 0.20s
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
	t58 = r_i_i_C(3) + pkin(7) + pkin(6);
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
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0; 0, -t63, -t47, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:59
	% EndTime: 2019-12-29 20:07:59
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (179->38), mult. (217->49), div. (0->0), fcn. (148->6), ass. (0->35)
	t189 = qJ(2) + qJ(3);
	t187 = cos(t189);
	t218 = r_i_i_C(3) + qJ(4);
	t202 = t218 * t187;
	t186 = sin(t189);
	t184 = t186 * qJD(4);
	t188 = qJD(2) + qJD(3);
	t190 = sin(qJ(2));
	t217 = pkin(2) * qJD(2);
	t210 = t190 * t217;
	t221 = pkin(3) - r_i_i_C(2);
	t226 = (-t186 * t221 + t202) * t188 + (r_i_i_C(1) + pkin(7) + pkin(6)) * qJD(1) + t184 - t210;
	t191 = sin(qJ(1));
	t214 = t188 * t191;
	t209 = t187 * t214;
	t193 = cos(qJ(1));
	t212 = qJD(1) * t193;
	t224 = t186 * t212 + t209;
	t220 = pkin(2) * t190;
	t216 = t187 * t188;
	t215 = t188 * t186;
	t213 = qJD(1) * t191;
	t211 = qJD(4) * t187;
	t208 = t193 * t216;
	t205 = t186 * t213;
	t207 = pkin(3) * t205 + r_i_i_C(2) * t208 + t193 * t211;
	t203 = t218 * t186;
	t201 = t224 * r_i_i_C(2) + t191 * t211 + t212 * t202;
	t200 = -r_i_i_C(2) * t186 - t202;
	t198 = -t221 * t215 + t218 * t216 + t184;
	t197 = (-pkin(3) * t187 - t203) * t188;
	t192 = cos(qJ(2));
	t196 = qJD(1) * (-t192 * pkin(2) - t187 * t221 - pkin(1) - t203);
	t195 = -t192 * t217 + t197;
	t1 = [-t226 * t191 + t193 * t196, t195 * t193 + (t200 + t220) * t213 + t207, t193 * t197 + t200 * t213 + t207, -t205 + t208, 0; t191 * t196 + t226 * t193, (-pkin(3) * t186 - t220) * t212 + t195 * t191 + t201, -pkin(3) * t209 + (-pkin(3) * t212 - t214 * t218) * t186 + t201, t224, 0; 0, t198 - t210, t198, t215, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:08:04
	% EndTime: 2019-12-29 20:08:05
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (334->60), mult. (461->87), div. (0->0), fcn. (348->8), ass. (0->55)
	t259 = qJ(2) + qJ(3);
	t257 = cos(t259);
	t260 = sin(qJ(5));
	t263 = cos(qJ(5));
	t311 = r_i_i_C(1) * t260 + r_i_i_C(2) * t263 + qJ(4);
	t315 = t257 * t311;
	t258 = qJD(2) + qJD(3);
	t302 = t257 * t258;
	t251 = qJ(4) * t302;
	t256 = sin(t259);
	t293 = pkin(3) + pkin(8) + r_i_i_C(3);
	t282 = t293 * t258;
	t261 = sin(qJ(2));
	t303 = pkin(2) * qJD(2);
	t292 = t261 * t303;
	t314 = (qJD(4) - t282) * t256 + (pkin(4) + pkin(7) + pkin(6)) * qJD(1) + t251 - t292;
	t294 = qJD(5) * t263;
	t285 = t257 * t294;
	t313 = r_i_i_C(1) * t285 + qJD(4) * t257;
	t280 = qJD(5) * t256 + qJD(1);
	t310 = t263 * t280;
	t309 = t280 * t260;
	t307 = pkin(2) * t261;
	t301 = t258 * t260;
	t300 = t258 * t263;
	t265 = cos(qJ(1));
	t299 = t258 * t265;
	t262 = sin(qJ(1));
	t298 = qJD(1) * t262;
	t297 = qJD(1) * t265;
	t295 = qJD(5) * t260;
	t291 = t257 * t300;
	t290 = t262 * t302;
	t289 = t257 * t299;
	t287 = t256 * t298;
	t286 = t257 * t295;
	t284 = t293 * t256;
	t283 = t293 * t257;
	t279 = -qJD(1) * t256 - qJD(5);
	t278 = t313 * t262 + t297 * t315;
	t277 = t313 * t265 + t293 * t287;
	t276 = t279 * t265;
	t273 = t311 * t256;
	t272 = -r_i_i_C(2) * t295 - t282;
	t264 = cos(qJ(2));
	t271 = qJD(1) * (-t264 * pkin(2) - qJ(4) * t256 - pkin(1) - t283);
	t270 = t279 * t262 + t289;
	t269 = t257 * r_i_i_C(1) * t301 + r_i_i_C(2) * t291 + t251 + (r_i_i_C(1) * t294 + qJD(4) + t272) * t256;
	t268 = -r_i_i_C(2) * t286 + (-t283 - t273) * t258;
	t267 = -t264 * t303 + t268;
	t239 = t270 * t260 + t265 * t310;
	t238 = t270 * t263 - t265 * t309;
	t237 = -t262 * t310 + (t276 - t290) * t260;
	t236 = t263 * t276 + (-t291 + t309) * t262;
	t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t314 * t262 + t265 * t271, (t307 - t315) * t298 + t267 * t265 + t277, -t273 * t299 + (t272 * t265 - t298 * t311) * t257 + t277, -t287 + t289, t238 * r_i_i_C(1) - t239 * r_i_i_C(2); t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t262 * t271 + t314 * t265, (-t284 - t307) * t297 + t267 * t262 + t278, t268 * t262 - t284 * t297 + t278, t256 * t297 + t290, -t236 * r_i_i_C(1) + t237 * r_i_i_C(2); 0, t269 - t292, t269, t258 * t256, (-t256 * t301 + t285) * r_i_i_C(2) + (t256 * t300 + t286) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end