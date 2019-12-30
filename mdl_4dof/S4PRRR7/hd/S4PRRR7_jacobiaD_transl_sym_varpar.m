% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PRRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(4));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(8));
	t50 = sin(pkin(8));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(4)) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:35:00
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t179 = sin(pkin(8));
	t181 = cos(pkin(8));
	t184 = sin(qJ(2));
	t182 = cos(pkin(4));
	t186 = cos(qJ(2));
	t194 = t182 * t186;
	t200 = -t179 * t184 + t181 * t194;
	t199 = pkin(6) + r_i_i_C(3);
	t180 = sin(pkin(4));
	t183 = sin(qJ(3));
	t197 = t180 * t183;
	t185 = cos(qJ(3));
	t196 = t180 * t185;
	t195 = t182 * t184;
	t192 = t183 * r_i_i_C(1) + t185 * r_i_i_C(2);
	t191 = t185 * r_i_i_C(1) - t183 * r_i_i_C(2) + pkin(2);
	t175 = t179 * t186 + t181 * t195;
	t190 = t179 * t194 + t181 * t184;
	t189 = t179 * t195 - t181 * t186;
	t188 = qJD(3) * t192;
	t187 = qJD(2) * t191;
	t172 = t190 * qJD(2);
	t170 = t200 * qJD(2);
	t1 = [0, -t199 * t172 + t189 * t187 + t190 * t188, t192 * t172 + ((-t179 * t197 + t185 * t189) * r_i_i_C(1) + (-t179 * t196 - t183 * t189) * r_i_i_C(2)) * qJD(3), 0; 0, t199 * t170 - t175 * t187 - t200 * t188, -t192 * t170 + ((-t175 * t185 + t181 * t197) * r_i_i_C(1) + (t175 * t183 + t181 * t196) * r_i_i_C(2)) * qJD(3), 0; 0, (-t186 * t188 + (-t191 * t184 + t199 * t186) * qJD(2)) * t180, -t192 * t186 * t180 * qJD(2) + ((-t182 * t183 - t184 * t196) * r_i_i_C(1) + (-t182 * t185 + t184 * t197) * r_i_i_C(2)) * qJD(3), 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:35:01
	% EndTime: 2019-12-29 12:35:02
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
	t290 = sin(qJ(3));
	t293 = cos(qJ(3));
	t289 = sin(qJ(4));
	t292 = cos(qJ(4));
	t304 = r_i_i_C(1) * t292 - r_i_i_C(2) * t289;
	t302 = pkin(3) + t304;
	t323 = pkin(7) + r_i_i_C(3);
	t325 = (t302 * t290 - t323 * t293) * qJD(3);
	t303 = t289 * r_i_i_C(1) + t292 * r_i_i_C(2);
	t320 = cos(pkin(4));
	t287 = sin(pkin(4));
	t319 = t287 * t290;
	t318 = t287 * t293;
	t294 = cos(qJ(2));
	t317 = t287 * t294;
	t291 = sin(qJ(2));
	t316 = qJD(2) * t291;
	t315 = qJD(2) * t294;
	t314 = qJD(4) * t289;
	t313 = qJD(4) * t292;
	t312 = qJD(4) * t293;
	t311 = t287 * t316;
	t310 = t287 * t315;
	t309 = t291 * t320;
	t308 = t294 * t320;
	t286 = sin(pkin(8));
	t306 = t286 * t309;
	t288 = cos(pkin(8));
	t305 = t288 * t308;
	t278 = t286 * t294 + t288 * t309;
	t301 = -t278 * t290 - t288 * t318;
	t300 = -t278 * t293 + t288 * t319;
	t280 = t288 * t294 - t306;
	t299 = -t280 * t290 + t286 * t318;
	t270 = t280 * t293 + t286 * t319;
	t298 = qJD(4) * t303;
	t279 = t286 * t308 + t288 * t291;
	t297 = -t291 * t319 + t320 * t293;
	t282 = t320 * t290 + t291 * t318;
	t296 = -t323 * t290 - t302 * t293 - pkin(2);
	t295 = t303 * t312 + t325;
	t277 = t286 * t291 - t305;
	t276 = -qJD(2) * t306 + t288 * t315;
	t275 = t279 * qJD(2);
	t274 = t278 * qJD(2);
	t273 = -qJD(2) * t305 + t286 * t316;
	t272 = t297 * qJD(3) + t293 * t310;
	t266 = t299 * qJD(3) - t275 * t293;
	t264 = t301 * qJD(3) - t273 * t293;
	t1 = [0, (-t275 * t289 + t280 * t313) * r_i_i_C(1) + (-t275 * t292 - t280 * t314) * r_i_i_C(2) - t275 * pkin(6) + t296 * t276 + t295 * t279, t323 * t266 - t299 * t298 + t302 * (-t270 * qJD(3) + t275 * t290), (-t266 * t289 + t276 * t292) * r_i_i_C(1) + (-t266 * t292 - t276 * t289) * r_i_i_C(2) + ((-t270 * t292 - t279 * t289) * r_i_i_C(1) + (t270 * t289 - t279 * t292) * r_i_i_C(2)) * qJD(4); 0, (-t273 * t289 + t278 * t313) * r_i_i_C(1) + (-t273 * t292 - t278 * t314) * r_i_i_C(2) - t273 * pkin(6) + t296 * t274 + t295 * t277, t323 * t264 - t301 * t298 + t302 * (t300 * qJD(3) + t273 * t290), (-t264 * t289 + t274 * t292) * r_i_i_C(1) + (-t264 * t292 - t274 * t289) * r_i_i_C(2) + ((-t277 * t289 + t292 * t300) * r_i_i_C(1) + (-t277 * t292 - t289 * t300) * r_i_i_C(2)) * qJD(4); 0, ((t296 * qJD(2) + t304 * qJD(4)) * t291 + (qJD(2) * pkin(6) - t325 + t303 * (qJD(2) - t312)) * t294) * t287, t323 * t272 - t297 * t298 + t302 * (-t282 * qJD(3) - t290 * t310), (-t272 * t289 + t292 * t311) * r_i_i_C(1) + (-t272 * t292 - t289 * t311) * r_i_i_C(2) + ((-t282 * t292 + t289 * t317) * r_i_i_C(1) + (t282 * t289 + t292 * t317) * r_i_i_C(2)) * qJD(4);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end