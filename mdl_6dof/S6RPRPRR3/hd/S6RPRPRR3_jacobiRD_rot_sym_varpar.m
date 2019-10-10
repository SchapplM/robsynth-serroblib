% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(10);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t37, 0, 0, 0, 0, 0; -t36, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:43
	% EndTime: 2019-10-10 00:49:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(10);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0, 0; -t35, 0, -t36, 0, 0, 0; 0, 0, -t44, 0, 0, 0; t36, 0, t35, 0, 0, 0; t34, 0, t37, 0, 0, 0; 0, 0, -t43, 0, 0, 0; -qJD(1) * t38, 0, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (45->18), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->18)
	t186 = sin(pkin(11));
	t189 = cos(qJ(3));
	t199 = t186 * t189;
	t187 = cos(pkin(11));
	t198 = t187 * t189;
	t188 = sin(qJ(3));
	t197 = qJD(1) * t188;
	t196 = qJD(1) * t189;
	t195 = qJD(3) * t188;
	t194 = qJD(3) * t189;
	t185 = qJ(1) + pkin(10);
	t183 = sin(t185);
	t193 = t183 * t195;
	t184 = cos(t185);
	t192 = t184 * t195;
	t191 = t183 * t194 + t184 * t197;
	t190 = t183 * t197 - t184 * t194;
	t1 = [t187 * t193 + (-t183 * t186 - t184 * t198) * qJD(1), 0, t190 * t187, 0, 0, 0; -t187 * t192 + (-t183 * t198 + t184 * t186) * qJD(1), 0, -t191 * t187, 0, 0, 0; 0, 0, -t187 * t195, 0, 0, 0; -t186 * t193 + (-t183 * t187 + t184 * t199) * qJD(1), 0, -t190 * t186, 0, 0, 0; t186 * t192 + (t183 * t199 + t184 * t187) * qJD(1), 0, t191 * t186, 0, 0, 0; 0, 0, t186 * t195, 0, 0, 0; -t191, 0, -t183 * t196 - t192, 0, 0, 0; -t190, 0, t184 * t196 - t193, 0, 0, 0; 0, 0, t194, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (168->27), mult. (173->43), div. (0->0), fcn. (173->6), ass. (0->33)
	t253 = pkin(11) + qJ(5);
	t251 = cos(t253);
	t254 = qJ(1) + pkin(10);
	t252 = cos(t254);
	t276 = t251 * t252;
	t255 = sin(qJ(3));
	t275 = qJD(1) * t255;
	t256 = cos(qJ(3));
	t274 = qJD(1) * t256;
	t273 = qJD(3) * t255;
	t272 = qJD(3) * t256;
	t271 = qJD(5) * t255;
	t270 = qJD(5) * t256;
	t269 = t252 * t275;
	t268 = t251 * t273;
	t250 = sin(t254);
	t267 = t250 * t273;
	t266 = t252 * t273;
	t249 = sin(t253);
	t265 = t249 * t271;
	t264 = t251 * t271;
	t263 = -qJD(1) + t270;
	t262 = -qJD(5) + t274;
	t261 = t263 * t249;
	t260 = t250 * t272 + t269;
	t259 = t250 * t275 - t252 * t272;
	t258 = -t251 * t272 + t265;
	t257 = t262 * t250 + t266;
	t248 = -t262 * t276 + (t261 + t268) * t250;
	t247 = t263 * t251 * t250 + (t262 * t252 - t267) * t249;
	t246 = t257 * t251 + t252 * t261;
	t245 = t257 * t249 - t263 * t276;
	t1 = [t248, 0, t259 * t251 + t252 * t265, 0, t245, 0; -t246, 0, t258 * t250 - t251 * t269, 0, -t247, 0; 0, 0, -t249 * t270 - t268, 0, -t249 * t272 - t264, 0; t247, 0, -t259 * t249 + t252 * t264, 0, t246, 0; t245, 0, t260 * t249 + t250 * t264, 0, t248, 0; 0, 0, t249 * t273 - t251 * t270, 0, t258, 0; -t260, 0, -t250 * t274 - t266, 0, 0, 0; -t259, 0, t252 * t274 - t267, 0, 0, 0; 0, 0, t272, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (348->31), mult. (233->48), div. (0->0), fcn. (233->6), ass. (0->34)
	t302 = pkin(11) + qJ(5) + qJ(6);
	t298 = sin(t302);
	t304 = qJ(1) + pkin(10);
	t301 = cos(t304);
	t324 = t298 * t301;
	t299 = cos(t302);
	t300 = sin(t304);
	t323 = t299 * t300;
	t303 = qJD(5) + qJD(6);
	t305 = sin(qJ(3));
	t322 = t303 * t305;
	t306 = cos(qJ(3));
	t321 = t303 * t306;
	t320 = qJD(1) * t305;
	t319 = qJD(1) * t306;
	t318 = qJD(3) * t305;
	t317 = qJD(3) * t306;
	t316 = t298 * t322;
	t315 = t299 * t322;
	t314 = t300 * t318;
	t313 = t301 * t318;
	t312 = -qJD(1) + t321;
	t311 = -t303 + t319;
	t310 = t298 * t312;
	t309 = t300 * t317 + t301 * t320;
	t308 = t300 * t320 - t301 * t317;
	t307 = t311 * t300 + t313;
	t294 = -t299 * t317 + t316;
	t293 = -t298 * t317 - t315;
	t292 = t300 * t310 + (-t311 * t301 + t314) * t299;
	t291 = -t298 * t314 - t303 * t324 - qJD(1) * t323 + (qJD(1) * t324 + t303 * t323) * t306;
	t290 = t307 * t299 + t301 * t310;
	t289 = -t312 * t301 * t299 + t307 * t298;
	t1 = [t292, 0, t308 * t299 + t301 * t316, 0, t289, t289; -t290, 0, -t309 * t299 + t300 * t316, 0, -t291, -t291; 0, 0, -t298 * t321 - t299 * t318, 0, t293, t293; t291, 0, -t308 * t298 + t301 * t315, 0, t290, t290; t289, 0, t309 * t298 + t300 * t315, 0, t292, t292; 0, 0, t298 * t318 - t299 * t321, 0, t294, t294; -t309, 0, -t300 * t319 - t313, 0, 0, 0; -t308, 0, t301 * t319 - t314, 0, 0, 0; 0, 0, t317, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end