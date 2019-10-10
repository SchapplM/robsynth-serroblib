% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:27
	% EndTime: 2019-10-10 00:51:27
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
	% StartTime: 2019-10-10 00:51:28
	% EndTime: 2019-10-10 00:51:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (27->8), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t170 = sin(qJ(3));
	t175 = qJD(1) * t170;
	t171 = cos(qJ(3));
	t174 = qJD(1) * t171;
	t173 = qJD(3) * t170;
	t172 = qJD(3) * t171;
	t169 = qJ(1) + pkin(10);
	t168 = cos(t169);
	t167 = sin(t169);
	t166 = -t167 * t173 + t168 * t174;
	t165 = t167 * t172 + t168 * t175;
	t164 = t167 * t174 + t168 * t173;
	t163 = -t167 * t175 + t168 * t172;
	t1 = [-qJD(1) * t167, 0, 0, 0, 0, 0; qJD(1) * t168, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t166, 0, t163, 0, 0, 0; t164, 0, t165, 0, 0, 0; 0, 0, t173, 0, 0, 0; -t165, 0, -t164, 0, 0, 0; t163, 0, t166, 0, 0, 0; 0, 0, t172, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:28
	% EndTime: 2019-10-10 00:51:29
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (109->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t242 = sin(qJ(5));
	t243 = sin(qJ(3));
	t257 = qJD(5) * t243;
	t251 = qJD(1) + t257;
	t244 = cos(qJ(5));
	t245 = cos(qJ(3));
	t258 = qJD(3) * t245;
	t252 = t244 * t258;
	t263 = t251 * t242 - t252;
	t253 = t242 * t258;
	t262 = t251 * t244 + t253;
	t261 = qJD(1) * t243;
	t260 = qJD(1) * t245;
	t259 = qJD(3) * t243;
	t256 = qJD(5) * t245;
	t255 = t242 * t260;
	t254 = t244 * t260;
	t250 = -qJD(5) - t261;
	t249 = t250 * t242;
	t248 = t250 * t244;
	t247 = t242 * t256 + t244 * t259;
	t246 = -t242 * t259 + t244 * t256;
	t241 = qJ(1) + pkin(10);
	t240 = cos(t241);
	t239 = sin(t241);
	t238 = t239 * t249 + t262 * t240;
	t237 = t239 * t248 - t263 * t240;
	t236 = -t262 * t239 + t240 * t249;
	t235 = t263 * t239 + t240 * t248;
	t1 = [t236, 0, -t239 * t255 + t246 * t240, 0, t237, 0; t238, 0, t246 * t239 + t240 * t255, 0, -t235, 0; 0, 0, t244 * t257 + t253, 0, t247, 0; t235, 0, -t239 * t254 - t247 * t240, 0, -t238, 0; t237, 0, -t247 * t239 + t240 * t254, 0, t236, 0; 0, 0, -t242 * t257 + t252, 0, t246, 0; t239 * t259 - t240 * t260, 0, t239 * t261 - t240 * t258, 0, 0, 0; -t239 * t260 - t240 * t259, 0, -t239 * t258 - t240 * t261, 0, 0, 0; 0, 0, -t259, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:51:29
	% EndTime: 2019-10-10 00:51:29
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (265->30), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->32)
	t307 = qJ(1) + pkin(10);
	t302 = sin(t307);
	t303 = cos(t307);
	t310 = cos(qJ(3));
	t320 = qJD(3) * t310;
	t309 = sin(qJ(3));
	t323 = qJD(1) * t309;
	t328 = -t302 * t320 - t303 * t323;
	t308 = qJ(5) + qJ(6);
	t305 = cos(t308);
	t306 = qJD(5) + qJD(6);
	t325 = t306 * t309;
	t312 = qJD(1) + t325;
	t304 = sin(t308);
	t315 = t304 * t320;
	t327 = t305 * t312 + t315;
	t326 = t302 * t304;
	t324 = t306 * t310;
	t322 = qJD(1) * t310;
	t321 = qJD(3) * t309;
	t319 = t302 * t323;
	t317 = t304 * t322;
	t316 = t305 * t322;
	t314 = t305 * t320;
	t311 = -t306 - t323;
	t296 = t304 * t324 + t305 * t321;
	t295 = -t304 * t321 + t305 * t324;
	t294 = t327 * t303 - t304 * t319 - t306 * t326;
	t293 = t311 * t305 * t302 + (-t304 * t312 + t314) * t303;
	t292 = t311 * t304 * t303 - t327 * t302;
	t291 = t312 * t326 + (-t303 * t306 + t328) * t305;
	t1 = [t292, 0, t295 * t303 - t302 * t317, 0, t293, t293; t294, 0, t295 * t302 + t303 * t317, 0, -t291, -t291; 0, 0, t305 * t325 + t315, 0, t296, t296; t291, 0, -t296 * t303 - t302 * t316, 0, -t294, -t294; t293, 0, -t296 * t302 + t303 * t316, 0, t292, t292; 0, 0, -t304 * t325 + t314, 0, t295, t295; t302 * t321 - t303 * t322, 0, -t303 * t320 + t319, 0, 0, 0; -t302 * t322 - t303 * t321, 0, t328, 0, 0, 0; 0, 0, -t321, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end