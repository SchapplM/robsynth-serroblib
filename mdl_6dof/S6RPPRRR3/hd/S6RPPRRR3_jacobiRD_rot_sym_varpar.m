% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPPRRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
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
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
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
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (5->2), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t16 = qJ(1) + pkin(10);
	t17 = qJD(1) * sin(t16);
	t14 = qJD(1) * cos(t16);
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; t17, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t17, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t42 = sin(qJ(4));
	t47 = qJD(1) * t42;
	t43 = cos(qJ(4));
	t46 = qJD(1) * t43;
	t45 = qJD(4) * t42;
	t44 = qJD(4) * t43;
	t41 = qJ(1) + pkin(10);
	t40 = cos(t41);
	t39 = sin(t41);
	t38 = -t39 * t45 + t40 * t46;
	t37 = t39 * t44 + t40 * t47;
	t36 = t39 * t46 + t40 * t45;
	t35 = -t39 * t47 + t40 * t44;
	t1 = [t35, 0, 0, t38, 0, 0; t37, 0, 0, t36, 0, 0; 0, 0, 0, -t44, 0, 0; -t36, 0, 0, -t37, 0, 0; t38, 0, 0, t35, 0, 0; 0, 0, 0, t45, 0, 0; -qJD(1) * t40, 0, 0, 0, 0, 0; -qJD(1) * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:50
	% EndTime: 2019-10-10 00:04:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t239 = sin(qJ(5));
	t240 = sin(qJ(4));
	t258 = qJD(1) * t240;
	t247 = qJD(5) + t258;
	t260 = t239 * t247;
	t241 = cos(qJ(5));
	t259 = t241 * t247;
	t242 = cos(qJ(4));
	t257 = qJD(1) * t242;
	t256 = qJD(4) * t240;
	t255 = qJD(4) * t242;
	t254 = qJD(5) * t240;
	t253 = qJD(5) * t242;
	t252 = t239 * t257;
	t251 = t241 * t257;
	t250 = t239 * t255;
	t249 = t241 * t255;
	t248 = -qJD(1) - t254;
	t246 = t239 * t253 + t241 * t256;
	t245 = t239 * t256 - t241 * t253;
	t244 = t248 * t241 - t250;
	t243 = t248 * t239 + t249;
	t238 = qJ(1) + pkin(10);
	t237 = cos(t238);
	t236 = sin(t238);
	t235 = t243 * t236 + t237 * t259;
	t234 = t244 * t236 - t237 * t260;
	t233 = -t236 * t259 + t243 * t237;
	t232 = t236 * t260 + t244 * t237;
	t1 = [t233, 0, 0, -t246 * t236 + t237 * t251, t234, 0; t235, 0, 0, t236 * t251 + t246 * t237, -t232, 0; 0, 0, 0, t239 * t254 - t249, t245, 0; t232, 0, 0, t245 * t236 - t237 * t252, -t235, 0; t234, 0, 0, -t236 * t252 - t245 * t237, t233, 0; 0, 0, 0, t241 * t254 + t250, t246, 0; t236 * t257 + t237 * t256, 0, 0, t236 * t255 + t237 * t258, 0, 0; t236 * t256 - t237 * t257, 0, 0, t236 * t258 - t237 * t255, 0, 0; 0, 0, 0, -t256, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:51
	% EndTime: 2019-10-10 00:04:51
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (265->32), mult. (233->48), div. (0->0), fcn. (233->6), ass. (0->34)
	t304 = qJ(1) + pkin(10);
	t300 = cos(t304);
	t303 = qJD(5) + qJD(6);
	t306 = sin(qJ(4));
	t321 = qJD(1) * t306;
	t308 = t303 + t321;
	t325 = t308 * t300;
	t299 = sin(t304);
	t305 = qJ(5) + qJ(6);
	t301 = sin(t305);
	t324 = t299 * t301;
	t323 = t303 * t306;
	t307 = cos(qJ(4));
	t322 = t303 * t307;
	t320 = qJD(1) * t307;
	t319 = qJD(4) * t306;
	t318 = qJD(4) * t307;
	t317 = t301 * t323;
	t302 = cos(t305);
	t316 = t302 * t323;
	t315 = t301 * t320;
	t314 = t302 * t320;
	t313 = t301 * t318;
	t312 = t302 * t318;
	t311 = t299 * t318;
	t310 = t300 * t318;
	t309 = -qJD(1) - t323;
	t293 = t301 * t322 + t302 * t319;
	t292 = t301 * t319 - t302 * t322;
	t291 = -t299 * t317 - qJD(1) * t324 + (t311 + t325) * t302;
	t290 = -t301 * t325 + (t309 * t302 - t313) * t299;
	t289 = -t308 * t302 * t299 + (t309 * t301 + t312) * t300;
	t288 = -t301 * t310 + t308 * t324 + (-qJD(1) * t302 - t316) * t300;
	t1 = [t289, 0, 0, -t293 * t299 + t300 * t314, t290, t290; t291, 0, 0, t293 * t300 + t299 * t314, -t288, -t288; 0, 0, 0, -t312 + t317, t292, t292; t288, 0, 0, t292 * t299 - t300 * t315, -t291, -t291; t290, 0, 0, -t292 * t300 - t299 * t315, t289, t289; 0, 0, 0, t313 + t316, t293, t293; t299 * t320 + t300 * t319, 0, 0, t300 * t321 + t311, 0, 0; t299 * t319 - t300 * t320, 0, 0, t299 * t321 - t310, 0, 0; 0, 0, 0, -t319, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end