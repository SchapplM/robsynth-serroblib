% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(8);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0; -t35, 0, -t36, 0, 0; 0, 0, -t44, 0, 0; t36, 0, t35, 0, 0; t34, 0, t37, 0, 0; 0, 0, -t43, 0, 0; -qJD(1) * t38, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:18
	% EndTime: 2019-12-31 18:46:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (108->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t245 = cos(qJ(4));
	t246 = cos(qJ(3));
	t261 = qJD(1) * t246;
	t251 = -qJD(4) + t261;
	t264 = t245 * t251;
	t257 = qJD(4) * t246;
	t252 = -qJD(1) + t257;
	t243 = sin(qJ(4));
	t244 = sin(qJ(3));
	t260 = qJD(3) * t244;
	t254 = t243 * t260;
	t263 = t252 * t245 - t254;
	t262 = qJD(1) * t244;
	t259 = qJD(3) * t246;
	t258 = qJD(4) * t244;
	t256 = t243 * t262;
	t255 = t245 * t262;
	t253 = t245 * t260;
	t250 = t251 * t243;
	t249 = t243 * t258 - t245 * t259;
	t248 = t243 * t259 + t245 * t258;
	t247 = t252 * t243 + t253;
	t242 = qJ(1) + pkin(8);
	t241 = cos(t242);
	t240 = sin(t242);
	t239 = t247 * t240 - t241 * t264;
	t238 = t263 * t240 + t241 * t250;
	t237 = t240 * t264 + t247 * t241;
	t236 = t240 * t250 - t263 * t241;
	t1 = [t239, 0, t240 * t255 + t249 * t241, t236, 0; -t237, 0, t249 * t240 - t241 * t255, -t238, 0; 0, 0, -t243 * t257 - t253, -t248, 0; t238, 0, -t240 * t256 + t248 * t241, t237, 0; t236, 0, t248 * t240 + t241 * t256, t239, 0; 0, 0, -t245 * t257 + t254, t249, 0; -t240 * t259 - t241 * t262, 0, -t240 * t261 - t241 * t260, 0, 0; -t240 * t262 + t241 * t259, 0, -t240 * t260 + t241 * t261, 0, 0; 0, 0, t259, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:18
	% EndTime: 2019-12-31 18:46:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (108->23), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->30)
	t295 = sin(qJ(4));
	t298 = cos(qJ(3));
	t313 = qJD(1) * t298;
	t303 = -qJD(4) + t313;
	t316 = t295 * t303;
	t297 = cos(qJ(4));
	t315 = t297 * t303;
	t296 = sin(qJ(3));
	t314 = qJD(1) * t296;
	t312 = qJD(3) * t296;
	t311 = qJD(3) * t298;
	t310 = qJD(4) * t296;
	t309 = qJD(4) * t298;
	t308 = t295 * t314;
	t307 = t297 * t314;
	t306 = t295 * t312;
	t305 = t297 * t312;
	t304 = qJD(1) - t309;
	t302 = -t295 * t310 + t297 * t311;
	t301 = -t295 * t311 - t297 * t310;
	t300 = t304 * t297 + t306;
	t299 = t304 * t295 - t305;
	t294 = qJ(1) + pkin(8);
	t293 = cos(t294);
	t292 = sin(t294);
	t291 = t299 * t292 + t293 * t315;
	t290 = t300 * t292 - t293 * t316;
	t289 = -t292 * t315 + t299 * t293;
	t288 = t292 * t316 + t300 * t293;
	t1 = [-t291, 0, t292 * t307 - t302 * t293, t288, 0; t289, 0, -t302 * t292 - t293 * t307, t290, 0; 0, 0, -t295 * t309 - t305, t301, 0; -t292 * t311 - t293 * t314, 0, -t292 * t313 - t293 * t312, 0, 0; -t292 * t314 + t293 * t311, 0, -t292 * t312 + t293 * t313, 0, 0; 0, 0, t311, 0, 0; t290, 0, t292 * t308 + t301 * t293, t289, 0; -t288, 0, t301 * t292 - t293 * t308, t291, 0; 0, 0, t297 * t309 - t306, t302, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end