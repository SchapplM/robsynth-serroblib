% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRPRP10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
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
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:31
	% EndTime: 2019-10-10 00:44:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:32
	% EndTime: 2019-10-10 00:44:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t161 = sin(qJ(1));
	t168 = qJD(1) * t161;
	t163 = cos(qJ(1));
	t167 = qJD(1) * t163;
	t160 = sin(qJ(3));
	t166 = qJD(3) * t160;
	t162 = cos(qJ(3));
	t165 = qJD(3) * t162;
	t164 = qJD(3) * t163;
	t159 = t161 * t166 - t162 * t167;
	t158 = t160 * t167 + t161 * t165;
	t157 = t160 * t164 + t162 * t168;
	t156 = t160 * t168 - t162 * t164;
	t1 = [-t167, 0, 0, 0, 0, 0; -t168, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t156, 0, t159, 0, 0, 0; -t158, 0, -t157, 0, 0, 0; 0, 0, t165, 0, 0, 0; t157, 0, t158, 0, 0, 0; t159, 0, t156, 0, 0, 0; 0, 0, -t166, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:32
	% EndTime: 2019-10-10 00:44:33
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->27), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t235 = cos(qJ(5));
	t237 = cos(qJ(1));
	t259 = t235 * t237;
	t234 = sin(qJ(1));
	t258 = qJD(1) * t234;
	t236 = cos(qJ(3));
	t257 = qJD(1) * t236;
	t256 = qJD(1) * t237;
	t233 = sin(qJ(3));
	t255 = qJD(3) * t233;
	t254 = qJD(3) * t236;
	t253 = qJD(3) * t237;
	t232 = sin(qJ(5));
	t252 = qJD(5) * t232;
	t251 = qJD(5) * t233;
	t250 = qJD(5) * t236;
	t249 = t235 * t255;
	t248 = t235 * t251;
	t247 = t234 * t255;
	t246 = t234 * t254;
	t245 = t233 * t253;
	t244 = t236 * t253;
	t243 = qJD(1) + t250;
	t242 = qJD(5) + t257;
	t241 = t243 * t232;
	t240 = t233 * t256 + t246;
	t239 = t233 * t258 - t244;
	t238 = t242 * t234 + t245;
	t231 = t238 * t232 - t243 * t259;
	t230 = t238 * t235 + t237 * t241;
	t229 = t243 * t235 * t234 + (t242 * t237 - t247) * t232;
	t228 = -t242 * t259 + (t241 + t249) * t234;
	t1 = [t231, 0, t240 * t232 + t234 * t248, 0, t228, 0; -t229, 0, t239 * t232 - t237 * t248, 0, -t230, 0; 0, 0, -t232 * t255 + t235 * t250, 0, -t232 * t251 + t235 * t254, 0; t230, 0, t235 * t246 + (-t234 * t252 + t235 * t256) * t233, 0, t229, 0; t228, 0, -t235 * t244 + (t235 * t258 + t237 * t252) * t233, 0, t231, 0; 0, 0, -t232 * t250 - t249, 0, -t232 * t254 - t248, 0; -t239, 0, t236 * t256 - t247, 0, 0, 0; t240, 0, t234 * t257 + t245, 0, 0, 0; 0, 0, -t254, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:44:33
	% EndTime: 2019-10-10 00:44:33
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (49->27), mult. (173->49), div. (0->0), fcn. (173->6), ass. (0->34)
	t286 = cos(qJ(5));
	t287 = cos(qJ(3));
	t301 = qJD(5) * t287;
	t294 = qJD(1) + t301;
	t311 = t286 * t294;
	t283 = sin(qJ(5));
	t310 = t294 * t283;
	t285 = sin(qJ(1));
	t309 = qJD(1) * t285;
	t308 = qJD(1) * t287;
	t288 = cos(qJ(1));
	t307 = qJD(1) * t288;
	t284 = sin(qJ(3));
	t306 = qJD(3) * t284;
	t305 = qJD(3) * t287;
	t304 = qJD(3) * t288;
	t303 = qJD(5) * t283;
	t302 = qJD(5) * t284;
	t300 = t286 * t306;
	t299 = t286 * t302;
	t298 = t285 * t306;
	t297 = t285 * t305;
	t296 = t284 * t304;
	t295 = t287 * t304;
	t293 = -qJD(5) - t308;
	t292 = t293 * t288;
	t291 = t284 * t307 + t297;
	t290 = t284 * t309 - t295;
	t289 = t293 * t285 - t296;
	t282 = t289 * t283 + t288 * t311;
	t281 = t289 * t286 - t288 * t310;
	t280 = -t285 * t311 + (t292 + t298) * t283;
	t279 = t286 * t292 + (t300 + t310) * t285;
	t1 = [-t282, 0, t291 * t283 + t285 * t299, 0, t279, 0; t280, 0, t290 * t283 - t288 * t299, 0, t281, 0; 0, 0, -t283 * t306 + t286 * t301, 0, -t283 * t302 + t286 * t305, 0; -t290, 0, t287 * t307 - t298, 0, 0, 0; t291, 0, t285 * t308 + t296, 0, 0, 0; 0, 0, -t305, 0, 0, 0; t281, 0, -t286 * t297 + (t285 * t303 - t286 * t307) * t284, 0, t280, 0; -t279, 0, t286 * t295 + (-t286 * t309 - t288 * t303) * t284, 0, t282, 0; 0, 0, t283 * t301 + t300, 0, t283 * t305 + t299, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end