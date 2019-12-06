% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t20 = qJD(2) * sin(qJ(2));
	t19 = qJD(2) * cos(qJ(2));
	t16 = cos(pkin(7));
	t15 = sin(pkin(7));
	t1 = [0, -t16 * t19, 0, 0, 0; 0, -t15 * t19, 0, 0, 0; 0, -t20, 0, 0, 0; 0, t16 * t20, 0, 0, 0; 0, t15 * t20, 0, 0, 0; 0, -t19, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:25
	% EndTime: 2019-12-05 16:14:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->15), mult. (77->36), div. (0->0), fcn. (77->6), ass. (0->17)
	t154 = sin(qJ(3));
	t157 = cos(qJ(2));
	t167 = t154 * t157;
	t156 = cos(qJ(3));
	t166 = t156 * t157;
	t155 = sin(qJ(2));
	t165 = qJD(2) * t155;
	t164 = qJD(2) * t157;
	t163 = qJD(3) * t155;
	t162 = qJD(3) * t157;
	t152 = sin(pkin(7));
	t161 = t152 * t165;
	t153 = cos(pkin(7));
	t160 = t153 * t165;
	t159 = t154 * t163 - t156 * t164;
	t158 = t154 * t164 + t156 * t163;
	t1 = [0, t159 * t153, t154 * t160 + (-t152 * t154 - t153 * t166) * qJD(3), 0, 0; 0, t159 * t152, t154 * t161 + (-t152 * t166 + t153 * t154) * qJD(3), 0, 0; 0, -t154 * t162 - t156 * t165, -t158, 0, 0; 0, t158 * t153, t156 * t160 + (-t152 * t156 + t153 * t167) * qJD(3), 0, 0; 0, t158 * t152, t156 * t161 + (t152 * t167 + t153 * t156) * qJD(3), 0, 0; 0, t154 * t165 - t156 * t162, t159, 0, 0; 0, -t160, 0, 0, 0; 0, -t161, 0, 0, 0; 0, t164, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:25
	% EndTime: 2019-12-05 16:14:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (32->24), mult. (144->59), div. (0->0), fcn. (144->8), ass. (0->25)
	t241 = sin(qJ(3));
	t244 = cos(qJ(2));
	t258 = t241 * t244;
	t242 = sin(qJ(2));
	t243 = cos(qJ(3));
	t257 = t242 * t243;
	t256 = t243 * t244;
	t255 = qJD(2) * t242;
	t254 = qJD(2) * t244;
	t253 = qJD(3) * t242;
	t252 = qJD(3) * t244;
	t251 = t241 * t255;
	t250 = t243 * t255;
	t249 = t241 * t252;
	t248 = t241 * t253;
	t247 = t241 * t254 + t243 * t253;
	t237 = sin(pkin(8));
	t239 = cos(pkin(8));
	t246 = -t237 * t248 + (t237 * t256 - t239 * t242) * qJD(2);
	t245 = t239 * t248 + (-t237 * t242 - t239 * t256) * qJD(2);
	t240 = cos(pkin(7));
	t238 = sin(pkin(7));
	t236 = t240 * t251 + (-t238 * t241 - t240 * t256) * qJD(3);
	t235 = t238 * t251 + (-t238 * t256 + t240 * t241) * qJD(3);
	t1 = [0, t245 * t240, t236 * t239, 0, 0; 0, t245 * t238, t235 * t239, 0, 0; 0, -t239 * t249 + (t237 * t244 - t239 * t257) * qJD(2), -t247 * t239, 0, 0; 0, t246 * t240, -t236 * t237, 0, 0; 0, t246 * t238, -t235 * t237, 0, 0; 0, t237 * t249 + (t237 * t257 + t239 * t244) * qJD(2), t247 * t237, 0, 0; 0, -t247 * t240, -t240 * t250 + (t238 * t243 - t240 * t258) * qJD(3), 0, 0; 0, -t247 * t238, -t238 * t250 + (-t238 * t258 - t240 * t243) * qJD(3), 0, 0; 0, t243 * t252 - t251, t243 * t254 - t248, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:25
	% EndTime: 2019-12-05 16:14:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (30->19), mult. (144->59), div. (0->0), fcn. (144->8), ass. (0->25)
	t281 = sin(qJ(3));
	t284 = cos(qJ(2));
	t298 = t281 * t284;
	t282 = sin(qJ(2));
	t283 = cos(qJ(3));
	t297 = t282 * t283;
	t296 = t283 * t284;
	t295 = qJD(2) * t282;
	t294 = qJD(2) * t284;
	t293 = qJD(3) * t282;
	t292 = qJD(3) * t284;
	t291 = t281 * t295;
	t290 = t283 * t295;
	t289 = t281 * t292;
	t288 = t281 * t293;
	t287 = -t281 * t294 - t283 * t293;
	t277 = sin(pkin(8));
	t279 = cos(pkin(8));
	t286 = t277 * t288 + (-t277 * t296 + t279 * t282) * qJD(2);
	t285 = t279 * t288 + (-t277 * t282 - t279 * t296) * qJD(2);
	t280 = cos(pkin(7));
	t278 = sin(pkin(7));
	t276 = t280 * t291 + (-t278 * t281 - t280 * t296) * qJD(3);
	t275 = t278 * t291 + (-t278 * t296 + t280 * t281) * qJD(3);
	t1 = [0, t285 * t280, t276 * t279, 0, 0; 0, t285 * t278, t275 * t279, 0, 0; 0, -t279 * t289 + (t277 * t284 - t279 * t297) * qJD(2), t287 * t279, 0, 0; 0, t287 * t280, -t280 * t290 + (t278 * t283 - t280 * t298) * qJD(3), 0, 0; 0, t287 * t278, -t278 * t290 + (-t278 * t298 - t280 * t283) * qJD(3), 0, 0; 0, t283 * t292 - t291, t283 * t294 - t288, 0, 0; 0, t286 * t280, t276 * t277, 0, 0; 0, t286 * t278, t275 * t277, 0, 0; 0, -t277 * t289 + (-t277 * t297 - t279 * t284) * qJD(2), t287 * t277, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end