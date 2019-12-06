% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:28
	% EndTime: 2019-12-05 17:32:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:29
	% EndTime: 2019-12-05 17:32:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:29
	% EndTime: 2019-12-05 17:32:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->5)
	t43 = qJD(1) * sin(qJ(1));
	t42 = qJD(1) * cos(qJ(1));
	t39 = cos(pkin(7));
	t38 = sin(pkin(7));
	t1 = [0, 0, 0, 0, 0; t39 * t43, 0, 0, 0, 0; -t39 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t38 * t43, 0, 0, 0, 0; t38 * t42, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0; -t43, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:29
	% EndTime: 2019-12-05 17:32:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t91 = cos(pkin(7));
	t92 = sin(qJ(1));
	t96 = t91 * t92;
	t93 = cos(qJ(1));
	t95 = t91 * t93;
	t94 = qJD(1) * sin(pkin(7));
	t90 = cos(pkin(8));
	t88 = sin(pkin(8));
	t1 = [0, 0, 0, 0, 0; (-t88 * t93 + t90 * t96) * qJD(1), 0, 0, 0, 0; (-t88 * t92 - t90 * t95) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; (-t88 * t96 - t90 * t93) * qJD(1), 0, 0, 0, 0; (t88 * t95 - t90 * t92) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0; t92 * t94, 0, 0, 0, 0; -t93 * t94, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:29
	% EndTime: 2019-12-05 17:32:29
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (46->25), div. (0->0), fcn. (46->8), ass. (0->15)
	t139 = cos(pkin(7));
	t140 = sin(qJ(1));
	t146 = t139 * t140;
	t141 = cos(qJ(1));
	t145 = t139 * t141;
	t144 = qJD(1) * sin(pkin(7));
	t143 = t140 * t144;
	t142 = t141 * t144;
	t138 = cos(pkin(8));
	t137 = cos(pkin(9));
	t135 = sin(pkin(8));
	t134 = sin(pkin(9));
	t133 = (-t135 * t140 - t138 * t145) * qJD(1);
	t132 = (-t135 * t141 + t138 * t146) * qJD(1);
	t1 = [0, 0, 0, 0, 0; t132 * t137 + t134 * t143, 0, 0, 0, 0; t133 * t137 - t134 * t142, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t132 * t134 + t137 * t143, 0, 0, 0, 0; -t133 * t134 - t137 * t142, 0, 0, 0, 0; 0, 0, 0, 0, 0; (t135 * t146 + t138 * t141) * qJD(1), 0, 0, 0, 0; (-t135 * t145 + t138 * t140) * qJD(1), 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:32:30
	% EndTime: 2019-12-05 17:32:30
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (80->27), mult. (266->58), div. (0->0), fcn. (292->10), ass. (0->37)
	t281 = cos(pkin(7));
	t277 = sin(pkin(8));
	t285 = cos(qJ(1));
	t292 = t285 * t277;
	t280 = cos(pkin(8));
	t283 = sin(qJ(1));
	t293 = t283 * t280;
	t287 = t281 * t293 - t292;
	t268 = t287 * qJD(1);
	t276 = sin(pkin(9));
	t279 = cos(pkin(9));
	t278 = sin(pkin(7));
	t290 = qJD(1) * t278;
	t289 = t283 * t290;
	t261 = t268 * t279 + t276 * t289;
	t291 = t285 * t280;
	t294 = t283 * t277;
	t274 = t281 * t291 + t294;
	t295 = t278 * t276;
	t265 = t274 * t279 + t285 * t295;
	t286 = t281 * t294 + t291;
	t267 = t286 * qJD(1);
	t273 = t281 * t292 - t293;
	t282 = sin(qJ(5));
	t284 = cos(qJ(5));
	t304 = (t265 * t284 + t273 * t282) * qJD(5) - t261 * t282 + t267 * t284;
	t301 = t261 * t284 + t267 * t282 + (t265 * t282 - t273 * t284) * qJD(5);
	t296 = t277 * t278;
	t288 = t285 * t290;
	t271 = t278 * t280 * t279 - t281 * t276;
	t270 = t274 * qJD(1);
	t269 = t273 * qJD(1);
	t264 = -t279 * t287 - t283 * t295;
	t263 = -t270 * t279 - t276 * t288;
	t260 = t263 * t284 - t269 * t282 + (-t264 * t282 - t284 * t286) * qJD(5);
	t259 = -t263 * t282 - t269 * t284 + (-t264 * t284 + t282 * t286) * qJD(5);
	t1 = [0, 0, 0, 0, (-t271 * t284 - t282 * t296) * qJD(5); t301, 0, 0, 0, t259; t260, 0, 0, 0, -t304; 0, 0, 0, 0, (t271 * t282 - t284 * t296) * qJD(5); t304, 0, 0, 0, -t260; t259, 0, 0, 0, t301; 0, 0, 0, 0, 0; t268 * t276 - t279 * t289, 0, 0, 0, 0; -t270 * t276 + t279 * t288, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end