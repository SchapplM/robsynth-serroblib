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
% Datum: 2019-10-24 10:39
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-24 10:39:20
	% EndTime: 2019-10-24 10:39:20
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (80->27), mult. (266->58), div. (0->0), fcn. (292->10), ass. (0->37)
	t282 = cos(pkin(7));
	t278 = sin(pkin(8));
	t286 = cos(qJ(1));
	t292 = t286 * t278;
	t281 = cos(pkin(8));
	t284 = sin(qJ(1));
	t293 = t284 * t281;
	t288 = t282 * t293 - t292;
	t269 = t288 * qJD(1);
	t277 = sin(pkin(9));
	t280 = cos(pkin(9));
	t279 = sin(pkin(7));
	t291 = qJD(1) * t279;
	t290 = t284 * t291;
	t262 = t269 * t280 + t277 * t290;
	t294 = t284 * t278;
	t295 = t281 * t286;
	t275 = t282 * t295 + t294;
	t297 = t277 * t279;
	t266 = t275 * t280 + t286 * t297;
	t287 = t282 * t294 + t295;
	t268 = t287 * qJD(1);
	t274 = t282 * t292 - t293;
	t283 = sin(qJ(5));
	t285 = cos(qJ(5));
	t305 = (t266 * t285 + t274 * t283) * qJD(5) - t262 * t283 + t268 * t285;
	t302 = t262 * t285 + t268 * t283 + (t266 * t283 - t274 * t285) * qJD(5);
	t296 = t278 * t279;
	t289 = t286 * t291;
	t272 = t279 * t280 * t281 - t277 * t282;
	t271 = t275 * qJD(1);
	t270 = t274 * qJD(1);
	t265 = -t288 * t280 - t284 * t297;
	t264 = -t271 * t280 - t277 * t289;
	t261 = t264 * t285 - t270 * t283 + (-t265 * t283 - t285 * t287) * qJD(5);
	t260 = -t264 * t283 - t270 * t285 + (-t265 * t285 + t283 * t287) * qJD(5);
	t1 = [0, 0, 0, 0, (-t272 * t285 - t283 * t296) * qJD(5); t302, 0, 0, 0, t260; t261, 0, 0, 0, -t305; 0, 0, 0, 0, (t272 * t283 - t285 * t296) * qJD(5); t305, 0, 0, 0, -t261; t260, 0, 0, 0, t302; 0, 0, 0, 0, 0; t269 * t277 - t280 * t290, 0, 0, 0, 0; -t271 * t277 + t280 * t289, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end